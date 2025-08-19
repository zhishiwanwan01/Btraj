/**
 * @file bezier_base.cpp
 * @brief Bezier曲线基础类的实现文件
 * 
 * 该文件实现了Bezier曲线的基础功能，包括：
 * - Cholesky分解用于优化计算
 * - 参数设置和预计算
 * - 不同阶数的映射矩阵生成
 * - 组合数计算和系数向量生成
 * 
 * 支持3-12阶的Bezier曲线，用于UAV轨迹规划
 */

#include "bezier_base.h"

/**
 * @brief 对矩阵Q进行Cholesky分解
 * @param Q 输入的正定矩阵
 * @return 返回Q的平方根矩阵F，使得Q = F' * F
 * 
 * 使用LDLT分解来计算Cholesky分解，避免数值不稳定问题
 */
MatrixXd Bernstein::CholeskyDecomp(MatrixXd Q) // return square root F of Q; Q = F' * F
{
	MatrixXd F, Ft;
	Eigen::LDLT< MatrixXd > ldlt(Q);
    F = ldlt.matrixL();
    F = ldlt.transpositionsP().transpose() * F;
    F *= ldlt.vectorD().array().sqrt().matrix().asDiagonal();
	Ft = F.transpose();

	return Ft;
}

/**
 * @brief 设置Bezier曲线参数并预计算相关矩阵
 * @param poly_order_min 多项式最小阶数
 * @param poly_order_max 多项式最大阶数  
 * @param min_order 最小导数阶数（用于成本函数）
 * @return 1表示成功，-1表示参数错误
 * 
 * 该函数预计算所有需要的矩阵和系数：
 * - 映射矩阵M：将Bezier系数映射到单项式多项式
 * - 成本矩阵Q：用于优化的Hessian矩阵
 * - 组合系数：位置、速度、加速度、jerk的Bernstein系数
 */
int Bernstein::setParam(int poly_order_min, int poly_order_max, double min_order)
{
	int ret = (poly_order_min >=3 && poly_order_max <= 13) ? 1 : -1; // Now we only applied for poly order from 3 ~ 12 ( num of control from 4 ~ 13 )

	_order_min = poly_order_min;
	_order_max = poly_order_max;
	_min_order = min_order; 

	/**
	 * @brief 计算阶乘的lambda函数
	 * @param n 输入整数
	 * @return n的阶乘
	 */
	const static auto factorial = [](int n)
	{
        int fact = 1;

        for(int i = n; i > 0 ; i--)
          fact *= i;

        return fact; 
	};

	/**
	 * @brief 计算组合数的lambda函数
	 * @param n 总数
	 * @param k 选择数
	 * @return C(n,k) = n!/(k!(n-k)!)，即从n个元素中选k个的组合数
	 */
	const static auto combinatorial = [](int n, int k) // for calculate n choose k combination problem
	{
		return factorial(n) / (factorial(k) * factorial(n - k));
	};

	CList.clear();
	CvList.clear();
	CaList.clear();
	CjList.clear();
	
	// 为每个多项式阶数预计算相关矩阵和系数
	for(int order = 0; order <= _order_max; order++)
	{	
		MatrixXd M;   // 映射矩阵，用于将Bezier曲线系数映射到单项式多项式

		MatrixXd Q, Q_l, Q_u; // 目标函数中每个块的成本Hessian矩阵，无缩放，仅元素
		MatrixXd MQM; 	      // 目标函数中每个块的M' * Q * M，无缩放，仅元素
		
		VectorXd C;   // 位置系数向量，记录Bernstein系数的预计算组合数'n choose k'
		VectorXd C_v; // 速度系数向量
		VectorXd C_a; // 加速度系数向量
		VectorXd C_j; // jerk系数向量

		int poly_num1D = order + 1; 
		M.resize(order + 1, order + 1);
		Q_l = MatrixXd::Zero(order + 1, order + 1);
		Q_u = MatrixXd::Zero(order + 1, order + 1);
		Q   = MatrixXd::Zero(order + 1, order + 1);
		
		C.resize  (order + 1);
		C_v.resize(order    );
		
		if(order > 1)
			C_a.resize(order - 1);

		if(order > 2)
			C_j.resize(order - 2);

		int min_order_l = floor(_min_order);
		int min_order_u = ceil (_min_order);

		// 计算下界成本矩阵Q_l（基于最小导数阶数的下界）
		if( poly_num1D > min_order_l )
		{
			for( int i = min_order_l; i < poly_num1D; i ++ )
			{
		    	for(int j = min_order_l; j < poly_num1D; j ++ )
		    	{
		      		int coeff = 1.0;
		            int _d = min_order_l - 1;
		                
		            while(_d >= 0)
		            {
		                coeff = coeff * (i - _d) * (j - _d);
		                _d -= 1;
		            }

		            Q_l(i,j) = double(coeff) / double(i + j - 2 * min_order_l + 1);
		    	}
		  	}
		}

		// 计算上界成本矩阵Q_u（基于最小导数阶数的上界）
		if( poly_num1D > min_order_u )
		{
			for( int i = min_order_u; i < poly_num1D; i ++ )
			{
		    	for(int j = min_order_u; j < poly_num1D; j ++ )
		    	{
		      		int coeff = 1.0;
		            int _d = min_order_u - 1;
		                
		            while(_d >= 0)
		            {
		                coeff = coeff * (i - _d) * (j - _d);
		                _d -= 1;
		            }

		            Q_u(i,j) = double(coeff) / double(i + j - 2 * min_order_u + 1);
		    	}
		  	}
		}

		// 根据min_order的整数部分选择或插值成本矩阵
		if(min_order_l == min_order_u)
			Q = Q_u;
		else
			Q = (_min_order - min_order_l) * Q_u + (min_order_u - _min_order) * Q_l;

		/**
		 * @brief 为不同阶数设置映射矩阵M
		 * 
		 * 映射矩阵将Bezier控制点系数映射到单项式多项式系数
		 * 每个case对应一个特定阶数的Bezier曲线
		 * 矩阵元素通过Bezier基函数的展开式计算得出
		 */
	 	switch(order)
		{	
			case 0: 
			{
				M << 1;

				break;

			}
			case 1: 
			{
				M << -1,  0,
					 -1,  1;
				break;

			}
			case 2:
			{
				M << -1,  0,  0,
					 -2,  2,  0,
					  1, -2,  1;
				break;

			}
			case 3: 
			{
				M << -1,  0,  0,  0,
					 -3,  3,  0,  0,
					  3, -6,  3,  0,
					 -1,  3, -3,  1;	
				break;

			}
			case 4:
			{
				M <<  1,   0,   0,   0,  0,
					 -4,   4,   0,   0,  0,
					  6, -12,   6,   0,  0,
					 -4,  12, -12,   4,  0,
					  1,  -4,   6,  -4,  1;
				break;
			}
			case 5:
			{
				M << 1,   0,   0,   0,  0,  0,
					-5,   5,   0,   0,  0,  0,
					10, -20,  10,   0,  0,  0,
				   -10,  30, -30,  10,  0,  0,
				     5, -20,  30, -20,  5,  0,
				    -1,   5, -10,  10, -5,  1;
				break;
			}
			case 6:
			{	

				M << 1,   0,   0,   0,   0,  0,  0,
					-6,   6,   0,   0,   0,  0,  0,
					15, -30,  15,   0,   0,  0,  0,
				   -20,  60, -60,  20,   0,  0,  0,
				    15, -60,  90, -60,  15,  0,  0,
				    -6,  30, -60,  60, -30,  6,  0,
				     1,  -6,  15, -20,  15, -6,  1;
				break;
			}
			case 7:
			{
				M << 1,    0,    0,    0,    0,   0,   0,   0,
				    -7,    7,    0,    0,    0,   0,   0,   0,
				    21,   42,   21,    0,    0,   0,   0,   0,
				   -35,  105, -105,   35,    0,   0,   0,   0, 
				    35, -140,  210, -140,   35,   0,   0,   0,
				   -21,  105, -210,  210, -105,  21,   0,   0,
				     7,  -42,  105, -140,  105, -42,   7,   0,
				    -1,    7,  -21,   35,  -35,  21,  -7,   1;
				break;
			}
			case 8:
			{
				M << 1,    0,    0,    0,    0,    0,   0,   0,   0,
				    -8,    8,    0,    0,    0,    0,   0,   0,   0,
				    28,  -56,   28,    0,    0,    0,   0,   0,   0,
				   -56,  168, -168,   56,    0,    0,   0,   0,   0, 
				    70, -280,  420, -280,   70,    0,   0,   0,   0,
				   -56,  280, -560,  560, -280,   56,   0,   0,   0,
				    28, -168,  420, -560,  420, -168,  28,   0,   0,
				    -8,   56, -168,  280, -280,  168, -56,   8,   0,
				     1,   -8,   28,  -56,   70,  -56,  28,  -8,   1;
				break;
			}
			case 9:
			{
				M << 1,    0,     0,     0,     0,    0,    0,     0,     0,    0,
				    -9,    9,     0,     0,     0,    0,    0,     0,     0,    0, 
				    36,  -72,    36,     0,     0,    0,    0,     0,     0,    0, 
				   -84,  252,  -252,    84,     0,    0,    0,     0,     0,    0, 
				   126, -504,   756,  -504,   126,    0,    0,     0,     0,    0,
				  -126,  630, -1260,  1260,  -630,  126,    0,     0,     0,    0,
				    84, -504,  1260, -1680,  1260, -504,   84,     0,     0,    0,
				   -36,  252,  -756,  1260, -1260,  756, -252,    36,     0,    0,
				     9,  -72,   252,  -504,   630, -504,  252,   -72,     9,    0,
				    -1,    9,   -36,    84,  -126,  126,  -84,    36,    -9,    1;
				break;
			}
			case 10:
			{
				M <<  1,     0,     0,     0,      0,     0,    0,     0,     0,    0,   0,
				    -10,    10,     0,     0,      0,     0,    0,     0,     0,    0,   0,
				     45,   -90,    45,     0,      0,     0,    0,     0,     0,    0,   0,
				   -120,   360,  -360,   120,      0,     0,    0,     0,     0,    0,   0,
				    210,  -840,  1260,  -840,    210,     0,    0,     0,     0,    0,   0,
				   -252,  1260, -2520,  2520,  -1260,   252,    0,     0,     0,    0,   0,
				    210, -1260,  3150, -4200,   3150, -1260,  210,     0,     0,    0,   0,
				   -120,  840,  -2520,  4200,  -4200,  2520, -840,   120,     0,    0,   0,
				     45, -360,   1260, -2520,   3150, -2520, 1260,  -360,    45,    0,   0,
				    -10,   90,   -360,   840,  -1260,  1260, -840,   360,   -90,   10,   0,
				      1,  -10,     45,  -120,    210,  -252,  210,  -120,    45,  -10,   1;
				break;
			}
			case 11:
			{
				M <<  1,     0,    0,      0,      0,      0,     0,     0,     0,    0,   0,  0,
				    -11,    11,    0,      0,      0,      0,     0,     0,     0,    0,   0,  0,
				     55,  -110,   55,      0,      0,      0,     0,     0,     0,    0,   0,  0,
				   -165,   495, -495,    165,      0,      0,     0,     0,     0,    0,   0,  0,
				    330, -1320, 1980,  -1320,    330,      0,     0,     0,     0,    0,   0,  0,
				   -462,  2310, -4620,  4620,  -2310,    462,     0,     0,     0,    0,   0,  0,
				    462, -2772,  6930, -9240,   6930,  -2772,   462,     0,     0,    0,   0,  0,
				   -330,  2310, -6930, 11550, -11550,   6930, -2310,   330,     0,    0,   0,  0,
				    165, -1320,  4620, -9240,  11550,  -9240,  4620, -1320,   165,    0,   0,  0,
				    -55,   495, -1980,  4620,  -6930,   6930, -4620,  1980,  -495,   55,   0,  0,
				     11,  -110,   495, -1320,   2310,  -2772,  2310, -1320,   495, -110,  11,  0,
				     -1,    11,   -55,   165,   -330,    462,  -462,   330,  -165,   55, -11,  1;
				break;
			}
			case 12:
			{
				M <<  1,     0,     0,      0,      0,      0,     0,     0,     0,    0,    0,   0,   0,
				    -12,    12,     0,      0,      0,      0,     0,     0,     0,    0,    0,   0,   0,
				     66,  -132,    66,      0,      0,      0,     0,     0,     0,    0,    0,   0,   0,
				   -220,   660,  -660,    220,      0,      0,     0,     0,     0,    0,    0,   0,   0,
				    495, -1980,  2970,  -1980,    495,      0,     0,     0,     0,    0,    0,   0,   0, 
				   -792,  3960, -7920,   7920,  -3960,    792,     0,     0,     0,    0,    0,   0,   0,
				    924, -5544, 13860, -18480,  13860,  -5544,   924,     0,     0,    0,    0,   0,   0,
				   -792,  5544,-16632,  27720, -27720,  16632, -5544,   792,     0,    0,    0,   0,   0,
				    495, -3960, 13860, -27720,  34650, -27720, 13860, -3960,   495,    0,    0,   0,   0,
				   -220,  1980, -7920,  18480, -27720,  27720,-18480,  7920, -1980,  220,    0,   0,   0,
				     66,  -660,  2970,  -7920,  13860, -16632, 13860, -7920,  2970, -660,   66,   0,   0,
				    -12,   132,  -660,   1980,  -3960,   5544, -5544,  3960, -1980,  660, -132,  12,   0,
				      1,   -12,    66,   -220,    495,   -792,   924,  -792,   495, -220,   66, -12,   1;
				break;
			}
		}
		
		// 存储映射矩阵
		MList.push_back(M);
		// 计算映射后的成本矩阵块
		MQM = M.transpose() * Q * M; // Get the cost block after mapping the coefficients

		// 计算Cholesky分解和相关矩阵
		MatrixXd F  = CholeskyDecomp(Q);
		MatrixXd FM = F * M;

		MQMList.push_back(MQM);
		FMList.push_back(FM);

		// 计算Bernstein系数的组合数
		int n = order;
		for(int k = 0; k <= n; k ++ )
		{
			C(k)   = combinatorial(n, k);     // 位置的Bernstein系数
			
			if( k <= (n - 1) )
				C_v(k) = combinatorial(n - 1, k);  // 速度的Bernstein系数
			if( k <= (n - 2) )
				C_a(k) = combinatorial(n - 2, k);  // 加速度的Bernstein系数
			if( k <= (n - 3) )
				C_j(k) = combinatorial(n - 3, k);  // jerk的Bernstein系数
		}

		// 存储所有系数向量
		CList.push_back(C);
		CvList.push_back(C_v);
		CaList.push_back(C_a);
		CjList.push_back(C_j);

	}
	return ret;
};