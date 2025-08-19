/**
 * @file bezier_base.h
 * @brief Bezier曲线基础类的头文件
 * 
 * 该头文件为基于Bernstein基函数的轨迹生成优化问题提供基础数学支持，包括：
 * 1. 映射矩阵：将Bernstein基函数系数（即控制点）映射到单项式基函数，支持3-12阶
 * 2. 模数列表：预计算给定阶数的Bernstein基函数的常数模数（组合数'n choose k'）
 *    以节省频繁调用这些数值的计算成本
 * 
 * 该类应在轨迹生成器调用前初始化为实例。
 * 提供了多个初始化器，实例根据给定的控制点阶数进行初始化。
 * 
 * @author UAV轨迹规划项目组
 * @date 2023
 */

#ifndef _BEZIER_BASE_H_
#define _BEZIER_BASE_H_

#include <stdio.h>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;

/**
 * @class Bernstein
 * @brief Bernstein基函数类，用于Bezier曲线轨迹生成
 * 
 * 该类实现了基于Bernstein基函数的数学运算支持，主要用于UAV轨迹规划中的
 * Bezier曲线生成和优化。提供映射矩阵、成本矩阵和组合系数的预计算功能。
 */
class Bernstein
{
	private:
		
		vector<MatrixXd> MQMList, MList;  ///< 映射矩阵列表：M为映射矩阵，MQM为映射后的成本矩阵
		vector<MatrixXd> FMList;          ///< Cholesky分解后的矩阵列表：F*M
		vector<VectorXd> CList, CvList, CaList, CjList;  ///< Bernstein系数列表：位置、速度、加速度、jerk

		int _order_min, _order_max;   ///< 每段多项式的最小和最大阶数，也是每段使用的控制点数量
		double _min_order;            ///< 最小化的导数阶数：1-速度，2-加速度，3-jerk，4-snap    

	public:
		/**
		 * @brief 默认构造函数
		 */
		Bernstein(){} // Empty constructor

		/**
		 * @brief 参数化构造函数
		 * @param poly_order_min 多项式最小阶数
		 * @param poly_order_max 多项式最大阶数
		 * @param min_order 最小化的导数阶数
		 */
		Bernstein(int poly_order_min, int poly_order_max, double min_order)
		{
			_order_min = poly_order_min;
			_order_max = poly_order_max;
			_min_order = min_order; 
		}

		/**
		 * @brief 析构函数
		 */
		~Bernstein(){}

		/**
		 * @brief 设置Bernstein基函数参数并预计算相关矩阵
		 * @param poly_order_min 多项式最小阶数
		 * @param poly_order_max 多项式最大阶数
		 * @param min_order 最小化的导数阶数
		 * @return 1表示成功，-1表示参数错误
		 */
		int setParam(int poly_order_min, int poly_order_max, double min_order);
		
		/**
		 * @brief 对矩阵Q进行Cholesky分解
		 * @param Q 输入的正定矩阵
		 * @return Q的平方根矩阵F，使得Q = F' * F
		 */
		MatrixXd CholeskyDecomp(MatrixXd Q); // return square root F of Q; Q = F' * F

		// Getter函数组
		/**
		 * @brief 获取映射矩阵列表
		 * @return 映射矩阵M的向量列表
		 */
		vector<MatrixXd> getM(){ return MList; }
		
		/**
		 * @brief 获取映射后成本矩阵列表
		 * @return MQM矩阵（M^T * Q * M）的向量列表
		 */
		vector<MatrixXd> getMQM(){ return MQMList; }
		
		/**
		 * @brief 获取Cholesky分解后的矩阵列表
		 * @return FM矩阵（F * M）的向量列表
		 */
		vector<MatrixXd> getFM(){ return FMList; }
		
		/**
		 * @brief 获取位置的Bernstein系数列表
		 * @return 位置系数向量列表
		 */
		vector<VectorXd> getC(){ return CList; }
		
		/**
		 * @brief 获取速度的Bernstein系数列表
		 * @return 速度系数向量列表
		 */
		vector<VectorXd> getC_v(){ return CvList; }
		
		/**
		 * @brief 获取加速度的Bernstein系数列表
		 * @return 加速度系数向量列表
		 */
		vector<VectorXd> getC_a(){ return CaList; }
		
		/**
		 * @brief 获取jerk的Bernstein系数列表
		 * @return jerk系数向量列表
		 */
		vector<VectorXd> getC_j(){ return CjList; }
};

#endif