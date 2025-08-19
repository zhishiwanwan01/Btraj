/*
 * 基于MOSEK二次规划的贝塞尔轨迹生成器
 * 
 * 本文件实现了Btraj规划器的核心轨迹生成算法。
 * 通过飞行立方体走廊生成平滑的贝塞尔曲线轨迹，
 * 同时强制执行速度、加速度和连续性约束。
 * 
 * 主函数：BezierPloyCoeffGeneration
 * - 建立带有空间走廊约束的QP问题
 * - 在轨迹段连接点强制连续性（位置、速度、加速度）
 * - 可选的速度和加速度幅值限制
 * - 最小化轨迹导数（速度、加速度、急动度或snap）
 * - 使用MOSEK内点法求解器进行优化
 */

#include "trajectory_generator.h"
using namespace std;    
using namespace Eigen;

/*
 * MOSEK求解器回调函数
 * 
 * 该回调函数由MOSEK在优化过程中用于输出求解器
 * 进度信息和调试消息。
 */
static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
}

/*
 * ===================================================================
 * 主要贝塞尔轨迹生成函数
 * ===================================================================
 * 
 * 该函数建立并求解二次规划问题，通过空间约束走廊
 * 生成平滑的贝塞尔曲线轨迹。
 * 
 * 参数：
 * - corridor: 由立方体约束定义的飞行走廊
 * - MQM: 针对选定导数阶数的预计算最小化矩阵
 * - pos/vel/acc: 起点和终点的边界条件
 * - maxVel/maxAcc: 动态约束限制
 * - traj_order: 贝塞尔曲线阶数（通常为5-10）
 * - minimize_order: 要最小化的导数（1=速度, 2=加速度, 3=急动度, 4=snap）
 * - margin: 走廊立方体内的安全裕度
 * - isLimitVel/isLimitAcc: 是否强制执行速度/加速度约束
 * 
 * 返回值：通过obj返回优化目标值，通过PolyCoeff返回系数矩阵
 */
int TrajectoryGenerator::BezierPloyCoeffGeneration(
            const vector<Cube> &corridor,
            const MatrixXd &MQM,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,
            const double minimize_order,
            const double margin,
            const bool & isLimitVel,
            const bool & isLimitAcc,
            double & obj,
            MatrixXd & PolyCoeff)  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
{   
/*
 * ===================================================================
 * 问题规模和参数计算
 * ===================================================================
 */
#define ENFORCE_VEL  isLimitVel // 是否添加额外约束以确保速度可行性
#define ENFORCE_ACC  isLimitAcc // 是否添加额外约束以确保加速度可行性

    // 从走廊中提取时间缩放因子
    double initScale = corridor.front().t;
    double lstScale  = corridor.back().t;
    int segment_num  = corridor.size();

    // 控制点维度
    int n_poly = traj_order + 1;                    // 每段每维度的控制点数量
    int s1d1CtrlP_num = n_poly;                     // 一段一维的控制点数量
    int s1CtrlP_num   = 3 * s1d1CtrlP_num;          // 一段三维的控制点数量

    // 等式约束数量：边界条件 + 连续性
    int equ_con_s_num = 3 * 3; // 起点在x, y, z轴的p, v, a
    int equ_con_e_num = 3 * 3; // 终点在x, y, z轴的p, v, a
    int equ_con_continuity_num = 3 * 3 * (segment_num - 1);
    int equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // 每个段连接点在x, y, z轴的p, v, a
    
    // 不等式约束数量：速度和加速度限制
    int vel_con_num = 3 *  traj_order * segment_num;
    int acc_con_num = 3 * (traj_order - 1) * segment_num;

    if( !ENFORCE_VEL )
        vel_con_num = 0;

    if( !ENFORCE_ACC )
        acc_con_num = 0;

    int high_order_con_num = vel_con_num + acc_con_num; 
    //int high_order_con_num = 0; //3 * traj_order * segment_num;

    // 总问题维度
    int con_num   = equ_con_num + high_order_con_num;   // 总约束数量
    int ctrlP_num = segment_num * s1CtrlP_num;          // 总控制点数量

    double x_var[ctrlP_num];
    double primalobj;

    MSKrescodee  r; 
    vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
    
/*
 * ===================================================================
 * 约束边界设置
 * ===================================================================
 */
    
    // 速度幅值约束：-maxVel <= vel <= +maxVel 
    if(ENFORCE_VEL)
    {
        for(int i = 0; i < vel_con_num; i++)
        {
            pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxVel,  + maxVel) );
            con_bdk.push_back(cb_ie);   
        }
    }

    // 加速度幅值约束：-maxAcc <= acc <= +maxAcc
    if(ENFORCE_ACC)
    {
        for(int i = 0; i < acc_con_num; i++)
        {
            pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxAcc,  maxAcc) ); 
            con_bdk.push_back(cb_ie);   
        }
    }

    // 等式约束：起点/终点边界条件和连续性
    // 顺序：起点位置(3), 起点速度(3), 起点加速度(3), 终点位置(3), 终点速度(3), 终点加速度(3), 连续性(...)
    for(int i = 0; i < equ_con_num; i ++ ){ 
        double beq_i;
        if(i < 3)                    beq_i = pos(0, i);        // 起点位置 (x,y,z)
        else if (i >= 3  && i < 6  ) beq_i = vel(0, i - 3);    // 起点速度 (x,y,z)
        else if (i >= 6  && i < 9  ) beq_i = acc(0, i - 6);    // 起点加速度 (x,y,z)
        else if (i >= 9  && i < 12 ) beq_i = pos(1, i - 9 );   // 终点位置 (x,y,z)
        else if (i >= 12 && i < 15 ) beq_i = vel(1, i - 12);   // 终点速度 (x,y,z)
        else if (i >= 15 && i < 18 ) beq_i = acc(1, i - 15);   // 终点加速度 (x,y,z)
        else beq_i = 0.0;                                       // 连续性约束 (= 0)

        pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); 
        con_bdk.push_back(cb_eq);
    }

/*
 * ===================================================================
 * 变量边界设置 - 空间走廊约束
 * ===================================================================
 */
    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 

    // 对于每个走廊段，将控制点约束在空间边界内
    for(int k = 0; k < segment_num; k++)
    {   
        Cube cube_     = corridor[k];
        double scale_k = cube_.t;

        // 对于每个维度 (x, y, z)
        for(int i = 0; i < 3; i++ )
        {   
            // 对于该段该维度的每个控制点
            for(int j = 0; j < n_poly; j ++ )
            {   
                pair<MSKboundkeye, pair<double, double> > vb_x;

                double lo_bound, up_bound;
                // 除了第一段（为了精确保持起始位置）外，其他段都应用安全裕度
                if(k > 0)
                {
                    lo_bound = (cube_.box[i].first  + margin) / scale_k;
                    up_bound = (cube_.box[i].second - margin) / scale_k;
                }
                else
                {
                    lo_bound = (cube_.box[i].first)  / scale_k;
                    up_bound = (cube_.box[i].second) / scale_k;
                }

                vb_x  = make_pair( MSK_BK_RA, make_pair( lo_bound, up_bound ) ); 

                var_bdk.push_back(vb_x);
            }
        } 
    }

/*
 * ===================================================================
 * MOSEK优化问题设置
 * ===================================================================
 */
    MSKint32t  j,i; 
    MSKenv_t   env; 
    MSKtask_t  task; 
    
    // 创建MOSEK环境和优化任务
    r = MSK_makeenv( &env, NULL ); 
    r = MSK_maketask(env,con_num, ctrlP_num, &task); 

    // 配置MOSEK求解器参数以获得性能和精度
    MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);                           // 单线程执行
    MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-2);           // 凸性检查容差
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS,  1e-4);                // 对偶可行性容差
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS,  1e-4);                // 原始可行性容差
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-4);                // 不可行性容差
    //MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 5e-2 );            // 相对间隙容差（已注释）
    
    //r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);        // 链接调试输出回调（已注释）
    
    // 向优化问题添加约束和变量
    if ( r == MSK_RES_OK ) 
      r = MSK_appendcons(task,con_num);     // 添加约束槽
    if ( r == MSK_RES_OK ) 
      r = MSK_appendvars(task,ctrlP_num);   // 添加变量槽（初始为零）

    // 设置变量边界（控制点空间约束）
    for(j = 0; j<ctrlP_num && r == MSK_RES_OK; ++j){ 
        if (r == MSK_RES_OK) 
            r = MSK_putvarbound(task, 
                                j,                            // 变量索引
                                var_bdk[j].first,             // 边界键
                                var_bdk[j].second.first,      // 下界数值
                                var_bdk[j].second.second );   // 上界数值
    } 
    
    // 设置约束边界（速度、加速度和等式约束）
    for( i = 0; i < con_num && r == MSK_RES_OK; i++ ) {
        r = MSK_putconbound(task, 
                            i,                            // 约束索引
                            con_bdk[i].first,             // 边界键
                            con_bdk[i].second.first,      // 下界数值
                            con_bdk[i].second.second );   // 上界数值
    }

/*
 * ===================================================================
 * 约束矩阵构建（线性部分A）
 * ===================================================================
 * 
 * 本节为QP问题Ax = b构建约束矩阵A。
 * 约束按以下顺序逐行添加：
 * 1. 速度幅值约束（如果启用）
 * 2. 加速度幅值约束（如果启用）
 * 3. 起点边界条件（位置、速度、加速度）
 * 4. 终点边界条件（位置、速度、加速度）
 * 5. 段连接点的连续性约束（位置、速度、加速度）
 */
    int row_idx = 0;
    
    // 1. 速度幅值约束：|velocity| <= maxVel
    // 使用有限差分：velocity = n * (P_{i+1} - P_i)，其中n为轨迹阶数
    if(ENFORCE_VEL)
    {   
        for(int k = 0; k < segment_num ; k ++ )
        {   
            for(int i = 0; i < 3; i++)              // 对于每个维度 (x, y, z)
            {  
                for(int p = 0; p < traj_order; p++) // 对于每个速度约束点
                {
                    int nzi = 2;                    // 每个约束2个非零元素
                    MSKint32t asub[nzi];
                    double aval[nzi];

                    // 速度有限差分的系数：n * (P_{i+1} - P_i)
                    aval[0] = -1.0 * traj_order;    // -n * P_i
                    aval[1] =  1.0 * traj_order;    // +n * P_{i+1}

                    // 控制点P_i和P_{i+1}的变量索引
                    asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    
                    asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;    

                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                }
            }
        }
    }

    // 2. 加速度幅值约束：|acceleration| <= maxAcc
    // 使用二阶有限差分：acceleration = n*(n-1)/T * (P_{i-1} - 2*P_i + P_{i+1})
    if(ENFORCE_ACC)
    {
        for(int k = 0; k < segment_num ; k ++ )
        {
            for(int i = 0; i < 3; i++)                      // 对于每个维度 (x, y, z)
            { 
                for(int p = 0; p < traj_order - 1; p++)     // 对于每个加速度约束点
                {    
                    int nzi = 3;                            // 每个约束3个非零元素
                    MSKint32t asub[nzi];
                    double aval[nzi];

                    // 加速度有限差分的系数：n*(n-1)/T * (P_{i-1} - 2*P_i + P_{i+1})
                    double acc_coeff = traj_order * (traj_order - 1) / corridor[k].t;
                    aval[0] =  1.0 * acc_coeff;             // +系数 * P_{i-1}
                    aval[1] = -2.0 * acc_coeff;             // -2*系数 * P_i
                    aval[2] =  1.0 * acc_coeff;             // +系数 * P_{i+1}
                    
                    // 控制点P_{i-1}, P_i, P_{i+1}的变量索引
                    asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    
                    asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;    
                    asub[2] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 2;    
                    
                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                }
            }
        }
    }
    // 3. 起点边界条件
    {
        // 起点位置约束：P_0 = start_position (由initScale缩放)
        for(int i = 0; i < 3; i++)          // 对于每个维度 (x, y, z)
        {        
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = 1.0 * initScale;      // 该段的缩放因子
            asub[0] = i * s1d1CtrlP_num;    // 维度i的第一个控制点
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        
        // 起点速度约束：velocity = n * (P_1 - P_0) = start_velocity
        for(int i = 0; i < 3; i++)          // 对于每个维度 (x, y, z)
        {       
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = - 1.0 * traj_order;   // -n * P_0
            aval[1] =   1.0 * traj_order;   // +n * P_1
            asub[0] = i * s1d1CtrlP_num;
            asub[1] = i * s1d1CtrlP_num + 1;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);   
            row_idx ++;
        }
        
        // 起点加速度约束：acceleration = n*(n-1)/T * (P_0 - 2*P_1 + P_2) = start_acceleration
        for(int i = 0; i < 3; i++)          // 对于每个维度 (x, y, z)
        {       
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            double acc_coeff = traj_order * (traj_order - 1) / initScale;
            aval[0] =   1.0 * acc_coeff;    // +系数 * P_0
            aval[1] = - 2.0 * acc_coeff;    // -2*系数 * P_1
            aval[2] =   1.0 * acc_coeff;    // +系数 * P_2
            asub[0] = i * s1d1CtrlP_num;
            asub[1] = i * s1d1CtrlP_num + 1;
            asub[2] = i * s1d1CtrlP_num + 2;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }      

    // 4. End point boundary conditions
    {   
        // End position constraints: P_final = end_position (scaled by lstScale)
        for(int i = 0; i < 3; i++)          // for each dimension (x, y, z)
        {       
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            // Index of the last control point for dimension i (control points are stored as [x0,x1,...,xn, y0,y1,...,yn, z0,z1,...,zn] for last segment)
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] = 1.0 * lstScale;       // scale factor for last segment
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        
        // End velocity constraints: velocity = n * (P_final - P_{final-1}) = end_velocity
        for(int i = 0; i < 3; i++)          // for each dimension (x, y, z)
        { 
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;  // P_{final-1}
            asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;      // P_final
            aval[0] = - 1.0;                // -P_{final-1}
            aval[1] =   1.0;                // +P_final
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        
        // End acceleration constraints: acceleration = n*(n-1)/T * (P_{final-2} - 2*P_{final-1} + P_final) = end_acceleration
        for(int i = 0; i < 3; i++)          // for each dimension (x, y, z)
        { 
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 2;  // P_{final-2}
            asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;  // P_{final-1}
            asub[2] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;      // P_final
            double acc_coeff = 1.0 / lstScale;
            aval[0] =   1.0 * acc_coeff;    // +coefficient * P_{final-2}
            aval[1] = - 2.0 * acc_coeff;    // -2*coefficient * P_{final-1}
            aval[2] =   1.0 * acc_coeff;    // +coefficient * P_final
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    // 5. Continuity constraints at segment joints 
    // Ensures position, velocity, and acceleration are continuous between adjacent segments
    {
        int sub_shift = 0;                      // offset for indexing control points of different segments
        double val0, val1;
        for(int k = 0; k < (segment_num - 1); k ++ )    // for each joint between segments
        {   
            double scale_k = corridor[k].t;     // time scale of current segment k
            double scale_n = corridor[k+1].t;   // time scale of next segment k+1
            
            // Position continuity: T_k * P_last_k = T_{k+1} * P_first_{k+1}
            val0 = scale_k;
            val1 = scale_n;
            for(int i = 0; i < 3; i++)          // for each dimension (x, y, z)
            {  
                int nzi = 2;
                MSKint32t asub[nzi];
                double aval[nzi];

                // This segment's last control point (scaled)
                aval[0] = 1.0 * val0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 1;

                // Next segment's first control point (scaled) 
                aval[1] = -1.0 * val1;
                asub[1] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
            
            // Velocity continuity: n*(P_last_k - P_{last-1}_k) = n*(P_{first+1}_{k+1} - P_first_{k+1})
            for(int i = 0; i < 3; i++)          // for each dimension (x, y, z)
            {  
                int nzi = 4;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                // This segment's velocity at end: n * (P_last - P_{last-1})
                aval[0] = -1.0;                 // -P_{last-1}
                aval[1] =  1.0;                 // +P_last
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 2;    
                asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
                
                // Next segment's velocity at start: n * (P_{first+1} - P_first)
                aval[2] =  1.0;                 // +P_first
                aval[3] = -1.0;                 // -P_{first+1}
                asub[2] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
                asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
            
            // Acceleration continuity: 
            // n*(n-1)/T_k * (P_{last-2} - 2*P_{last-1} + P_last) = n*(n-1)/T_{k+1} * (P_first - 2*P_{first+1} + P_{first+2})
            val0 = 1.0 / scale_k;               // inverse time scale for current segment
            val1 = 1.0 / scale_n;               // inverse time scale for next segment
            for(int i = 0; i < 3; i++)          // for each dimension (x, y, z)
            {  
                int nzi = 6;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                // This segment's acceleration at end: n*(n-1)/T_k * (P_{last-2} - 2*P_{last-1} + P_last)
                aval[0] =  1.0  * val0;         // +coefficient/T_k * P_{last-2}
                aval[1] = -2.0  * val0;         // -2*coefficient/T_k * P_{last-1}
                aval[2] =  1.0  * val0;         // +coefficient/T_k * P_last
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 3;    
                asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 2;   
                asub[2] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
                
                // Next segment's acceleration at start: n*(n-1)/T_{k+1} * (P_first - 2*P_{first+1} + P_{first+2})
                aval[3] =  -1.0  * val1;        // -coefficient/T_{k+1} * P_first
                aval[4] =   2.0  * val1;        // +2*coefficient/T_{k+1} * P_{first+1}
                aval[5] =  -1.0  * val1;        // -coefficient/T_{k+1} * P_{first+2}
                asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
                asub[4] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;
                asub[5] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 2;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }

            sub_shift += s1CtrlP_num;           // move to next segment's control points
        }
    }

/*
 * ===================================================================
 * 二次目标函数设置
 * ===================================================================
 * 
 * 设置二次目标函数：最小化 (1/2) * x^T * Q * x
 * 其中Q由预计算的最小化矩阵MQM构建。
 * 目标函数最小化指定的导数阶数（速度、加速度、急动度或snap）
 * 在所有轨迹段上，并按时间适当缩放。
 */
    
    // 通过整数阶之间插值处理分数最小化阶
    int min_order_l = floor(minimize_order);    // 下整数界
    int min_order_u = ceil (minimize_order);    // 上整数界

    // 计算二次目标矩阵Q中非零元素的数量
    // Q是对称的，所以我们只存储上三角部分
    int NUMQNZ = 0;
    for(int i = 0; i < segment_num; i ++)
    {   
        int NUMQ_blk = (traj_order + 1);        // 每段每维度的控制点数量
        NUMQNZ      += 3 * NUMQ_blk * (NUMQ_blk + 1) / 2;   // 3个维度 × 上三角元素
    }
    MSKint32t  qsubi[NUMQNZ], qsubj[NUMQNZ];   // Q矩阵的行和列索引
    double     qval[NUMQNZ];                   // Q矩阵的值
    
    // Build the quadratic objective matrix Q
    {    
        int sub_shift = 0;                      // offset for current segment's control points
        int idx = 0;                           // index into Q matrix arrays
        for(int k = 0; k < segment_num; k ++)  // for each trajectory segment
        {
            double scale_k = corridor[k].t;     // time scale for this segment
            for(int p = 0; p < 3; p ++ )        // for each dimension (x, y, z)
                for( int i = 0; i < s1d1CtrlP_num; i ++ )      // for each control point i
                    for( int j = 0; j < s1d1CtrlP_num; j ++ )  // for each control point j  
                        if( i >= j )                           // only upper triangular (Q is symmetric)
                        {
                            qsubi[idx] = sub_shift + p * s1d1CtrlP_num + i;   // row index in Q
                            qsubj[idx] = sub_shift + p * s1d1CtrlP_num + j;   // column index in Q
                            
                            // Scale the pre-computed MQM matrix by appropriate time factors
                            // Time scaling: T^(2*order-3) where order is the derivative being minimized
                            if(min_order_l == min_order_u)
                                // Integer minimize order case
                                qval[idx]  = MQM(i, j) /(double)pow(scale_k, 2 * min_order_u - 3);
                            else
                                // Fractional minimize order: interpolate between adjacent integer orders
                                qval[idx] = ( (minimize_order - min_order_l) / (double)pow(scale_k, 2 * min_order_u - 3)
                                            + (min_order_u - minimize_order) / (double)pow(scale_k, 2 * min_order_l - 3) ) * MQM(i, j);
                            idx ++ ;
                        }

            sub_shift += s1CtrlP_num;          // move to next segment
        }
    }
         
    ros::Time time_end1 = ros::Time::now();

    // Set the quadratic objective in MOSEK
    if ( r== MSK_RES_OK )
         r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval); 
    
    // Set objective sense to minimization
    if ( r==MSK_RES_OK ) 
         r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
    
/*
 * ===================================================================
 * 优化求解和结果处理
 * ===================================================================
 */
    
    bool solve_ok = false;
    if ( r==MSK_RES_OK ) 
      { 
        // Solve the optimization problem using MOSEK's interior-point method
        MSKrescodee trmcode; 
        r = MSK_optimizetrm(task,&trmcode); 
        MSK_solutionsummary (task,MSK_STREAM_LOG);      // Print solution summary for debugging
          
        if ( r==MSK_RES_OK ) 
        { 
          // Check the solution status to determine if optimization was successful
          MSKsolstae solsta; 
          MSK_getsolsta (task,MSK_SOL_ITR,&solsta); 
           
          switch(solsta) 
          { 
            case MSK_SOL_STA_OPTIMAL:           // Optimal solution found
            case MSK_SOL_STA_NEAR_OPTIMAL:      // Near-optimal solution found (acceptable)
              
            // Extract the optimal control point values
            r = MSK_getxx(task, 
                          MSK_SOL_ITR,          // Request the interior-point solution
                          x_var); 
            
            // Extract the optimal objective value
            r = MSK_getprimalobj(
                task,
                MSK_SOL_ITR,
                &primalobj);

            obj = primalobj;                    // Return objective value to caller
            solve_ok = true;                    // Mark successful solution
            
            break; 
            
            // Handle infeasible problem cases
            case MSK_SOL_STA_DUAL_INFEAS_CER: 
            case MSK_SOL_STA_PRIM_INFEAS_CER: 
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER: 
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:   
              printf("Primal or dual infeasibility certificate found.\n"); 
              break; 
               
            case MSK_SOL_STA_UNKNOWN: 
              printf("The status of the solution could not be determined.\n"); 
              //solve_ok = true; // debug
              break; 
            default: 
              printf("Other solution status."); 
              break; 
          } 
        } 
        else 
        { 
          printf("Error while optimizing.\n"); 
        } 
      }
     
      // Handle any MOSEK errors
      if (r != MSK_RES_OK) 
      { 
        char symname[MSK_MAX_STR_LEN]; 
        char desc[MSK_MAX_STR_LEN]; 
         
        printf("An error occurred while optimizing.\n");      
        MSK_getcodedesc (r, 
                         symname, 
                         desc); 
        printf("Error %s - '%s'\n",symname,desc); 
      } 
    
    // Clean up MOSEK resources
    MSK_deletetask(&task); 
    MSK_deleteenv(&env); 

    // Log optimization timing
    ros::Time time_end2 = ros::Time::now();
    ROS_WARN("time consume in optimize is :");
    cout<<time_end2 - time_end1<<endl;

    // Check if optimization was successful
    if(!solve_ok){
      ROS_WARN("In solver, falied ");
      return -1;
    }

/*
 * ===================================================================
 * 结果转换和输出
 * ===================================================================
 */
    
    // 将解从MOSEK格式转换为Eigen格式
    VectorXd d_var(ctrlP_num);
    for(int i = 0; i < ctrlP_num; i++)
        d_var(i) = x_var[i];
    
    // 将解重塑为系数矩阵：
    // 每行代表一段，列为 [x0,x1,...,xn, y0,y1,...,yn, z0,z1,...,zn]
    PolyCoeff = MatrixXd::Zero(segment_num, 3 *(traj_order + 1) );

    int var_shift = 0;
    for(int i = 0; i < segment_num; i++ )       // 对于每个轨迹段
    {
        for(int j = 0; j < 3 * n_poly; j++)     // 对于每个系数 (x,y,z维度)
            PolyCoeff(i , j) = d_var(j + var_shift);

        var_shift += 3 * n_poly;               // 移动到下一段的系数
    }   

    return 1;                                  // 成功
}