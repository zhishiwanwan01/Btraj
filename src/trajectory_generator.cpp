/**
 * @file trajectory_generator.cpp
 * @brief 基于MOSEK二次规划的贝塞尔轨迹生成器
 *
 * 本文件实现了Btraj规划器的核心轨迹优化算法。该算法将轨迹规划问题
 * 转化为约束二次规划(Constrained Quadratic Programming, CQP)问题，
 * 通过MOSEK商业优化求解器求解，生成平滑且满足动力学约束的贝塞尔轨迹。
 *
 * ===== 二次规划问题数学形式 =====
 * minimize:   (1/2) * x^T * Q * x + c^T * x     (目标函数: 最小化snap/jerk)
 * subject to: A_eq * x = b_eq                   (等式约束: 边界条件+连续性)
 *             A_ineq * x ≤ b_ineq               (不等式约束: 速度+加速度限制)
 *             lb ≤ x ≤ ub                       (盒约束: 空间走廊限制)
 *
 * 其中 x 为贝塞尔控制点向量: [seg0_x, seg0_y, seg0_z, seg1_x, seg1_y, ...]
 *
 * ===== 主要功能模块 =====
 * 1. 问题规模计算: 确定变量数量和约束数量
 * 2. 约束边界设置: 定义MOSEK约束边界类型和数值
 * 3. 约束矩阵构建: 构建线性约束矩阵A和右端向量b
 * 4. 目标函数设置: 构建二次目标矩阵Q
 * 5. MOSEK求解: 调用优化求解器获得最优解
 * 6. 结果转换: 将MOSEK解转换为贝塞尔系数矩阵
 *
 * ===== MOSEK优化器优势 =====
 * - 内点法求解器: 数值稳定，收敛快速
 * - 商业级精度: 支持高精度求解
 * - 大规模优化: 处理数千变量的问题
 * - 凸优化保证: QP问题具有唯一全局最优解
 *
 * @author UAV轨迹规划团队
 * @date 2023
 */

#include "trajectory_generator.h"
using namespace std;
using namespace Eigen;

/**
 * @brief MOSEK求解器日志回调函数
 *
 * 该回调函数由MOSEK优化器在求解过程中调用，用于输出：
 * - 迭代进度信息 (迭代次数、目标函数值、约束违反度等)
 * - 算法收敛状态 (原始/对偶间隙、可行性等)
 * - 调试和诊断信息
 *
 * MOSEK内点法求解器输出示例:
 * "Iter     Primal obj        Dual obj       ..."
 * "   0   1.000000e+00    0.000000e+00    ..."
 * "   1   5.234567e-01    4.987654e-01    ..."
 *
 * @param handle MOSEK环境句柄(未使用)
 * @param str 需要输出的字符串消息
 */
static void MSKAPI printstr(void* handle, MSKCONST char str[]) {
  printf("%s", str);  // 直接输出MOSEK日志到控制台
}

/**
 * @brief 主要贝塞尔轨迹生成函数 - 二次规划求解核心
 *
 * 该函数是整个轨迹优化的核心，将轨迹规划问题转化为标准的
 * 约束二次规划(CQP)形式，并调用MOSEK商业求解器获得最优解。
 *
 * ===== 算法流程概述 =====
 * 1. 问题建模: 将贝塞尔轨迹优化转化为QP标准形式
 * 2. 约束构建: 设置边界条件、连续性、速度/加速度限制
 * 3. 目标函数: 构建最小化snap/jerk的二次目标矩阵
 * 4. MOSEK求解: 调用内点法求解器获得全局最优解
 * 5. 结果提取: 将优化解转换为贝塞尔系数矩阵格式
 *
 * ===== 数学问题描述 =====
 * 决策变量 x: 所有段所有维度的贝塞尔控制点
 *   x = [P0_x0, P0_x1, ..., P0_xn, P0_y0, ..., P0_zn, P1_x0, P1_x1, ...]
 *   维度: (segment_num × 3 × (order+1)) × 1
 *
 * 目标函数: minimize (1/2) * x^T * Q * x
 *   Q矩阵: 由MQM预计算矩阵和时间缩放构成
 *   目的: 最小化轨迹的高阶导数(通常是snap或jerk)
 *
 * @param corridor 飞行走廊立方体序列，定义空间约束和时间分配
 * @param MQM 预计算的最小化矩阵(来自bezier_base)，用于构建Q矩阵
 * @param pos 边界位置条件 [起点位置; 终点位置] (2×3矩阵)
 * @param vel 边界速度条件 [起点速度; 终点速度] (2×3矩阵)
 * @param acc 边界加速度条件 [起点加速度; 终点加速度] (2×3矩阵)
 * @param maxVel 速度幅值上限 (m/s)，用于动力学约束
 * @param maxAcc 加速度幅值上限 (m/s²)，用于动力学约束
 * @param traj_order 贝塞尔曲线阶数 (通常5-12)，影响轨迹平滑度
 * @param minimize_order 最小化的导数阶数 (1=速度, 2=加速度, 3=jerk, 4=snap)
 * @param margin 走廊内安全裕度 (m)，避免贴边飞行
 * @param isLimitVel 是否启用速度约束 (影响约束矩阵规模)
 * @param isLimitAcc 是否启用加速度约束 (影响约束矩阵规模)
 * @param obj 输出参数: 优化目标函数的最优值
 * @param PolyCoeff 输出参数: 贝塞尔系数矩阵 (segment_num × 3×(order+1))
 *
 * @return int 1=成功, -1=失败 (MOSEK求解失败或参数错误)
 */
int TrajectoryGenerator::BezierPloyCoeffGeneration(
    const vector<Cube>& corridor,
    const MatrixXd& MQM,
    const MatrixXd& pos,
    const MatrixXd& vel,
    const MatrixXd& acc,
    const double maxVel,
    const double maxAcc,
    const int traj_order,
    const double minimize_order,
    const double margin,
    const bool& isLimitVel,
    const bool& isLimitAcc,
    double& obj,
    MatrixXd& PolyCoeff)  // define the order to which we minimize.   1 --
                          // velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap
{
  /**
   * ===================================================================
   * 第1步: 问题规模计算和MOSEK变量初始化
   * ===================================================================
   *
   * 本节计算二次规划问题的规模，包括：
   * - 决策变量数量 (贝塞尔控制点总数)
   * - 约束条件数量 (等式约束 + 不等式约束)
   * - MOSEK API所需的数据结构初始化
   */
#define ENFORCE_VEL isLimitVel  // 宏定义：是否启用速度约束(影响问题规模)
#define ENFORCE_ACC isLimitAcc  // 宏定义：是否启用加速度约束(影响问题规模)

  // ===== 时间和几何参数提取 =====
  double initScale = corridor.front().t;  // 第一段轨迹的时间缩放因子
  double lstScale = corridor.back().t;    // 最后一段轨迹的时间缩放因子
  int segment_num = corridor.size();      // 轨迹段数量

  // ===== 贝塞尔控制点维度计算 =====
  int n_poly = traj_order + 1;          // n_poly = number of polynomial coefficients: 每段每维度控制点数 (n阶=n+1个点)
  int s1d1CtrlP_num = n_poly;           // s1d1CtrlP_num = single segment 1-dimensional control point number: 一段一维控制点数量
  int s1CtrlP_num = 3 * s1d1CtrlP_num;  // s1CtrlP_num = single segment control point number: 一段三维控制点数量 (x,y,z)

  // ===== 等式约束数量统计 =====
  // 边界条件约束: 起点和终点的位置、速度、加速度
  // ? 为什么约束数量等于维度和导数的乘积
  int equ_con_s_num = 3 * 3;  // equ_con_s_num = equality constraint start number: 起点: 3维度 × 3阶导数(pos,vel,acc) = 9个约束
  int equ_con_e_num = 3 * 3;  // equ_con_e_num = equality constraint end number: 终点: 3维度 × 3阶导数(pos,vel,acc) = 9个约束

  // 连续性约束: 相邻段在连接点处的C²连续性
  int equ_con_continuity_num =
      3 * 3 * (segment_num - 1);  // equ_con_continuity_num = equality constraint continuity number: 3维×3阶导数×(段数-1)

  // 总等式约束数 = 边界条件 + 连续性
  int equ_con_num = equ_con_s_num + equ_con_e_num + equ_con_continuity_num;  // equ_con_num = equality constraint number

  // ===== 不等式约束数量统计 =====
  // 速度幅值约束: |vel| ≤ maxVel (每个轨迹点每个维度)
  int vel_con_num = 3 * traj_order * segment_num;  // vel_con_num = velocity constraint number: 3维×n阶×段数

  // 加速度幅值约束: |acc| ≤ maxAcc (每个轨迹点每个维度)
  int acc_con_num = 3 * (traj_order - 1) * segment_num;  // acc_con_num = acceleration constraint number: 3维×(n-1)阶×段数

  // 根据用户设置决定是否包含动力学约束
  if (!ENFORCE_VEL)
    vel_con_num = 0;  // 禁用速度约束

  if (!ENFORCE_ACC)
    acc_con_num = 0;  // 禁用加速度约束

  int high_order_con_num = vel_con_num + acc_con_num;  // high_order_con_num = high order constraint number: 总不等式约束数

  // ===== 问题总规模 =====
  // * Ac = b 中 A 的行数 = 等式约束 + 不等式约束
  int con_num = equ_con_num + high_order_con_num;  // con_num = constraint number: 总约束数量
  // * Ac = b 中 x 的维度，也就是A的列数 = 控制点的数量
  int ctrlP_num = segment_num * s1CtrlP_num;  // ctrlP_num = control point number: 总决策变量数量

  // ===== MOSEK数据结构初始化 =====
  double x_var[ctrlP_num];  // x_var = x variable: 存储最优解的数组
  double primalobj;         // primalobj = primal objective: 存储目标函数最优值
  MSKrescodee r;            // r = response code: MOSEK返回状态码

  // 约束边界数据结构: 存储每个约束的类型和边界值
  // 格式: <边界类型, <下界值, 上界值>>
  vector<pair<MSKboundkeye, pair<double, double> > > con_bdk;  // con_bdk = constraint boundary key: "k" comes from MSKboundkeye (MOSEK bound key enum)

  /**
   * ===================================================================
   * 第2步: MOSEK约束边界设置
   * ===================================================================
   *
   * 本节设置所有约束的边界类型和数值，MOSEK支持以下边界类型：
   * - MSK_BK_FX: 固定约束 (constraint = value)
   * - MSK_BK_RA: 范围约束 (lower ≤ constraint ≤ upper)
   * - MSK_BK_UP: 上界约束 (constraint ≤ upper)
   * - MSK_BK_LO: 下界约束 (constraint ≥ lower)
   *
   * 约束按以下顺序添加:
   * 1. 速度幅值约束 (不等式, MSK_BK_RA)
   * 2. 加速度幅值约束 (不等式, MSK_BK_RA)
   * 3. 边界条件约束 (等式, MSK_BK_FX)
   * 4. 连续性约束 (等式, MSK_BK_FX)
   */

  // ===== 速度幅值约束设置: -maxVel ≤ vel ≤ +maxVel =====
  if (ENFORCE_VEL) {
    for (int i = 0; i < vel_con_num; i++) {
      // MSK_BK_RA: 双边界约束，限制速度在[-maxVel, +maxVel]范围内
      // 这确保无人机不会超过物理速度限制，保证轨迹的动力学可行性
      pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair(  // cb_ie = constraint bound inequality
          MSK_BK_RA,                     // 约束类型: 范围约束
          make_pair(-maxVel, +maxVel));  // 约束范围: [-maxVel, +maxVel]
      con_bdk.push_back(cb_ie);
    }
  }

  // ===== 加速度幅值约束设置: -maxAcc ≤ acc ≤ +maxAcc =====
  if (ENFORCE_ACC) {
    for (int i = 0; i < acc_con_num; i++) {
      // MSK_BK_RA: 双边界约束，限制加速度在[-maxAcc, +maxAcc]范围内
      // 这确保无人机不会超过最大推力限制，避免电机饱和
      pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair(  // cb_ie = constraint bound inequality
          MSK_BK_RA,                    // 约束类型: 范围约束
          make_pair(-maxAcc, maxAcc));  // 约束范围: [-maxAcc, +maxAcc]
      con_bdk.push_back(cb_ie);
    }
  }

  // ===== 等式约束设置: 边界条件 + 连续性约束 =====
  // 约束顺序严格按照: 起点位置(3) → 起点速度(3) → 起点加速度(3) →
  //                 终点位置(3) → 终点速度(3) → 终点加速度(3) → 连续性约束(...)
  for (int i = 0; i < equ_con_num; i++) {
    double beq_i;  // beq_i = b equality i: 等式约束的右端值

    // 根据约束索引确定约束类型和数值
    if (i < 3)
      beq_i = pos(0, i);  // 起点位置 (x,y,z)
    else if (i >= 3 && i < 6)
      beq_i = vel(0, i - 3);  // 起点速度 (x,y,z)
    else if (i >= 6 && i < 9)
      beq_i = acc(0, i - 6);  // 起点加速度 (x,y,z)
    else if (i >= 9 && i < 12)
      beq_i = pos(1, i - 9);  // 终点位置 (x,y,z)
    else if (i >= 12 && i < 15)
      beq_i = vel(1, i - 12);  // 终点速度 (x,y,z)
    else if (i >= 15 && i < 18)
      beq_i = acc(1, i - 15);  // 终点加速度 (x,y,z)
    else
      beq_i = 0.0;  // 连续性约束: 段间差值=0

    // MSK_BK_FX: 固定值约束，要求 constraint = beq_i (上下界相等)
    pair<MSKboundkeye, pair<double, double> > cb_eq =  // cb_eq = constraint bound equality
        make_pair(MSK_BK_FX,                 // 约束类型: 固定值约束
                  make_pair(beq_i, beq_i));  // 约束值: 上下界均为beq_i
    con_bdk.push_back(cb_eq);
  }

  /**
   * ===================================================================
   * 第3步: 变量边界设置 - 空间走廊约束
   * ===================================================================
   *
   * 本节设置决策变量(贝塞尔控制点)的边界约束，确保轨迹始终在
   * 安全飞行走廊内。这是"盒约束"的实现: lb ≤ x ≤ ub
   *
   * 空间约束的数学意义:
   * - 每个控制点必须在对应走廊立方体内
   * - 控制点坐标经过时间缩放: actual_pos = scale_k * control_point
   * - 因此约束形式为: (box_min/scale_k) ≤ control_point ≤ (box_max/scale_k)
   */
  vector<pair<MSKboundkeye, pair<double, double> > > var_bdk;  // var_bdk = variable boundary key: "k" comes from MSKboundkeye (MOSEK bound key enum)

  // 遍历每个轨迹段，为其控制点设置空间边界
  for (int k = 0; k < segment_num; k++) {
    Cube cube_ = corridor[k];  // 当前段的走廊立方体
    double scale_k = cube_.t;  // 当前段的时间缩放因子

    // 对于每个空间维度 (x, y, z)
    for (int i = 0; i < 3; i++) {
      // 对于该段该维度的每个贝塞尔控制点
      for (int j = 0; j < n_poly; j++) {
        pair<MSKboundkeye, pair<double, double> > vb_x;  // vb_x = variable bound x

        double lo_bound, up_bound;  // lo_bound/up_bound = lower/upper bound: 控制点的上下界

        // 安全裕度处理: 第一段保持精确边界，其他段应用安全裕度
        if (k > 0) {
          // 其他段: 应用安全裕度，避免贴边飞行
          // actual_constraint: (box_min + margin) ≤ scale_k * x ≤ (box_max -
          // margin) 转化为控制点约束: lo_bound ≤ x ≤ up_bound
          lo_bound = (cube_.box[i].first + margin) / scale_k;
          up_bound = (cube_.box[i].second - margin) / scale_k;
        } else {
          // 第一段: 不应用安全裕度，保证精确到达起点
          lo_bound = (cube_.box[i].first) / scale_k;
          up_bound = (cube_.box[i].second) / scale_k;
        }

        // 设置范围约束: lo_bound ≤ control_point ≤ up_bound
        vb_x =
            make_pair(MSK_BK_RA,                       // 边界类型: 范围约束
                      make_pair(lo_bound, up_bound));  // 边界值: [下界, 上界]

        var_bdk.push_back(vb_x);  // 添加到变量边界向量
      }
    }
  }

  /**
   * ===================================================================
   * 第4步: MOSEK优化环境初始化和参数配置
   * ===================================================================
   *
   * 本节创建MOSEK优化任务并配置求解器参数，包括：
   * 1. 创建MOSEK环境和任务实例
   * 2. 设置内点法求解器的数值精度参数
   * 3. 分配问题的约束和变量空间
   * 4. 将之前设置的边界条件传递给MOSEK
   *
   * MOSEK内点法求解器特点:
   * - 数值稳定：采用预条件共轭梯度法
   * - 收敛快速：通常10-30次迭代即可收敛
   * - 高精度：支持1e-8级别的最优性条件检查
   */
  MSKint32t j, i;  // j, i = indices: MOSEK循环索引变量
  MSKenv_t env;    // env = environment: MOSEK环境句柄
  MSKtask_t task;  // task = optimization task: MOSEK优化任务句柄

  // ===== 创建MOSEK优化环境和任务 =====
  r = MSK_makeenv(&env, NULL);                       // 创建MOSEK环境
  r = MSK_maketask(env, con_num, ctrlP_num, &task);  // 创建优化任务

  // ===== 配置MOSEK内点法求解器参数 =====
  // 这些参数控制求解器的数值精度和性能平衡
  MSK_putintparam(task, MSK_IPAR_NUM_THREADS, 1);  // 线程数：单线程避免并发开销
  MSK_putdouparam(task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL,
                  1e-2);                                    // 凸性检查相对容差
  MSK_putdouparam(task, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-4);   // 对偶可行性容差
  MSK_putdouparam(task, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-4);   // 原始可行性容差
  MSK_putdouparam(task, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-4);  // 不可行性检测容差
  // MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 5e-2 );      //
  // 相对间隙容差（已注释，使用默认值）

  // r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);  //
  // 启用求解器日志输出（已注释）

  // ===== 分配优化问题的空间 =====
  if (r == MSK_RES_OK)
    r = MSK_appendcons(task, con_num);  // 分配约束空间 (con_num个约束)
  if (r == MSK_RES_OK)
    r = MSK_appendvars(task,
                       ctrlP_num);  // 分配变量空间 (ctrlP_num个变量，初始值为0)

  // ===== 设置变量边界约束 (盒约束: lb ≤ x ≤ ub) =====
  for (j = 0; j < ctrlP_num && r == MSK_RES_OK; ++j) {
    if (r == MSK_RES_OK)
      r = MSK_putvarbound(task,
                          j,                          // 变量索引 (第j个控制点)
                          var_bdk[j].first,           // 边界类型 (MSK_BK_RA)
                          var_bdk[j].second.first,    // 下界值 (空间走廊下界)
                          var_bdk[j].second.second);  // 上界值 (空间走廊上界)
  }

  // ===== 设置约束边界 (线性约束: Ax ≤ b 或 Ax = b) =====
  for (i = 0; i < con_num && r == MSK_RES_OK; i++) {
    r = MSK_putconbound(task,
                        i,                 // 约束索引
                        con_bdk[i].first,  // 边界类型 (MSK_BK_FX或MSK_BK_RA)
                        con_bdk[i].second.first,    // 下界值
                        con_bdk[i].second.second);  // 上界值
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
  int row_idx = 0;  // row_idx = row index: 当前处理的约束矩阵行索引

  // 1. 速度幅值约束：|velocity| <= maxVel
  // 使用有限差分：velocity = n * (P_{i+1} - P_i)，其中n为轨迹阶数
  if (ENFORCE_VEL) {
    for (int k = 0; k < segment_num; k++) {
      for (int i = 0; i < 3; i++)  // 对于每个维度 (x, y, z)
      {
        for (int p = 0; p < traj_order; p++)  // 对于每个速度约束点
        {
          int nzi = 2;  // nzi = number of zeros i: 每个约束2个非零元素
          MSKint32t asub[nzi];  // asub = A matrix subscript indices
          double aval[nzi];     // aval = A matrix values

          // 速度有限差分的系数：n * (P_{i+1} - P_i)
          aval[0] = -1.0 * traj_order;  // -n * P_i
          aval[1] = 1.0 * traj_order;   // +n * P_{i+1}

          // 控制点P_i和P_{i+1}的变量索引
          asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;
          asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;

          r = MSK_putarow(task, row_idx, nzi, asub, aval);
          row_idx++;
        }
      }
    }
  }

  // 2. 加速度幅值约束：|acceleration| <= maxAcc
  // 使用二阶有限差分：acceleration = n*(n-1)/T * (P_{i-1} - 2*P_i + P_{i+1})
  if (ENFORCE_ACC) {
    for (int k = 0; k < segment_num; k++) {
      for (int i = 0; i < 3; i++)  // 对于每个维度 (x, y, z)
      {
        for (int p = 0; p < traj_order - 1; p++)  // 对于每个加速度约束点
        {
          int nzi = 3;  // nzi = number of zeros i: 每个约束3个非零元素
          MSKint32t asub[nzi];  // asub = A matrix subscript indices
          double aval[nzi];     // aval = A matrix values

          // 加速度有限差分的系数：n*(n-1)/T * (P_{i-1} - 2*P_i + P_{i+1})
          double acc_coeff = traj_order * (traj_order - 1) / corridor[k].t;  // acc_coeff = acceleration coefficient
          aval[0] = 1.0 * acc_coeff;   // +系数 * P_{i-1}
          aval[1] = -2.0 * acc_coeff;  // -2*系数 * P_i
          aval[2] = 1.0 * acc_coeff;   // +系数 * P_{i+1}

          // 控制点P_{i-1}, P_i, P_{i+1}的变量索引
          asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;
          asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;
          asub[2] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 2;

          r = MSK_putarow(task, row_idx, nzi, asub, aval);
          row_idx++;
        }
      }
    }
  }
  // 3. 起点边界条件
  {
    // 起点位置约束：P_0 = start_position (由initScale缩放)
    for (int i = 0; i < 3; i++)  // 对于每个维度 (x, y, z)
    {
      int nzi = 1;  // nzi = number of zeros i
      MSKint32t asub[nzi];  // asub = A matrix subscript indices
      double aval[nzi];     // aval = A matrix values
      aval[0] = 1.0 * initScale;    // 该段的缩放因子
      asub[0] = i * s1d1CtrlP_num;  // 维度i的第一个控制点
      r = MSK_putarow(task, row_idx, nzi, asub, aval);
      row_idx++;
    }

    // 起点速度约束：velocity = n * (P_1 - P_0) = start_velocity
    for (int i = 0; i < 3; i++)  // 对于每个维度 (x, y, z)
    {
      int nzi = 2;
      MSKint32t asub[nzi];
      double aval[nzi];
      aval[0] = -1.0 * traj_order;  // -n * P_0
      aval[1] = 1.0 * traj_order;   // +n * P_1
      asub[0] = i * s1d1CtrlP_num;
      asub[1] = i * s1d1CtrlP_num + 1;
      r = MSK_putarow(task, row_idx, nzi, asub, aval);
      row_idx++;
    }

    // 起点加速度约束：acceleration = n*(n-1)/T * (P_0 - 2*P_1 + P_2) =
    // start_acceleration
    for (int i = 0; i < 3; i++)  // 对于每个维度 (x, y, z)
    {
      int nzi = 3;
      MSKint32t asub[nzi];
      double aval[nzi];
      double acc_coeff = traj_order * (traj_order - 1) / initScale;  // acc_coeff = acceleration coefficient
      aval[0] = 1.0 * acc_coeff;   // +系数 * P_0
      aval[1] = -2.0 * acc_coeff;  // -2*系数 * P_1
      aval[2] = 1.0 * acc_coeff;   // +系数 * P_2
      asub[0] = i * s1d1CtrlP_num;
      asub[1] = i * s1d1CtrlP_num + 1;
      asub[2] = i * s1d1CtrlP_num + 2;
      r = MSK_putarow(task, row_idx, nzi, asub, aval);
      row_idx++;
    }
  }

  // 4. End point boundary conditions
  {
    // End position constraints: P_final = end_position (scaled by lstScale)
    for (int i = 0; i < 3; i++)  // for each dimension (x, y, z)
    {
      int nzi = 1;
      MSKint32t asub[nzi];
      double aval[nzi];
      // Index of the last control point for dimension i (control points are
      // stored as [x0,x1,...,xn, y0,y1,...,yn, z0,z1,...,zn] for last segment)
      asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
      aval[0] = 1.0 * lstScale;  // scale factor for last segment
      r = MSK_putarow(task, row_idx, nzi, asub, aval);
      row_idx++;
    }

    // End velocity constraints: velocity = n * (P_final - P_{final-1}) =
    // end_velocity
    for (int i = 0; i < 3; i++)  // for each dimension (x, y, z)
    {
      int nzi = 2;
      MSKint32t asub[nzi];
      double aval[nzi];
      asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;  // P_{final-1}
      asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;      // P_final
      aval[0] = -1.0;                                         // -P_{final-1}
      aval[1] = 1.0;                                          // +P_final
      r = MSK_putarow(task, row_idx, nzi, asub, aval);
      row_idx++;
    }

    // End acceleration constraints: acceleration = n*(n-1)/T * (P_{final-2} -
    // 2*P_{final-1} + P_final) = end_acceleration
    for (int i = 0; i < 3; i++)  // for each dimension (x, y, z)
    {
      int nzi = 3;
      MSKint32t asub[nzi];
      double aval[nzi];
      asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 2;  // P_{final-2}
      asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;  // P_{final-1}
      asub[2] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;      // P_final
      double acc_coeff = 1.0 / lstScale;  // acc_coeff = acceleration coefficient
      aval[0] = 1.0 * acc_coeff;   // +coefficient * P_{final-2}
      aval[1] = -2.0 * acc_coeff;  // -2*coefficient * P_{final-1}
      aval[2] = 1.0 * acc_coeff;   // +coefficient * P_final
      r = MSK_putarow(task, row_idx, nzi, asub, aval);
      row_idx++;
    }
  }

  // 5. Continuity constraints at segment joints
  // Ensures position, velocity, and acceleration are continuous between
  // adjacent segments
  {
    int sub_shift =
        0;  // sub_shift = subscript shift: offset for indexing control points of different segments
    double val0, val1;  // val0, val1 = values for scaling coefficients
    for (int k = 0; k < (segment_num - 1);
         k++)  // for each joint between segments
    {
      double scale_k = corridor[k].t;      // time scale of current segment k
      double scale_n = corridor[k + 1].t;  // time scale of next segment k+1

      // Position continuity: T_k * P_last_k = T_{k+1} * P_first_{k+1}
      val0 = scale_k;
      val1 = scale_n;
      for (int i = 0; i < 3; i++)  // for each dimension (x, y, z)
      {
        int nzi = 2;
        MSKint32t asub[nzi];
        double aval[nzi];

        // This segment's last control point (scaled)
        aval[0] = 1.0 * val0;
        asub[0] = sub_shift + (i + 1) * s1d1CtrlP_num - 1;

        // Next segment's first control point (scaled)
        aval[1] = -1.0 * val1;
        asub[1] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;
        r = MSK_putarow(task, row_idx, nzi, asub, aval);
        row_idx++;
      }

      // Velocity continuity: n*(P_last_k - P_{last-1}_k) = n*(P_{first+1}_{k+1}
      // - P_first_{k+1})
      for (int i = 0; i < 3; i++)  // for each dimension (x, y, z)
      {
        int nzi = 4;
        MSKint32t asub[nzi];
        double aval[nzi];

        // This segment's velocity at end: n * (P_last - P_{last-1})
        aval[0] = -1.0;  // -P_{last-1}
        aval[1] = 1.0;   // +P_last
        asub[0] = sub_shift + (i + 1) * s1d1CtrlP_num - 2;
        asub[1] = sub_shift + (i + 1) * s1d1CtrlP_num - 1;

        // Next segment's velocity at start: n * (P_{first+1} - P_first)
        aval[2] = 1.0;   // +P_first
        aval[3] = -1.0;  // -P_{first+1}
        asub[2] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;
        asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;

        r = MSK_putarow(task, row_idx, nzi, asub, aval);
        row_idx++;
      }

      // Acceleration continuity:
      // n*(n-1)/T_k * (P_{last-2} - 2*P_{last-1} + P_last) = n*(n-1)/T_{k+1} *
      // (P_first - 2*P_{first+1} + P_{first+2})
      val0 = 1.0 / scale_k;        // inverse time scale for current segment
      val1 = 1.0 / scale_n;        // inverse time scale for next segment
      for (int i = 0; i < 3; i++)  // for each dimension (x, y, z)
      {
        int nzi = 6;
        MSKint32t asub[nzi];
        double aval[nzi];

        // This segment's acceleration at end: n*(n-1)/T_k * (P_{last-2} -
        // 2*P_{last-1} + P_last)
        aval[0] = 1.0 * val0;   // +coefficient/T_k * P_{last-2}
        aval[1] = -2.0 * val0;  // -2*coefficient/T_k * P_{last-1}
        aval[2] = 1.0 * val0;   // +coefficient/T_k * P_last
        asub[0] = sub_shift + (i + 1) * s1d1CtrlP_num - 3;
        asub[1] = sub_shift + (i + 1) * s1d1CtrlP_num - 2;
        asub[2] = sub_shift + (i + 1) * s1d1CtrlP_num - 1;

        // Next segment's acceleration at start: n*(n-1)/T_{k+1} * (P_first -
        // 2*P_{first+1} + P_{first+2})
        aval[3] = -1.0 * val1;  // -coefficient/T_{k+1} * P_first
        aval[4] = 2.0 * val1;   // +2*coefficient/T_{k+1} * P_{first+1}
        aval[5] = -1.0 * val1;  // -coefficient/T_{k+1} * P_{first+2}
        asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;
        asub[4] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;
        asub[5] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 2;

        r = MSK_putarow(task, row_idx, nzi, asub, aval);
        row_idx++;
      }

      sub_shift += s1CtrlP_num;  // sub_shift = subscript shift: move to next segment's control points
    }
  }

  /**
   * ===================================================================
   * 第5步: 二次目标函数矩阵构建 (核心优化目标)
   * ===================================================================
   *
   * 本节构建二次规划的目标函数: minimize (1/2) * x^T * Q * x
   *
   * ===== 目标函数的物理意义 =====
   * Q矩阵编码了轨迹平滑性的数学定义:
   * - minimize_order = 3: 最小化jerk (急动度) → 平滑轨迹
   * - minimize_order = 4: 最小化snap (snap) → 更平滑轨迹
   * - 更高阶导数 → 更平滑但计算复杂度更高
   *
   * ===== MQM矩阵的来源 =====
   * MQM矩阵来自bezier_base.cpp的预计算，包含:
   * - Bernstein基函数的导数积分
   * - 时间区间[0,1]上的内积计算
   * - 不同阶数贝塞尔曲线的通用形式
   *
   * ===== 时间缩放的必要性 =====
   * 由于轨迹段时间不同，需要按时间缩放:
   * - 物理量级: 更短时间段 → 更大导数值 → 更大惩罚
   * - 缩放公式: Q_scaled = Q / T^(2*order-3)
   * - 确保不同时间段具有一致的优化权重
   */

  // ===== 处理分数阶最小化 (支持如2.5阶导数) =====
  int min_order_l = floor(minimize_order);  // 最小化阶数的下界 (如2.5 → 2)
  int min_order_u = ceil(minimize_order);   // 最小化阶数的上界 (如2.5 → 3)

  // ===== 计算Q矩阵稀疏表示所需空间 =====
  // MOSEK使用三元组格式: (row, col, value) 存储稀疏矩阵
  // 由于Q对称，只需存储上三角部分，减少内存使用
  int NUMQNZ = 0;  // NUMQNZ = number of Q matrix nonzeros: Q矩阵非零元素总数
  for (int i = 0; i < segment_num; i++) {
    int NUMQ_blk = (traj_order + 1);              // NUMQ_blk = number of Q matrix block: 每段每维的控制点数
    NUMQNZ += 3 * NUMQ_blk * (NUMQ_blk + 1) / 2;  // 3维 × 上三角元素数
  }

  // MOSEK稀疏矩阵存储数组
  MSKint32t qsubi[NUMQNZ], qsubj[NUMQNZ];  // qsubi/qsubj = Q matrix subscript i/j: Q矩阵的行索引和列索引
  double qval[NUMQNZ];                     // qval = Q matrix values: Q矩阵对应位置的数值

  // ===== 构建二次目标矩阵Q =====
  {
    int sub_shift = 0;  // sub_shift = subscript shift: 当前段控制点在全局数组中的偏移
    int idx = 0;        // idx = index: Q矩阵三元组数组的当前索引

    // 遍历每个轨迹段
    for (int k = 0; k < segment_num; k++) {
      double scale_k = corridor[k].t;  // 当前段的时间缩放因子

      // 遍历每个空间维度 (x, y, z)
      for (int p = 0; p < 3; p++)
        // 遍历控制点对 (i,j)，构建Q矩阵的每个元素
        for (int i = 0; i < s1d1CtrlP_num; i++)    // 行控制点索引
          for (int j = 0; j < s1d1CtrlP_num; j++)  // 列控制点索引
            if (i >= j)                            // 仅存储上三角 (Q对称)
            {
              // 计算全局变量索引: 段偏移 + 维度偏移 + 控制点索引
              qsubi[idx] = sub_shift + p * s1d1CtrlP_num + i;  // 全局行索引
              qsubj[idx] = sub_shift + p * s1d1CtrlP_num + j;  // 全局列索引

              // ===== 时间缩放的Q矩阵元素计算 =====
              // 基本原理: ∫[derivative^2]dt 需要按时间段长度缩放
              // 时间缩放因子: T^(2*minimize_order - 3)
              if (min_order_l == min_order_u) {
                // 整数阶情况: 直接应用缩放
                qval[idx] = MQM(i, j) / pow(scale_k, 2 * min_order_u - 3);
              } else {
                // 分数阶情况: 在相邻整数阶间线性插值
                // weight_u = (minimize_order - floor) / T^(2*ceil-3)
                // weight_l = (ceil - minimize_order) / T^(2*floor-3)
                // Q = (weight_u + weight_l) * MQM
                qval[idx] = ((minimize_order - min_order_l) /
                                 pow(scale_k, 2 * min_order_u - 3) +
                             (min_order_u - minimize_order) /
                                 pow(scale_k, 2 * min_order_l - 3)) *
                            MQM(i, j);
              }
              idx++;
            }

      sub_shift += s1CtrlP_num;  // sub_shift = subscript shift: 移动到下一段控制点
    }
  }

  ros::Time time_end1 = ros::Time::now();  // 记录Q矩阵构建完成时间

  // ===== 将Q矩阵传递给MOSEK优化器 =====
  if (r == MSK_RES_OK)
    r = MSK_putqobj(task, NUMQNZ, qsubi, qsubj, qval);  // 设置二次目标函数

  // ===== 设置目标函数方向为最小化 =====
  if (r == MSK_RES_OK)
    r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);

  /**
   * ===================================================================
   * 第6步: MOSEK求解器执行和结果处理
   * ===================================================================
   *
   * 本节调用MOSEK内点法求解器求解二次规划问题，并处理各种可能的求解状态。
   *
   * MOSEK内点法求解过程:
   * 1. 初始化: 寻找初始可行点 (满足约束的起始点)
   * 2. 迭代优化: 沿着中心路径向最优解移动
   * 3. 收敛检查: 检查KKT条件和对偶间隙
   * 4. 终止: 达到最优解或检测到不可行
   *
   * 可能的求解状态:
   * - OPTIMAL: 找到全局最优解 ✓
   * - NEAR_OPTIMAL: 找到近似最优解 ✓
   * - INFEASIBLE: 约束冲突，无可行解 ✗
   * - UNBOUNDED: 目标函数无下界 ✗
   * - UNKNOWN: 数值问题或迭代限制 ⚠
   */

  bool solve_ok = false;  // 求解成功标志
  if (r == MSK_RES_OK) {
    // ===== 调用MOSEK内点法求解器 =====
    MSKrescodee trmcode;                        // trmcode = termination code: 终止状态码
    r = MSK_optimizetrm(task, &trmcode);        // 执行优化 (trm = terminate)
    MSK_solutionsummary(task, MSK_STREAM_LOG);  // 输出求解摘要到控制台

    if (r == MSK_RES_OK) {
      // ===== 检查求解状态并提取结果 =====
      MSKsolstae solsta;                          // solsta = solution status: 解的状态枚举
      MSK_getsolsta(task, MSK_SOL_ITR, &solsta);  // 获取内点解状态

      switch (solsta) {
        case MSK_SOL_STA_OPTIMAL:       // 最优解: 满足所有KKT条件
        case MSK_SOL_STA_NEAR_OPTIMAL:  // 近似最优解: 在容差范围内
        {
          // ===== 提取最优控制点值 =====
          r = MSK_getxx(task,
                        MSK_SOL_ITR,  // 请求内点解
                        x_var);       // 存储到x_var数组

          // ===== 提取最优目标函数值 =====
          r = MSK_getprimalobj(task, MSK_SOL_ITR, &primalobj);

          obj = primalobj;  // 返回给调用者
          solve_ok = true;  // 标记求解成功

          break;
        }

        // ===== 处理不可行情况 =====
        case MSK_SOL_STA_DUAL_INFEAS_CER:       // 对偶不可行证书
        case MSK_SOL_STA_PRIM_INFEAS_CER:       // 原始不可行证书
        case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:  // 近似对偶不可行
        case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  // 近似原始不可行
          printf("Primal or dual infeasibility certificate found.\n");
          printf("检测到不可行问题: 约束条件相互冲突或走廊过窄\n");
          break;

        case MSK_SOL_STA_UNKNOWN:  // 未知状态
          printf("The status of the solution could not be determined.\n");
          printf("求解状态未知: 可能由数值问题或迭代限制造成\n");
          // solve_ok = true; // debug模式可能接受未知解
          break;

        default:
          printf("Other solution status: %d\n", solsta);
          break;
      }
    } else {
      printf("Error while optimizing. Return code: %d\n", r);
    }
  }

  // ===== MOSEK错误处理 =====
  if (r != MSK_RES_OK) {
    char symname[MSK_MAX_STR_LEN];  // 错误符号名称
    char desc[MSK_MAX_STR_LEN];     // 错误描述信息

    printf("An error occurred while optimizing.\n");
    MSK_getcodedesc(r, symname, desc);  // 获取详细错误信息
    printf("MOSEK Error %s - '%s'\n", symname, desc);
  }

  // ===== 清理MOSEK资源 =====
  MSK_deletetask(&task);  // 删除优化任务，释放内存
  MSK_deleteenv(&env);    // 删除MOSEK环境，释放许可证

  // ===== 优化性能日志 =====
  ros::Time time_end2 = ros::Time::now();
  ROS_WARN("MOSEK求解耗时: %f 秒", (time_end2 - time_end1).toSec());

  // ===== 求解结果验证 =====
  if (!solve_ok) {
    ROS_WARN("MOSEK求解失败: 无法找到可行的轨迹解");
    return -1;  // 返回失败状态
  }

  /**
   * ===================================================================
   * 第7步: 结果转换和输出格式化
   * ===================================================================
   *
   * 本节将MOSEK求解器的结果转换为上层调用者期望的格式：
   * 1. 从C数组格式转换为Eigen格式
   * 2. 将一维控制点数组重塑为二维系数矩阵
   * 3. 按照[段×维度×控制点]的格式组织数据
   *
   * 输出矩阵格式 (PolyCoeff):
   * 行: 轨迹段索引 (0 to segment_num-1)
   * 列: [x0,x1,...,xn, y0,y1,...,yn, z0,z1,...,zn]
   *
   * 这种格式便于后续贝塞尔轨迹评估和控制器使用。
   */

  // ===== 格式转换: C数组 → Eigen向量 =====
  VectorXd d_var(ctrlP_num);  // 创建Eigen向量存储解
  for (int i = 0; i < ctrlP_num; i++)
    d_var(i) = x_var[i];  // 复制MOSEK解到Eigen格式

  // ===== 重塑为系数矩阵格式 =====
  // 目标格式: [段数 × 3*(阶数+1)] 矩阵
  // 每行代表一个轨迹段的所有控制点: [x0,x1,...,xn, y0,y1,...,yn, z0,z1,...,zn]
  PolyCoeff = MatrixXd::Zero(segment_num, 3 * (traj_order + 1));

  int var_shift = 0;                     // var_shift = variable shift: 当前段在全局变量数组中的起始位置
  for (int i = 0; i < segment_num; i++)  // 遍历每个轨迹段
  {
    for (int j = 0; j < 3 * n_poly; j++)       // 遍历该段的所有系数(x,y,z维度)
      PolyCoeff(i, j) = d_var(j + var_shift);  // 复制系数到输出矩阵

    var_shift += 3 * n_poly;  // var_shift = variable shift: 移动到下一段的系数起始位置
  }

  // ===== 成功返回 =====
  return 1;  // 返回成功状态码，PolyCoeff包含优化后的贝塞尔系数矩阵
}