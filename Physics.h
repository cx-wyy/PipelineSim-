#ifndef PHYSICS_H
#define PHYSICS_H

#include <cmath>
#include <iostream>

class Physics {
public:
    // --- 物理常数 ---
	static constexpr double g = 9.81;//m/s^2 重力加速度
    static constexpr double P_std_Pa = 101325.0; // 标准大气压 Pa
    static constexpr double T_std_C = 15.56;     // 60F 对应的摄氏度
    static constexpr double T_Kelvin = 273.15;   // 摄氏度转开尔文

    // --- 管道参数 ---
	static constexpr double Diameter = 0.508;        // m 内径
    static constexpr double Roughness = 4.5e-5;      // m 绝对粗糙度
    static constexpr double WallThickness = 0.008;   // m 壁厚
    static constexpr double E_modulus = 2.06e11;     // pa 钢管弹性模量

    // --- 热力学参数 (新增) ---
    static constexpr double Cp = 2000.0;       // 比热容 J/(kg*C) (典型油品)
    static constexpr double K_thermal = 2.0;   // 总传热系数 W/(m^2*C)
    static constexpr double T_ground = 10.0;   // 地温 (C)

    // --- 状态方程参数 (由初始化计算得出并锁定) ---
    static double rho_std_computed;
    static double B_T_computed;
    static double alpha_computed;

    // --- 粘温方程参数 (ASTM D341) ---
    static double ASTM_A;
    static double ASTM_B;

    /**
     * @brief API 11.1 初始化核心函数
     */
    static void initProperties(double rho_init_kgm3, double P_init_MPa, double T_init_C);

    /**
     * @brief 初始化粘温方程 (ASTM D341)
     * 根据两点数据拟合 A 和 B
     * @param T1_C 温度1 (C)
     * @param v1_cSt 运动粘度1 (cSt = mm^2/s)
     * @param T2_C 温度2 (C)
     * @param v2_cSt 运动粘度2 (cSt = mm^2/s)
     */
    static void initViscosity(double T1_C, double v1_cSt, double T2_C, double v2_cSt);

    /**
     * @brief 获取运动粘度
     * @param T_C 温度 (C)
     * @return 运动粘度 (m^2/s) Note: 内部计算用cSt但返回SI单位
     */
    static double getViscosity(double T_C);

    static double getDensity(double P_Pa, double T_C);
    static double getWaveSpeed(double rho);
    static double getFrictionFactor(double V, double D, double rho, double mu_visc);
};

#endif // PHYSICS_H