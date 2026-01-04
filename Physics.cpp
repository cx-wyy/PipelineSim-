#include "Physics.h"
#include <cmath>
#include <iostream>
#include <iomanip>

// 初始化静态成员
double Physics::rho_std_computed = 0.0;
double Physics::B_T_computed = 0.0;
double Physics::alpha_computed = 0.0;
double Physics::ASTM_A = 0.0;
double Physics::ASTM_B = 0.0;

void Physics::initProperties(double rho_init_kgm3, double P_init_Pa, double T_init_C) {
    std::cout << "--- API 11.1 Initialization Start ---" << std::endl;

    double T_F = T_init_C * 1.8 + 32.0;
    double P_psig = P_init_Pa * 145.038 * 1e-6;
    if (P_psig < 0) P_psig = 0;

    double A = -1.99470;
    double B = 0.00013427;
    double C = 793920.0;
    double D = 2326.0;

    double rho_60 = rho_init_kgm3;
    double tolerance = 1e-4;
    int max_iter = 100;

    double final_alpha_60 = 0.0;
    double final_Fp_psi = 0.0;

    for (int i = 0; i < max_iter; ++i) {
        double K0 = 0.0, K1 = 0.0, K2 = 0.0;
        if (rho_60 < 770.3520) {
            K0 = 192.4571; K1 = 0.2438; K2 = 0.0;
        }
        else if (rho_60 >= 770.3520 && rho_60 < 787.5195) {
            K0 = 1489.067; K1 = 0.0; K2 = -0.0018684;
        }
        else if (rho_60 >= 787.5195 && rho_60 < 838.3127) {
            K0 = 330.301; K1 = 0.0; K2 = 0.0;
        }
        else {
            K0 = 103.8720; K1 = 0.2701; K2 = 0.0;
        }

        double alpha_60_F = (K0 / (rho_60 * rho_60)) + (K1 / rho_60) + K2;
        double delta_T = T_F - 60.0;
        double exp_term_Ctl = -alpha_60_F * delta_T * (1.0 + 0.8 * alpha_60_F * delta_T);
        double Ctl = std::exp(exp_term_Ctl);

        double exp_term_Fp = A + B * T_F + (C + D * T_F) / (rho_60 * rho_60);
        double Fp_psi = 1.0e-5 * std::exp(exp_term_Fp);
        double Cpl = 1.0 / (1.0 - Fp_psi * P_psig);

        double rho_calc = rho_60 * Ctl * Cpl;
        double diff = rho_init_kgm3 - rho_calc;

        if (std::abs(diff) < tolerance) {
            final_alpha_60 = alpha_60_F;
            final_Fp_psi = Fp_psi;
            break;
        }
        rho_60 = rho_init_kgm3 / (Ctl * Cpl);
    }

    Physics::rho_std_computed = rho_60;
    Physics::alpha_computed = final_alpha_60 * 1.8;
    if (final_Fp_psi > 1e-15) {
        double B_T_psi = 1.0 / final_Fp_psi;
        Physics::B_T_computed = B_T_psi * 6894.76;
    }
    else {
        Physics::B_T_computed = 1.5e9;
    }

    std::cout << "--- API Initialization Complete ---" << std::endl;
    std::cout << "Locked Base Density (60F): " << Physics::rho_std_computed << " kg/m3" << std::endl;
}

void Physics::initViscosity(double T1_C, double v1_cSt, double T2_C, double v2_cSt) {
    // 转换温度为开尔文
    double T1_K = T1_C + T_Kelvin;
    double T2_K = T2_C + T_Kelvin;

    // Z = v + 0.7 (近似公式， v > 2 cSt)
    double Z1 = v1_cSt + 0.7;
    double Z2 = v2_cSt + 0.7;

    // 构建方程组:
    // lg(lg(Z1)) = A - B * lg(T1)
    // lg(lg(Z2)) = A - B * lg(T2)

    double Y1 = std::log10(std::log10(Z1));
    double Y2 = std::log10(std::log10(Z2));
    double X1 = std::log10(T1_K);
    double X2 = std::log10(T2_K);

    // 求解斜率 -B
    double slope = (Y2 - Y1) / (X2 - X1);
    ASTM_B = -slope;

    // 求解截距 A = Y + B*X
    ASTM_A = Y1 + ASTM_B * X1;

    std::cout << "--- Viscosity Parameters (ASTM D341) ---" << std::endl;
    std::cout << "A: " << ASTM_A << ", B: " << ASTM_B << std::endl;
}

double Physics::getViscosity(double T_C) {
    double T_K = T_C + T_Kelvin;

    // lg(lg(Z)) = A - B * lg(T)
    double log_T = std::log10(T_K);
    double log_log_Z = ASTM_A - ASTM_B * log_T;

    double log_Z = std::pow(10.0, log_log_Z);
    double Z = std::pow(10.0, log_Z);

    // Z = v + 0.7 => v = Z - 0.7 (cSt)
    double v_cSt = Z - 0.7;

    // 安全检查
    if (v_cSt < 1e-6) v_cSt = 1e-6;

    // 转换为 m^2/s (1 cSt = 1e-6 m^2/s)
    return v_cSt * 1.0e-6;
}

double Physics::getDensity(double P_Pa, double T_C) {
    double P_gauge_Pa = P_Pa - P_std_Pa;
    double term_T = alpha_computed * (T_C - T_std_C);
    double term_P = P_gauge_Pa / B_T_computed;
    double exponent = -term_T + term_P;
    return rho_std_computed * std::exp(exponent);
}

double Physics::getWaveSpeed(double rho) {
    double Ks = B_T_computed;
    double C1 = 0.85;
    double numerator = Ks / rho;
    double denominator = 1.0 + (Ks / E_modulus) * (Diameter / WallThickness) * C1;
    if (numerator < 0 || denominator < 0) return 1000.0;
    return std::sqrt(numerator / denominator);
}

double Physics::getFrictionFactor(double V, double D, double rho, double mu_visc) {
    // mu_visc 在这里输入是 动力粘度 Pa.s = rho * 运动粘度
    // 或者我们修改接口只传运动粘度? 为了兼容Colebrook公式，这里认为传入的是动力粘度
    // 如果传入的是动力粘度，Re = rho * V * D / mu

    if (std::abs(V) < 1e-6) return 0.02;
    double Re = (rho * std::abs(V) * D) / mu_visc;

    if (Re < 2300.0) return 64.0 / Re;

    double rel_rough = Roughness / D;
    double t1 = rel_rough / 3.7;
    double t2 = 5.74 / std::pow(Re, 0.9);
    double f = 0.25 / std::pow(std::log10(t1 + t2), 2);

    double tolerance = 1e-8;
    int max_iter = 50;
    for (int i = 0; i < max_iter; ++i) {
        double sqrt_f = std::sqrt(f);
        double term = rel_rough + 9.35 / (Re * sqrt_f);
        double rhs = 1.14 - 2.0 * std::log10(term);
        double f_new = 1.0 / (rhs * rhs);
        if (std::abs(f_new - f) < tolerance) return f_new;
        f = f_new;
    }
    return f;
}