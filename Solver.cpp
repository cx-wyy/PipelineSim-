#include "Solver.h"
#include "Physics.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

Solver::Solver(int numNodes, double length, double inletP, double outletP, double inletT)
    : N(numNodes), L(length), inletPressure(inletP), outletPressure(outletP), inletTemp(inletT) {
    dx = L / ((double)N - 1.0);
}

void Solver::initialize() {
    nodes.resize(N);
    // 假设初始为静止状态，压力线性分布
    for (int i = 0; i < N; ++i) {
        nodes[i].x = i * dx;
        nodes[i].z = 0;

        // 1. 压力初始化：线性下降
        double ratio = (double)i / ((double)N - 1.0);
        nodes[i].P = inletPressure - (inletPressure - outletPressure) * ratio;

        // 2. 速度初始化：保持你要求的 0 (静止启动)
        nodes[i].V = 0.0;

        // 3. 温度初始化：关键修改
        // 只有入口节点(i=0)是热的(40度)，管道内部其余部分是冷的(15度)
        if (i == 0) {
            nodes[i].T = inletTemp; // 40.0 C
        }
        else {
            nodes[i].T = 15.0;      // 管道内初始充填的是冷油
        }

        // 4. 根据当前的 P 和 T 计算初始 rho 和 viscosity
        nodes[i].rho = Physics::getDensity(nodes[i].P, nodes[i].T);
        nodes[i].visc = Physics::getViscosity(nodes[i].T);

        // 5. 记录旧值
        nodes[i].P_old = nodes[i].P;
        nodes[i].V_old = nodes[i].V;
        nodes[i].T_old = nodes[i].T;
    }
}

void Solver::solveTransient(double dt, double totalTime, const std::string& filename) {
    double time = 0;
    int steps = (int)(totalTime / dt);

    std::cout << "Starting Simulation (Fully Coupled P(MPa)-V-T)..." << std::endl;
    std::ofstream file(filename);
    if (!file.is_open()) return;

    // 输出表头
    file << "time,node_id,x,P(MPa),v,T(C),visc(cSt),rho\n";
    file << std::fixed << std::setprecision(6);

    for (int step = 0; step < steps; ++step) {
        time += dt;

        // 1. 保存旧状态
        for (int i = 0; i < N; ++i) {
            nodes[i].P_old = nodes[i].P;
            nodes[i].V_old = nodes[i].V;
            nodes[i].T_old = nodes[i].T;
        }

        // 2. Picard 迭代
        for (int iter = 0; iter < 4; ++iter) {
            solveStep(dt);
        }

        // 3. 输出
        if (step % 50 == 0) {
            std::cout << "Time: " << std::setprecision(2) << time << "s | "
                << "P_in: " << nodes[0].P << " MPa, "
                << "T_out: " << nodes[N - 1].T << " C, "
                << "V_mid: " << nodes[N / 2].V << " m/s" << std::endl;
        }

        if (step % 10 == 0) {
            for (int i = 0; i < N; ++i) {
                file << time << "," << i << "," << nodes[i].x << ","
                    << nodes[i].P << "," // 直接输出 MPa
                    << nodes[i].V << "," << nodes[i].T << ","
                    << nodes[i].visc * 1e6 << ","
                    << nodes[i].rho << "\n";
            }
        }
    }
    file.close();
    std::cout << "Simulation Complete." << std::endl;
}

void Solver::solveStep(double dt) {
    int dim = 3 * N;
    std::vector<std::vector<double>> A(dim, std::vector<double>(dim, 0.0));
    std::vector<double> b(dim, 0.0);

    // --- 边界条件 ---
    // P 单位为 MPa
    A[0][0] = 1.0; b[0] = inletPressure;
    A[1][2] = 1.0; b[1] = inletTemp;
    A[dim - 1][3 * (N - 1)] = 1.0; b[dim - 1] = outletPressure;

    int current_row = 2;

    for (int i = 0; i < N - 1; ++i) {
        Node& nL = nodes[i];
        Node& nR = nodes[i + 1];

        // 1. 物理参数平均
        double rho_avg = (nL.rho + nR.rho) / 2.0;
        double a_avg = Physics::getWaveSpeed(rho_avg);
        double a2 = a_avg * a_avg;
        double D = Physics::Diameter;

        // 2. 线性化点的速度选取
        double V_k_L = nL.V;
        double V_k_R = nR.V;
        double V_avg_k = (V_k_L + V_k_R) / 2.0;
        double abs_V_k = std::abs(V_avg_k);

        // 3. 计算摩擦系数
        double visc_kin = Physics::getViscosity((nL.T + nR.T) / 2.0);
        double visc_dyn = visc_kin * rho_avg;
        double f = Physics::getFrictionFactor(V_avg_k, D, rho_avg, visc_dyn);

        // V 系数 (动量方程)
        double Coeff_A_Mom = (f * abs_V_k) / D;
        double Val_b_Mom = (f * V_avg_k * abs_V_k) / (2.0 * D);

        // 能量方程源项
        double V3 = std::pow(abs_V_k, 3.0);
        double HeatSource_Fric = (f * V3) / (2.0 * D * Physics::Cp);
        double Coeff_Heat = (4.0 * Physics::K_thermal) / (rho_avg * Physics::Cp * D);

        // --- 索引计算 ---
        int idx_Pi = 3 * i, idx_Vi = 3 * i + 1, idx_Ti = 3 * i + 2;
        int idx_Pip1 = 3 * (i + 1), idx_Vip1 = 3 * (i + 1) + 1, idx_Tip1 = 3 * (i + 1) + 2;

        // -----------------------------
        // A. 连续性方程 (含 theta)
        // -----------------------------
        double C1 = 1.0 / (2.0 * dt);
        double C2 = (rho_avg * a2 * theta) / (dx * 1.0e6);

        A[current_row][idx_Pi] = C1;     A[current_row][idx_Pip1] = C1;
        A[current_row][idx_Vi] = -C2;    A[current_row][idx_Vip1] = C2;

        double C2_RHS = (rho_avg * a2 * (1.0 - theta)) / (dx * 1.0e6);
        b[current_row] = C1 * (nL.P_old + nR.P_old) - C2_RHS * (nR.V_old - nL.V_old);
        current_row++;

        // -----------------------------
        // B. 运动方程 (含 theta)
        // -----------------------------
        double M1 = 1.0 / (2.0 * dt);
        double M2 = (theta * 1.0e6) / (rho_avg * dx);

        A[current_row][idx_Vi] = M1 + Coeff_A_Mom;
        A[current_row][idx_Vip1] = M1 + Coeff_A_Mom;
        A[current_row][idx_Pi] = -M2;
        A[current_row][idx_Pip1] = M2;

        double gravity = Physics::g * std::sin(0.0);
        double old_grad_P = ((1.0 - theta) * 1.0e6 / (rho_avg * dx)) * (nR.P_old - nL.P_old);

        b[current_row] = M1 * (nL.V_old + nR.V_old) - old_grad_P - gravity + Val_b_Mom;
        current_row++;

        // -----------------------------
        // C. 能量方程 (为全隐式)
        // -----------------------------
        // 1. 时间项系数: (T_new - T_old) / dt
        double E_dt = 1.0 / dt;

        // 2. 对流项系数: V/dx
        // V_avg_k 是流速，方向通常为正 (i -> i+1)
        double E_dx = V_avg_k / dx;

        // 3. 环境换热系数: 
        // 迎风格式，只计算节点 i+1 的换热
        double E_env = (4.0 * Physics::K_thermal) / (nR.rho * Physics::Cp * D);

        // --- 填充矩阵 ---
        //利用当前行 (current_row) 来求解节点 i+1 (nR) 的温度
        // 公式变型为: (1/dt + V/dx + E_env)*T(i+1) - (V/dx)*T(i) = T_old(i+1)/dt + Sources

        // 系数对应 T(i) (idx_Ti) -> 移项到左边是负的
        A[current_row][idx_Ti] = -E_dx;

        // 系数对应 T(i+1) (idx_Tip1) -> 主对角线部分
        A[current_row][idx_Tip1] = E_dt + E_dx + E_env;

        // --- 填充 RHS (b向量) ---
        // 时间项只取 i+1 的旧值
        double term_time_old = E_dt * nR.T_old;

        // 环境温度项 (常数部分)
        double term_env_const = E_env * Physics::T_ground;

        // 摩擦生热 (使用 i+1 点的流速估算或保持段平均)
        // 这里沿用段平均值以保持能量守恒的一致性
        b[current_row] = term_time_old + HeatSource_Fric + term_env_const;

        current_row++;
    }

    std::vector<double> result(dim);
    gaussianElimination(A, b, result);

    for (int i = 0; i < N; ++i) {
        nodes[i].P = result[3 * i];
        nodes[i].V = result[3 * i + 1];
        nodes[i].T = result[3 * i + 2];
        nodes[i].rho = Physics::getDensity(nodes[i].P, nodes[i].T);
    }
}

void Solver::gaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x) {
    int n = b.size();
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) maxRow = k;
        }
        std::swap(A[i], A[maxRow]);
        std::swap(b[i], b[maxRow]);
        for (int k = i + 1; k < n; ++k) {
            double c = -A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                if (i == j) A[k][j] = 0;
                else A[k][j] += c * A[i][j];
            }
            b[k] += c * b[i];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }
}