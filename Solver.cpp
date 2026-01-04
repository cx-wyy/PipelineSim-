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
//初始化函数
void Solver::initialize() {
    nodes.resize(N);
// 假设初始为静止状态，压力线性分布
    for (int i = 0; i < N; ++i) {
        nodes[i].x = i * dx;
        nodes[i].z = 0;

        double ratio = (double)i / ((double)N - 1.0);
        nodes[i].P = inletPressure - (inletPressure - outletPressure) * ratio;
        nodes[i].V = 0.0; //初始静止

        // 温度初始化：假设全线初始为地温或者入口温度，这里暂设为入口温度
        nodes[i].T = inletTemp;

        nodes[i].rho = Physics::getDensity(nodes[i].P, nodes[i].T);
        nodes[i].visc = Physics::getViscosity(nodes[i].T);

        nodes[i].P_old = nodes[i].P;
        nodes[i].V_old = nodes[i].V;
        nodes[i].T_old = nodes[i].T;
    }
}

void Solver::solveTransient(double dt, double totalTime, const std::string& filename) {
    double time = 0;
    int steps = (int)(totalTime / dt);

    std::cout << "Starting Simulation (Fully Coupled P-V-T)..." << std::endl;
    std::ofstream file(filename);
    if (!file.is_open()) return;

    file << "time,node_id,x,P(Pa),P(MPa),v,T(C),visc(cSt),rho\n";
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
                << "P_in: " << nodes[0].P / 1e6 << " MPa, "
                << "T_out: " << nodes[N - 1].T << " C, "
                << "V_mid: " << nodes[N / 2].V << " m/s" << std::endl;
        }

        if (step % 10 == 0) {
            for (int i = 0; i < N; ++i) {
                file << time << "," << i << "," << nodes[i].x << ","
                    << nodes[i].P << "," << nodes[i].P / 1.0e6 << ","
                    << nodes[i].V << "," << nodes[i].T << ","
                    << nodes[i].visc * 1e6 << "," // 输出 cSt
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

    // --- 边界条件 (保持不变) ---
    A[0][0] = 1.0; b[0] = inletPressure;
    A[1][2] = 1.0; b[1] = inletTemp; // 假设 Solver 类里有 inletTemp
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

        // 2. 线性化点的速度选取 (关键修改)
        // 为了实现全隐式效果，这里不应只用 V_old (上一时刻)，
        // 而应该用当前节点中存储的 V (它是上一轮 Picard 迭代的结果)。
        // 第一轮迭代时，nodes[i].V 等于 nodes[i].V_old。
        double V_k_L = nL.V;
        double V_k_R = nR.V;
        double V_avg_k = (V_k_L + V_k_R) / 2.0; // 对应公式中的 V_average
        double abs_V_k = std::abs(V_avg_k);

        // 3. 计算摩擦系数 f (基于当前迭代速度)
        double visc_kin = Physics::getViscosity((nL.T + nR.T) / 2.0);
        double visc_dyn = visc_kin * rho_avg;
        double f = Physics::getFrictionFactor(V_avg_k, D, rho_avg, visc_dyn);

        // ==========================================
        // V_new * |V_new| ≈ 2 * |V_k| * V_new - |V_k| * V_k
        // ==========================================

        // 放入矩阵 A 的系数 (V_new 的系数): f/2D * 2|V_k|
        double Coeff_A_Mom = (f * abs_V_k) / D;

        // 放入向量 b 的项 (常数项移项到右边): - f/2D * (-|V_k|*V_k)
        double Val_b_Mom = (f * V_avg_k * abs_V_k) / (2.0 * D);

        // ==========================================
        // 能量方程：三次项
        // - f / (2 Cp d) * |V_new|^3
        // 这里使用 "半隐式" 处理：直接使用最新迭代的 V_k 计算值放入 RHS
        // 随着 Picard 迭代收敛，V_k -> V_new，从而满足方程
        // ==========================================
        double V3 = std::pow(abs_V_k, 3.0);
        double HeatSource_Fric = (f * V3) / (2.0 * D * Physics::Cp);

        // 热传导系数
        double Coeff_Heat = (4.0 * Physics::K_thermal) / (rho_avg * Physics::Cp * D);

        // --- 索引计算 ---
        int idx_Pi = 3 * i, idx_Vi = 3 * i + 1, idx_Ti = 3 * i + 2;
        int idx_Pip1 = 3 * (i + 1), idx_Vip1 = 3 * (i + 1) + 1, idx_Tip1 = 3 * (i + 1) + 2;

        // -----------------------------
        // A. 连续性方程 
        // -----------------------------
        double C1 = 1.0 / (2.0 * dt);
        double C2 = (rho_avg * a2 * theta) / dx;

        A[current_row][idx_Pi] = C1;     A[current_row][idx_Pip1] = C1;
        A[current_row][idx_Vi] = -C2;    A[current_row][idx_Vip1] = C2;
        b[current_row] = C1 * (nL.P_old + nR.P_old) - (rho_avg * a2 * (1.0 - theta) / dx) * (nR.V_old - nL.V_old);
        current_row++;

        // -----------------------------
        // B. 运动方程 
        // -----------------------------
        double M1 = 1.0 / (2.0 * dt);
        double M2 = theta / (rho_avg * dx);

        // V 的系数：时间项 + 摩擦项(泰勒展开线性部分)
        A[current_row][idx_Vi] = M1 + Coeff_A_Mom;
        A[current_row][idx_Vip1] = M1 + Coeff_A_Mom;

        // P 的系数
        A[current_row][idx_Pi] = -M2;
        A[current_row][idx_Pip1] = M2;

        double gravity = Physics::g * std::sin(0.0);
        double old_grad_P = ((1.0 - theta) / (rho_avg * dx)) * (nR.P_old - nL.P_old);

        // RHS: 旧时间步项 + 压力梯度旧值项 + 重力 + 摩擦常数项(泰勒展开常数部分)
        b[current_row] = M1 * (nL.V_old + nR.V_old) - old_grad_P - gravity + Val_b_Mom;
        current_row++;

        // -----------------------------
        // C. 能量方程 
        // -----------------------------
        double E_dt = 1.0 / (2.0 * dt);
        double E_dx = (V_avg_k * theta) / dx; //对流项系数也用当前迭代速度 V_avg_k
        double E_env = Coeff_Heat * 0.5;

        A[current_row][idx_Ti] = E_dt - E_dx + E_env;
        A[current_row][idx_Tip1] = E_dt + E_dx + E_env;

        double term_time_old = E_dt * (nL.T_old + nR.T_old);
        // 对流项旧值 (注意：旧时刻对流必须用旧时刻速度 V_old)
        double V_old_avg = (nL.V_old + nR.V_old) / 2.0;
        double term_conv_old = (V_old_avg * (1.0 - theta) / dx) * (nR.T_old - nL.T_old);
        double term_env_const = Coeff_Heat * Physics::T_ground;

        // RHS 加入 HeatSource_Fric (三次项)
        b[current_row] = term_time_old - term_conv_old + HeatSource_Fric + term_env_const;
        current_row++;
    }

    std::vector<double> result(dim);
    gaussianElimination(A, b, result);

    // 将计算结果更新到 nodes 中
    // 注意：如果是 Picard 迭代的第 1、2、3 次，这里更新的 V 会被下一轮循环
    // 用作 "V_k" 来重新计算泰勒系数，从而逼近真实解。
    for (int i = 0; i < N; ++i) {
        nodes[i].P = result[3 * i];
        nodes[i].V = result[3 * i + 1];
        nodes[i].T = result[3 * i + 2];
        nodes[i].rho = Physics::getDensity(nodes[i].P, nodes[i].T);
    }
}

void Solver::gaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x) {
    // 高斯消元
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