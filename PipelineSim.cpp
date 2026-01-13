#include <iostream>
#include "Physics.h"
#include "Solver.h"

int main() {
    // 1. 设置仿真参数
    double length = 500.0;
    int numNodes = 50;

    // 边界条件
    double inletP = 2.2;    // 2.2 MPa
    double outletP = 2.0;   // 2 MPa
    double inletT = 40.0;   // 入口温度 40C

    // --- API 初始化参数 ---
    double initial_Rho_Observed = 850.0; // kg/m3
    double initial_P_Observed = 3.5;     // MPa (这里本来就是MPa，保持不变)
    double initial_T_Observed = 15.0;    // C
    Physics::initProperties(initial_Rho_Observed, initial_P_Observed, initial_T_Observed);

    // --- 粘温方程初始化 ---
    Physics::initViscosity(20.0, 100.0, 50.0, 20.0);

    // 3. 创建求解器
    Solver solver(numNodes, length, inletP, outletP, inletT);

    // 4. 初始化网格
    solver.initialize();

    // 5. 运行瞬变仿真
    // 时间步长 0.5s, 总时长 400s
    solver.solveTransient(0.5, 400.0, "final_result.csv");

    return 0;
}