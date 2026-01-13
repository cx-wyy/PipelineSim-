#ifndef NODE_H
#define NODE_H

struct Node {
    double x;       // 距离 (m)
    double z;       // 高程 (m)
    double P;       // 压力 (MPa)
    double V;       // 流速 (m/s)
    double T;       // 温度 (deg C)
    double rho;     // 密度 (kg/m^3)
    double visc;    // 运动粘度 (m^2/s)

    // 用于存储上一时间步的结果
    double P_old;
    double V_old;
    double T_old;

    Node() : x(0), z(0), P(0), V(0), T(15.0), rho(0), visc(0), P_old(0), V_old(0), T_old(15.0) {}
};

#endif // NODE_H