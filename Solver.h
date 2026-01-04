#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <string>  
#include "Node.h"

class Solver {
public:
    // 增加 inletT 参数
    Solver(int numNodes, double length, double inletP, double outletP, double inletT);

    void initialize();
    void solveTransient(double dt, double totalTime, const std::string& filename);
    std::vector<Node> get_nodes() { return nodes; }

private:
    int N;
    double L;
    double dx;
    double inletPressure;
    double outletPressure;
    double inletTemp; // 新增：入口温度

    const double theta = 0.55; // Preissmann 权重

    std::vector<Node> nodes;

    // 构建并求解线性方程组 (Ax = b)，矩阵大小为 3N
    void solveStep(double dt);
    void gaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& result);
};

#endif // SOLVER_H