#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <string>
#include "Node.h"

class Solver {
public:
    // P µ¥Î»Îª MPa
    Solver(int numNodes, double length, double inletP, double outletP, double inletT);

    void initialize();
    void solveTransient(double dt, double totalTime, const std::string& filename);
    std::vector<Node> get_nodes() { return nodes; }

private:
    int N;
    double L;
    double dx;
    double inletPressure; // MPa
    double outletPressure; // MPa
    double inletTemp;

    const double theta = 0.55;

    std::vector<Node> nodes;

    void solveStep(double dt);
    void gaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& result);
};

#endif // SOLVER_H