#include <iostream>
#include <math.h>
#include "cmake-build-debug/system_2l2_worker.hpp"

// 2 * x1 - cos(x2) = 0
// 2 * x2 - exp(x1) = 0

//graphically selected area
// 0.3 < x1 < 0.5
// pi / 6 < x2 < pi / 4


int main() {
    double EPS = 0.001;
    int iter = 0;
    //cin >> EPS;
    system_2l2_worker<double> solN(EPS);
    pair<double, double> sols = solN.newtonSolver(0.4, (M_PI / 5), iter);
    cout << "Newton method answers: x1 = "<< sols.first << ", x2 = " << sols.second << ". Iteration number: " << iter <<endl;

}