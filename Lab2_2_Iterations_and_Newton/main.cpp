#include <iostream>
#include <math.h>
#include "cmake-build-debug/system_2l2_worker.hpp"

// 2 * x1 - cos(x2) = 0
// 2 * x2 - exp(x1) = 0

//graphically selected area
// 0.3 < x1 < 0.5
// pi / 6 < x2 < pi / 4


int main() {
    double EPS = 0.00001;
    int iter = 0;
    //cin >> EPS;
    system_2l2_worker<double> solN(EPS);
    //pair<double, double> sols = solN.newtonSolverSimple(0.4, (M_PI / 5), iter);
    vector<double> sols = solN.newtonSolver(vector<double>{0.4, (M_PI / 5)}, iter);
    cout << "Newton method answers: x1 = "<< sols[0] << ", x2 = " << sols[1] << ". Iteration number: " << iter <<endl;

    //pair<double, double> sols = solN.newtonSolverSimple(0.4, (M_PI / 5), iter);
    pair<double, double> sols2 = solN.iterationSolver(0, 0.5, 0, M_PI / 2 , iter);
    cout << "Newton method answers: x1 = "<< sols[0] << ", x2 = " << sols[1] << ". Iteration number: " << iter <<endl;

}