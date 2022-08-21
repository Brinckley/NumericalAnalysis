#include <iostream>

#include "Diff.hpp"


int main() {
    double x_ = 2;
    vector<double> x = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    vector<double> y = { 0.0, 2.0, 3.4142, 4.7321, 6.0 };
//    double x_ = 1;
//    vector<double> x = { -1.0, 0.0, 1.0, 2.0, 3.0 };
//    vector<double> y = { -0.5, 0.0, 0.50, 0.86603, 1.0 };
    Diff d = Diff(x, y);
    cout << "First diff result : " << d.diff_first(x_) << endl;
    cout << "Second diff result : " << d.diff_second(x_) << endl;
}
