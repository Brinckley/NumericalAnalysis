#include "newton_lagrange_worker.hpp"

double function(double x) {
    return log(x) + x;
}

int main() {
    double _x = 0.8;

    vector<double> x {0.1, 0.5, 0.9, 1.3};
    int n = x.size();
    vector<double> y (n);
    for (int i = 0; i < n; ++i) {
        y[i] = function(x[i]);
    }
    newton_lagrange_worker lagrange1(x, y);
    polynom lagrange = lagrange1.lagrange_method();
    cout << "Largange poly 1 : " << lagrange << endl;
    cout << "Largange delta 1 : " << abs(lagrange.calculate(_x) - function(_x)) << endl;

    newton_lagrange_worker newton1(x, y);
    polynom newton = newton1.newton_method();
    cout << "Newton poly 1 : " << newton << endl;
    cout << "Newton's delta 1 : " << abs(newton.calculate(_x) - atan(_x)) << endl;


    x = vector<double>{0.1, 0.5, 1.1, 1.3};
    n = x.size();
    y = vector<double>(n);
    for (int i = 0; i < n; ++i) {
        y[i] = function(x[i]);
    }
    newton_lagrange_worker lagrange2(x, y);
    lagrange = lagrange2.lagrange_method();
    cout << "Largange poly 2 : " << lagrange << endl;
    cout << "Delta 2 : " << abs(lagrange.calculate(_x) - function(_x)) << endl;

    newton_lagrange_worker newton2(x, y);
    newton = newton2.newton_method();
    cout << "Newton poly 2 : " << newton << endl;
    cout << "Newton's delta 2 : " << abs(newton.calculate(_x) - atan(_x)) << endl;

}
