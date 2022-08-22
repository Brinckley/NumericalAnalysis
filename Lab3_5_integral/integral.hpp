//
// Created by alex on 20.08.22.
//

#ifndef LAB3_5_INTEGRAL_INTEGRAL_HPP
#define LAB3_5_INTEGRAL_INTEGRAL_HPP

#include "vector"
#include <cmath>

using namespace std;

class integral_worker {
public:
    integral_worker(double x0, double xk, double h) {
        this->x0 = x0;
        this->xk = xk;
        this->h = h;
    }

    double rectangle_method() {
        double xi0 = x0;
        double xi1 = x0 + h;
        double res = 0;

        do {
            res += h * my_function((xi0 + xi1) / 2);
            xi0 += h;
            xi1 += h;
        } while (xi0 < xk);

        return res;
    }

    double trapezoid_method() {
        double xi0 = x0;
        double xi1 = x0 + h;
        double res = 0;

        do {
            res += h * (my_function(xi1) + my_function(xi0));
            xi0 += h;
            xi1 += h;
        } while (xi0 < xk);

        return res / 2;
    }

    double simpson_method() {
        double xi0;
        double res = my_function(x0) + my_function(xk);

        xi0 = x0 + h / 2;
        do {
            res += 4 * my_function(xi0 / 2);
            xi0 += h;
        } while (xi0 < xk);

        xi0 = x0 + h;
        do {
            res += 2 * my_function(xi0);
            xi0 += h;
        } while (xi0 < xk);

        return res * h / 3;
    }

    static double Runge_Romberg_Richardson_method(double Fh, double Fkh, double k, double p) {
        return (Fh - Fkh) / (pow(k, p) - 1);
    }
private:
    double x0;
    double xk;
    double h;

    double my_function(double x) {
        return pow(x, 2) / (pow(x, 4) + 256);
    }
};


#endif //LAB3_5_INTEGRAL_INTEGRAL_HPP
