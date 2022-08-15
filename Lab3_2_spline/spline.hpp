//
// Created by alex on 15.08.22.
//

#ifndef LAB3_2_SPLINE_SPLINE_HPP
#define LAB3_2_SPLINE_SPLINE_HPP

#include "polynom.hpp"
#include "TridiagonalWorker.hpp"


class spline {
public:
    spline(vector<double> &x, vector<double> &y) {
        this->x = x;
        this->y = y;
        n = x.size();
        a.resize(n);
        b.resize(n);
        c.resize(n);
        d.resize(n);
    }

    void builder() {
        vector<double> h(n);       // h = xi - xi-1
        for (int i = 1; i < n; ++i) {
            h[i] = x[i] - x[i - 1];
        }

        // solving for c
        vector<double> ts_a(n - 2);
        vector<double> ts_b(n - 2);
        vector<double> ts_c(n - 2);
        vector<double> ts_d(n - 2);

        for (int i = 2; i < n; ++i) {     // for i 2...n - 1
            ts_a[i - 2] = h[i - 1];
            ts_b[i - 2] = 2.0 * (h[i - 1] + h[i]);
            ts_c[i - 2] = h[i];
            ts_d[i - 2] = 3.0 * ((y[i] - y[i - 1]) / h[i] - (y[i - 1] - y[i - 2]) / h[i - 1]);
        }
        ts_a[0] = 0;                       // for i = 1
        ts_c.back() = 0;                   // last and first in system are zeroes
        TridiagonalWorker<double> tridiagonalWorker(ts_a, ts_b, ts_c);         // system like 3.13

        vector<double> solution_c = tridiagonalWorker.TridiagonalSolver(ts_d);
        for (int i = 2; i < n; ++i) {
            c[i] = solution_c[i - 2];
            a[i] = y[i - 1];
        }
        for (int i = 1; i < n; ++i) {
            b[i] = (y[i] - y[i - 1]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        }

        // last condition
        c[1] = 0.0;
        b[n - 1] = (y[n - 1] - y[n - 2]) / h[n - 1] - (2.0 / 3.0) * h[n - 1] * c[n - 1];
        d[n - 1] = - c[n - 1] / (3 * h[n - 1]);
    }

    double solver(double _x) {
        for (size_t i = 1; i < n; ++i) {
            if (x[i - 1] <= _x && _x <= x[i]) {
                double x1 = _x - x[i - 1];
                double x2 = x1 * x1;
                double x3 = x2 * x1;
                return a[i] + b[i] * x1 + c[i] * x2 + d[i] * x3;
            }
        }
        return INFINITY;
    }
private:
    vector<double> a, b, c, d;
    int n;

    vector<double> x;
    vector<double> y;


};


#endif //LAB3_2_SPLINE_SPLINE_HPP
