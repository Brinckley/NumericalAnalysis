//
// Created by alex on 10.08.22.
//

#ifndef LAB3_1_LAGRANGE_NEWTON_NEWTON_LAGRANGE_WORKER_HPP
#define LAB3_1_LAGRANGE_NEWTON_NEWTON_LAGRANGE_WORKER_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include "polynom.hpp"

using namespace std;

class newton_lagrange_worker {
private:
    vector<double> x;
    vector<double> y;
    int n;

    vector<vector<double>> memo;
    vector<vector<bool>> flags;

    double f(int l, int r) {
        if (flags[l][r])
            return memo[l][r];

        flags[l][r] = true;
        double res;
        if (l + 1 == r) {
            res = (y[l] - y[r]) / (x[l] - x[r]);
        } else {
            res = (f(l, r - 1) - f(l + 1, r)) / (x[l] - x[r]);
        }
        return memo[l][r] = res;
    }
public:

    newton_lagrange_worker(const vector<double> &x, const vector<double> &y) {
        this->x = x;
        this->y = y;
        this->n = x.size();

        memo.resize(n, vector<double>(n));
        flags.resize(n, vector<bool>(n));
    };

    polynom lagrange_method() {
        polynom res(vector<double>({0}));
        for (size_t i = 0; i < n; ++i) {
            polynom li(vector<double>({1}));
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    continue;
                }
                polynom xij(vector<double>{-x[j], 1});
                li = li * xij;
                li = li / (x[i] - x[j]);
            }
            res = res + y[i] * li;
        }
        return res;
    }

    polynom newton_method() {
        polynom res(vector<double>({y[0]}));
        polynom li(vector<double>({-x[0], 1}));
        int r = 0;
        for (int i = 1; i < n; ++i) {
            res = res + f(0, ++r) * li;
            li = li * polynom(vector<double>({-x[i], 1}));
        }
        return res;
    }
};

#endif //LAB3_1_LAGRANGE_NEWTON_NEWTON_LAGRANGE_WORKER_HPP
