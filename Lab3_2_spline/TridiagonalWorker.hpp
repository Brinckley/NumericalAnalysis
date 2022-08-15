//
// Created by alex on 15.08.22.
//

#ifndef LAB3_2_SPLINE_TRIDIAGONALWORKER_HPP
#define LAB3_2_SPLINE_TRIDIAGONALWORKER_HPP

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

template<class T>
class TridiagonalWorker {
public:
    int n;
    vector<T> a;
    vector<T> b;
    vector<T> c;

    TridiagonalWorker(int n) {
        this->n = n;
        a = vector<T>(n, 0);
        b = vector<T>(n, 0);
        c = vector<T>(n, 0);
    }

    TridiagonalWorker(vector<T> a, vector<T> b, vector<T> c) {
        this->n = a.size();
        this->a = a;
        this->b = b;
        this->c = c;
    }

    void read() {
        cin >> b[0] >> c[0];
        for (int i = 1; i < n - 1; ++i) {
            cin >> a[i] >> b[i] >> c[i];
        }
        cin >> a[n - 1] >> b[n - 1];
    }

    vector<T> TridiagonalSolver(vector<T> d) {
        vector<T> P(n, 0), Q(n, 0);
        vector<T> x(n, 0);

        P[0] = - c[0] / b[0];
        Q[0] = d[0] / b[0];
        for(int i = 1; i < n; ++i) {
            P[i] = - c[i] / (b[i] + a[i] * P[i - 1]);
            Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1]);
        }
        x[n - 1] = Q[n - 1];
        for(int i = n - 2; i >= 0; --i) {
            x[i] = P[i] * x[i + 1] + Q[i];
        }
        return x;
    }



};


#endif //LAB3_2_SPLINE_TRIDIAGONALWORKER_HPP
