//
// Created by alex on 16.08.22.
//

#ifndef LAB3_3_SQUARES_SQUARES_HPP
#define LAB3_3_SQUARES_SQUARES_HPP

#include "matrix.hpp"

class squares {
public:
    squares(vector<double> &x, vector<double> &y, int p) {
        n = x.size();
        this->x = x;                 // all x values
        this->y = y;                 // all y values
        this->p = p;                 // m
    }

    void square_build() {
        matrix<double> lhs(n, p);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < p; ++j) {
                lhs[i][j] = pow(x[i], j);       // matrix of values of pows
            }
        }

        matrix<double> lhs_t = lhs.transpose();
        matrix<double> A = lhs_t * lhs;
        vector<double> R = lhs_t * y;

        vector<double> answer = matrix<double>::LU_solver_improved(A, R);
        coef = answer;
    }

    void printPoly() {
        if (coef.empty())
            return;
        for (int i = 0; i < p; ++i) {
            cout << coef[i] << " * " <<  "x^" << i;
            if (i != p - 1) {
                cout << " + ";
            } else {
                cout << endl;
            }
        }
    }

    double squares_result() {
        double PHI = 0;
        for (int i = 0; i < n; ++i) {
            double tmp = 0;
            for (int j = 0; j < p; ++j) {
                tmp += coef[j] * pow(x[i] , j);
            }
            PHI += pow(tmp - y[i], 2.0);         // calculating mistake function
        }

        return PHI;
    }

private:
    int n;
    vector<double> x;
    vector<double> y;
    int p;

    vector<double> coef;
};


#endif //LAB3_3_SQUARES_SQUARES_HPP
