//
// Created by alex on 16.08.22.
//

#ifndef LAB3_3_SQUARES_MATRIX_HPP
#define LAB3_3_SQUARES_MATRIX_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

template <class T>
class matrix {
public:
    int n;
    int m;

    vector<T>& operator[](int i) {
        return matrix_[i];
    }

    matrix(){n = 3; m = 3; matrix_ = vector<vector<T>>(3, vector<T>(3));}

    matrix<T>(int nt) {
        n = nt;
        m = nt;
        matrix_ = vector<vector<T>>(nt, vector<T>(nt));
        initializeZero();
    }

    matrix<T>(int nt, int mt) {
        n = nt;
        m = mt;
        matrix_ = vector<vector<T>>(nt, vector<T>(mt));
        initializeZero();
    }

    matrix<T>(const matrix<double> &matrix1) {
        n = matrix1.n;
        m = matrix1.m;
        matrix_ = vector<vector<T>>(n, vector<T>(m));
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix_[i][j] = matrix1.matrix_[i][j];
            }
        }
    }

    void initializeOne() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix_[i][j] = 0;
            }
            matrix_[i][i] = 1;
        }
    }

    void initializeZero() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix_[i][j] = 0;
            }
        }
    }

    void read() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                cin >> matrix_[i][j];
            }
        }
    }

    void print() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                cout << matrix_[i][j] << " ";
            }
            cout << "\n";
        }
    }

    friend matrix<T> operator*(const matrix<T> &left, const matrix<T> &right) {
        matrix<T> res(left.n, right.m);
        if(left.m != right.n)
            return res;
        for (int i = 0; i < left.n; ++i) {
            for (int j = 0; j < right.m; ++j) {
                for (int k = 0; k < left.m; ++k) {
                    res.matrix_[i][j] += left.matrix_[i][k] * right.matrix_[k][j];
                }
            }
        }
        return res;
    }

    friend bool operator==(const matrix<T> &left, const matrix<T> &right) {
        if(left.n != right.n || right.m != left.m)
            return false;

        for (int i = 0; i < left.n; ++i) {
            for (int j = 0; j < left.m; ++j) {
                if(left.matrix_[i][j] != right.matrix_[i][j])
                    return false;
            }
        }
        return true;
    }

    friend vector<T> operator*(const matrix<T> &left, const vector<T> &right) {
        vector<T> res(left.n);
        if (left.m != right.size())
            return res;
        for (size_t i = 0; i < left.n; ++i)
            for (size_t j = 0; j < left.m; ++j)
                res[i] += left.matrix_[i][j] * right[j];
        return res;
    }

    void swapR(int f, int s) {
        for(int i = 0; i < n; ++i)
            swap(matrix_[f][i], matrix_[s][i]);
    }

    void swapC(int f, int s) {
        for(int i = 0; i < n; ++i)
            swap(matrix_[i][f], matrix_[i][s]);
    }

    matrix<double> transpose() {
        matrix<double> result(m, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result[j][i] = matrix_[i][j];
            }
        }
        return result;
    }

    static void LU_setter(matrix<double> &A, matrix<double> &L, matrix<double> &U, matrix<double> &Permutation, int n, int &det) {
        L.initializeOne();
        Permutation.initializeOne();
        U = A;
        for(int i = 0; i < n; ++i) {
            double maxIndex = i;

            for(int j = i + 1; j < n; ++j) {
                if(abs(U[j][i]) > abs(U[maxIndex][i])) {
                    maxIndex = j;
                }
            }

            // setting row with the max element in column to the highest possible position
            if (maxIndex != i) {
                U.swapR(i, maxIndex);
                L.swapR(i, maxIndex);
                L.swapC(i, maxIndex);
                Permutation.swapC(i, maxIndex);
                det++;
            }

            for(int j = i + 1; j < n; ++j) {
                L[j][i] = U[j][i] / U[i][i]; // - (-a21/a11)
                for(int k = 0; k < n; ++k)
                    U[j][k] -= L[j][i] * U[i][k]; // working with the other elements in the row
            }
        }
    }

    static vector<double> LU_solver_improved(matrix<double> &A, vector<double> &R) {
        int n = R.size();
        int junk = 0;
        matrix<double> U(A);
        matrix<double> L(n);
        L.initializeZero();

        return LU_solver(A, R, L, U, n, junk);
    }

    static vector<double> LU_solver(matrix<double> &A, vector<double> &R, matrix<double> &L, matrix<double> &U, int n, int &det) {
        vector<double> b(n, 0);
        vector<double> y(n, 0);
        vector<double> x(n, 0);
        matrix<double> Permutation(n);
        det = 0;

        LU_setter(A, L, U, Permutation, n, det);

//        A.print();
//        cout << endl;
//        cout << "U___________" << endl;
//        U.print();
//        cout << endl;
//        cout << "L___________" << endl;
//        L.print();
//        cout << endl;
//        (L * U).print();

        if(det % 2 == 0)
            det = 1;
        else
            det = -1;

        for (int i = 0; i < n; ++i)
            det *= U[i][i];

        // applying permutations
        for(int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                b[i] += Permutation[j][i] * R[j];
            }
        }

        // Ly = b
        for(int i = 0; i < n; ++i) {
            double sum = 0;
            for(int k = 0; k < i; ++k) {
                sum += L[i][k] * y[k];
            }
            y[i] = b[i] - sum;
        }

        // Ux = y
        for(int i = n - 1; i >= 0; --i) {
            double sum = 0;
            for(int k = i + 1; k < n; ++k) {
                sum += U[i][k] * x[k];
            }
            x[i] = (y[i] - sum) / U[i][i];
        }

        return x;
    }

    static matrix<double> invertible(matrix<double> &A, vector<double> &R, matrix<double> &L, matrix<double> &U, int n) {
        matrix<double> result(n);
        vector<double> solutions(n, 0);
        int junk = 0;

        for(int i = 0; i < n; ++i) {
            vector<double> y(n, 0);
            y[i] = 1;
            vector<double> itr = LU_solver(A, y, L, U, n, junk);
            for(int k = 0; k < n; ++k)
                result[k][i] = itr[k];
        }

        return result;
    }
private:
    vector<vector<T>> matrix_;
};

#endif //LAB3_3_SQUARES_MATRIX_HPP
