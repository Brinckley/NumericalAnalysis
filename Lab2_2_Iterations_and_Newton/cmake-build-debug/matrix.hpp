#ifndef LAB2_2_ITERATIONS_AND_NEWTON_MATRIX_HPP
#define LAB2_2_ITERATIONS_AND_NEWTON_MATRIX_HPP

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;


 vector<double> operator-(const vector<double> &left, const vector<double> &right) {
    vector<double> res(left.size());
     if(left.size() != right.size()) return res;
     for(int i = 0; i < left.size(); ++i) {
         res[i] = left[i] - right[i];
     }
    return res;
}

template <class T>
class matrix {
public:
    int n;
    int m;

    vector<T>& operator[](int i) {
        return _matrix[i];
    }

    matrix(){ n = 3; m = 3; _matrix = vector<vector<T>>(3, vector<T>(3));}

    matrix(int nt) {
        n = nt;
        m = nt;
        _matrix = vector<vector<T>>(nt, vector<T>(nt));
        initializeZero();
    }

    matrix(int nt, int mt) {
        n = nt;
        m = mt;
        _matrix = vector<vector<T>>(nt, vector<T>(mt));
        initializeZero();
    }

    matrix(int size, vector<T> v) {
        n = size;
        m = size;
        _matrix = vector<vector<T>>(n, vector<T>(n));
        int k = 0;
        for(int i = 0; i < n ;++i) {
            for(int j = 0; j < n; ++j) {
                _matrix[i][j] = v[k];
                ++k;
            }
        }
    }

    void initializeOne() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                _matrix[i][j] = 0;
            }
            _matrix[i][i] = 1;
        }
    }

    void initializeZero() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                _matrix[i][j] = 0;
            }
        }
    }

    void read() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                cin >> _matrix[i][j];
            }
        }
    }

    void print() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                cout << _matrix[i][j] << " ";
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
                    res._matrix[i][j] += left._matrix[i][k] * right._matrix[k][j];
                }
            }
        }
        return res;
    }

    friend vector<T> operator*(const matrix<T> &left, const vector<T> &right) {
        vector<T> res(right.size(), 0);
        if(left.m != right.size())
            return res;
        for (int i = 0; i < left.n; ++i) {
            for (int j = 0; j < right.size(); ++j) {
                res[i] += left._matrix[i][j] * right[j];
            }
        }
        return res;
    }

    friend bool operator==(const matrix<T> &left, const matrix<T> &right) {
        if(left.n != right.n || right.m != left.m)
            return false;

        for (int i = 0; i < left.n; ++i) {
            for (int j = 0; j < left.m; ++j) {
                if(left._matrix[i][j] != right._matrix[i][j])
                    return false;
            }
        }
        return true;
    }

    void swapR(int f, int s) {
        for(int i = 0; i < n; ++i)
            swap(_matrix[f][i], _matrix[s][i]);
    }

    void swapC(int f, int s) {
        for(int i = 0; i < n; ++i)
            swap(_matrix[i][f], _matrix[i][s]);
    }

    T det2() {
        if(n == m && n == 2)
            return _matrix[0][0] * _matrix[1][1] - _matrix[0][1] * _matrix[1][0];
        return 0;
    }

    T normaMatrixC() {
        T max = 0;
        for(int j = 0; j < n; ++j) {
            max += abs(_matrix[0][j]);
        }

        for(int i = 1; i < n; ++i){
            T sum = 0;
            for(int j = 0; j < n; ++j) {
                sum += abs(_matrix[i][j]);
            }
            if(sum > max)
                max = sum;
        }
        return max;
    }

    void LU_setter(matrix<T> &A, matrix<T> &L, matrix<T> &U, matrix<T> &Permutation, int n, int &det) {
        L.initializeOne();
        Permutation.initializeOne();
        U = A;

        for(int i = 0; i < n; ++i) {
            double maxInColumn = 0;
            double maxIndex = -1;

            for(int m = i; m < n; ++m) {
                if(abs(U[m][i]) > maxInColumn) {
                    maxIndex = m;
                    maxInColumn = abs(U[m][i]);
                }
            }

            if(maxInColumn == 0)
                continue; //empty column case

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

    vector<T> LU_solver(matrix<T> &A, vector<T> &R, matrix<T> &L, matrix<T>  &U, int n, int &det) {
        vector<T> b(n, 0);
        vector<T> y(n, 0);
        vector<T> x(n, 0);
        matrix<T> Permutation(n);
        det = 0;

        LU_setter(A, L, U, Permutation, n, det);

        if(det % 2 == 0)
            det = 1;
        else
            det = -1;

        for (int i = 0; i < n; ++i) {
            det *= U[i][i];
        }

        // applying permutations
        for(int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                b[i] += Permutation[j][i] * R[j];
            }
        }

        // Ly = b
        for(int i = 0; i < n; ++i) {
            T sum = 0;
            for(int k = 0; k < i; ++k) {
                sum += L[i][k] * y[k];
            }
            y[i] = b[i] - sum;
        }

        // Ux = y
        for(int i = n - 1; i >= 0; --i) {
            T sum = 0;
            for(int k = i + 1; k < n; ++k) {
                sum += U[i][k] * x[k];
            }
            x[i] = (y[i] - sum) / U[i][i];
        }

        return x;
    }

    matrix<T> invertible(matrix<T> &A, int n) {
        matrix<T> L(n);
        matrix<T> U(n);
        matrix<T> result(n);
        matrix<T> solutions(n, 0);
        int junk = 0;

        for(int i = 0; i < n; ++i) {
            vector<T> y(n, 0);
            y[i] = 1;
            vector<T> itr = LU_solver(A, y, L, U, n, junk);
            for(int k = 0; k < n; ++k)
                result[k][i] = itr[k];
        }

        return result;
    }

private:
    vector<vector<T>> _matrix;
};


#endif //LAB2_2_ITERATIONS_AND_NEWTON_MATRIX_HPP
