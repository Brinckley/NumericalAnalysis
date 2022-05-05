#ifndef LAB2_2_ITERATIONS_AND_NEWTON_MATRIX_HPP
#define LAB2_2_ITERATIONS_AND_NEWTON_MATRIX_HPP

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

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

    T det2() {
        if(n == m && n == 2)
            return _matrix[0][0] * _matrix[1][1] - _matrix[0][1] * _matrix[1][0];
        return 0;
    }

private:
    vector<vector<T>> _matrix;
};


#endif //LAB2_2_ITERATIONS_AND_NEWTON_MATRIX_HPP
