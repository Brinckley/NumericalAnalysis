#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

template <class T>
class Matrix {
public:
    int n;
    int m;

    vector<T>& operator[](int i) {
        return matrix[i];
    }

    Matrix(){n = 3; m = 3; matrix = vector<vector<T>>(3, vector<T>(3));}

    Matrix(int nt) {
        n = nt;
        m = nt;
        matrix = vector<vector<T>>(nt, vector<T>(nt));
        initializeZero();
    }

    Matrix(int nt, int mt) {
        n = nt;
        m = mt;
        matrix = vector<vector<T>>(nt, vector<T>(mt));
        initializeZero();
    }

    void initializeOne() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix[i][j] = 0;
            }
            matrix[i][i] = 1;
        }
    }

    void initializeZero() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                matrix[i][j] = 0;
            }
        }
    }

    void read() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                cin >> matrix[i][j];
            }
        }
    }

    void print() {
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < m; ++j) {
                cout << matrix[i][j] << " ";
            }
            cout << "\n";
        }
    }

    friend Matrix<T> operator*(const Matrix<T> &left, const Matrix<T> &right) {
        Matrix<T> res(left.n, right.m);
        if(left.m != right.n)
            return res;
        for (int i = 0; i < left.n; ++i) {
            for (int j = 0; j < right.m; ++j) {
                for (int k = 0; k < left.m; ++k) {
                    res.matrix[i][j] += left.matrix[i][k] * right.matrix[k][j];
                }
            }
        }
        return res;
    }

    friend bool operator==(const Matrix<T> &left, const Matrix<T> &right) {
        if(left.n != right.n || right.m != left.m)
            return false;

        for (int i = 0; i < left.n; ++i) {
            for (int j = 0; j < left.m; ++j) {
                if(left.matrix[i][j] != right.matrix[i][j])
                    return false;
            }
        }
        return true;
    }

    double det(Matrix<T> &matrix, int n){
        if(n == 1)
            return matrix[0][0];
        if(n == 2)
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

        Matrix<T> tmp(n - 1);
        double result = 0;

        for(int i = 0; i < n; ++i) {
            int vertical = 0;
            for(int j = 1; j < n; ++j) {
                int horizontal = 0;
                for(int k = 0; k < n; ++k) {
                    if(i != k) {
                        tmp[vertical][horizontal] = matrix[j][k];
                        ++horizontal;
                    }
                }
                ++vertical;
            }
            result += pow(-1, (double)i + 2) * matrix[0][i] * det(tmp, n - 1);
        }
        return result;
    }

    void swapR(int f, int s) {
        for(int i = 0; i < n; ++i)
            swap(matrix[f][i], matrix[s][i]);
    }

    void swapC(int f, int s) {
        for(int i = 0; i < n; ++i)
            swap(matrix[i][f], matrix[i][s]);
    }

private:
    vector<vector<T>> matrix;
};

// part1  ------  LU algorithm
//  -5x1 - x2 - 3x3 - x4 = 18
//  -2x1 + 8x3 - 4x4 = -12
//  -7x1 - 2x2 + 2x3 - 2x4 = 6
//  2x1 - 4x2 - 4x3 + 4x4 = -12
/*
4
-5 -1 -3 -1
-2 0 8 -4
-7 -2 2 -2
2 -4 -4 4
18 -12 6 -12
 Matrix L:
1 0 0 0
4 1 0 0
-2 -1 1 0
-1.5 -0.5 -0.5 1
Matrix U:
2 -2 -2 -7
0 8 4 26
0 0 4 14
0 0 0 4.5
 */

void LU_setter(Matrix<double> &A, Matrix<double> &L, Matrix<double> &U, Matrix<double> &Permutation, int n) {
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
        }

        for(int j = i + 1; j < n; ++j) {
            L[j][i] = U[j][i] / U[i][i]; // - (-a21/a11)
            for(int k = 0; k < n; ++k)
                U[j][k] -= L[j][i] * U[i][k]; // working with the other elements in the row
        }

    }
}

vector<double> LU_solver(Matrix<double> &A, vector<double> &R, Matrix<double> &L, Matrix<double> &U, int n) {
    vector<double> b(n, 0);
    vector<double> y(n, 0);
    vector<double> x(n, 0);
    Matrix<double> Permutation(n);

    LU_setter(A, L, U, Permutation, n);

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

Matrix<double> invertible(Matrix<double> &A, vector<double> &R, Matrix<double> &L, Matrix<double> &U, int n) {
    Matrix<double> result(n);
    vector<double> solutions(n, 0);

    for(int i = 0; i < n; ++i) {
        vector<double> y(n, 0);
        y[i] = 1;
        vector<double> itr = LU_solver(A, y, L, U, n);
        for(int k = 0; k < n; ++k)
            result[k][i] = itr[k];
    }

    return result;
}

int main() {
    int n;
    cin >> n;
    cout.precision(3);
    Matrix<double> A(n), U(n), L(n), M(n);
    A.read();
    vector<double> R(n);
    for(int i = 0; i < n; ++i) {
        cin >> R[i];
    }
    vector<double> answer(n);
    answer = LU_solver(A, R, L, U, n);

    cout << "Matrix A: " << endl;
    A.print();
    cout << "Matrix L: " << endl;
    L.print(); // lower
    cout << "Matrix U: " << endl;
    U.print(); // upper

    cout << "Answer: ";
    for(int i = 0; i < n; ++i) {
        cout << answer[i] << " ";
    }

    cout << "\n";
    double detA = A.det(A, n);
    cout << "Determinant A: " << detA << endl;
    cout << "\nInvertible matrix:\n";
    Matrix<double> I(n);
    I = invertible(A, R, L, U, n);
    I.print();
}
