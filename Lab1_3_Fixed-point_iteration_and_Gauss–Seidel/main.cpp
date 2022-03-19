#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

double EPS = 0.001;

vector<double> operator-(const vector<double> &left, const vector<double> &right) {
    vector<double> res(left.size());
    for(int i = 0; i < left.size(); ++i) {
        res[i] = left[i] - right[i];
    }
    return res;
}

vector<double> operator+(const vector<double> &left, const vector<double> &right) {
    vector<double> res(left.size());
    for(int i = 0; i < left.size(); ++i) {
        res[i] = left[i] + right[i];
    }
    return res;
}

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

    friend vector<T> operator*(const Matrix<T> &left, const vector<T> &right) {
        vector<T> res(left.n, 0);
        for (int i = 0; i < left.n; ++i) {
            for (int j = 0; j < left.n; ++j) {
                res[i] += left.matrix[i][j] * right[j];
            }
        }
        return res;
    }

    double NormaMatrixC() {
        double max = 0;
        for(int j = 0; j < n; ++j) {
            max += abs(matrix[0][j]);
        }

        for(int i = 1; i < n; ++i){
            double sum = 0;
            for(int j = 0; j < n; ++j) {
                sum += abs(matrix[i][j]);
            }
            if(sum > max)
                max = sum;
        }
        return max;
    }

    double NormaC() {
        Matrix<T> C(n);
        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                C[i][j] = matrix[i][j];
            }
        }
        return C.NormaMatrixC();
    }

    double NormaVectorC(vector<T> vec) {
        double max = abs(vec[0]);
        for(int i = 1; i < n; ++i) {
            if(abs(vec[i]) > max) {
                max = abs(vec[i]);
            }
        }
        return max;
    }

    void JacobiSolver(Matrix<T> &a, vector<T> &b, Matrix<T> &alpha, vector<T> &beta) {
        for(int i = 0; i < n; ++i) {
            beta[i] = b[i] / a[i][i];
            for(int j = 0; j < n; ++j) {
                if(i != j) {
                    alpha[i][j] = - a[i][j] / a[i][i];
                    continue;
                }
                alpha[i][j] = 0;
            }
        }
    }

    vector<T> SolverFixedPointIteration(Matrix<T> &a, vector<T> &b, int &k) {
        Matrix<T> alpha(n);
        vector<T> beta(n);
        vector<vector<T>> x;
        k = 0;

        JacobiSolver(a, b, alpha, beta);
        double nalpha = alpha.NormaMatrixC();
        double eps_k = EPS + 1;

        x.push_back(beta);
        if(nalpha > 1) {
            cout << "False" << endl;
            return vector<T>();
        }
        while(eps_k > EPS) {
            x.push_back(beta + alpha * x[k]);
            ++k;
            eps_k = (nalpha / (1 - nalpha)) * NormaVectorC(x[k] - x[k - 1]);
        }
        return x[k];
    }

    vector<T> SolverGaussSeidel(Matrix<T> &a, vector<T> &b, int &k) {
        Matrix<T> alpha(n);
        vector<T> beta(n);
        vector<vector<T>> x;
        k = 0;

        JacobiSolver(a, b, alpha, beta);
        double normaC = alpha.NormaC();
        double nalpha = alpha.NormaMatrixC();
        double eps_k = EPS + 1;

        x.push_back(beta);
        while(eps_k > EPS) {
            vector<T> xn(n, 0);
            for(int i = 0; i < n; ++i) {
                xn[i] = beta[i];
                for(int j = 0; j < i; ++j) {
                    xn[i] += alpha[i][j] * xn[j];
                }
                for(int j = i; j < n; ++j) {
                    xn[i] += alpha[i][j] * x[k][j];
                }
            }
            x.push_back(xn);
            ++k;
            eps_k = (normaC / (1 - nalpha)) * NormaVectorC(x[k] - x[k - 1]);
        }
        return x[k];
    }


private:
    vector<vector<T>> matrix;
};

/*
0.001
4
21 -6 -9 -4
-6 20 -4 2
-2 -7 -20 3
4 9 6 24

127 -144 236 -5
 */

int main() {
    int n;
    cin >> EPS;
    cin >> n;
    cout.precision(3);
    Matrix<double> A(n);
    A.read();
    vector<double> B(n, 0);
    for(int i = 0; i < n; ++i)
        cin >> B[i];

    cout << "EPS = " << EPS << endl;
    int counterIterationFixedPoint = 0;
    cout << "\nFixed Point Iteration method" << endl;
    cout << "Solutions:  ";
    vector<double> x = A.SolverFixedPointIteration(A, B, counterIterationFixedPoint);
    for(int i = 0; i < n; ++i)
        cout << x[i] << " ";
    cout << "\nFixed Point Iteration iterations number: " << counterIterationFixedPoint << endl;

    cout << "\nGauss Seidel method" << endl;
    int counterSeidel = 0;
    cout << "Solutions:  ";
    vector<double> y = A.SolverGaussSeidel(A, B, counterSeidel);
    for(int i = 0; i < n; ++i)
        cout << x[i] << " ";
    cout << "\nGauss Seidel method iterations number: " << counterSeidel << endl;
}

/*
EPS = 0.001

Fixed Point Iteration method
Solutions:  1 -9 -8 5
Fixed Point Iteration iterations number: 16

Gauss Seidel method
Solutions:  1 -9 -8 5
Gauss Seidel method iterations number: 7
 */