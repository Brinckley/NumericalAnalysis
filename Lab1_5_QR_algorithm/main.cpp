#include <iostream>
#include <vector>
#include <math.h>
#include <complex>

using namespace std;

double EPS = 0.00001;
double ZERO = pow(10, -10);

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
        for(size_t i = 0; i < n; ++i) {
            for(size_t j = 0; j < m; ++j) {
                matrix[i][j] = 0;
            }
            matrix[i][i] = 1;
        }
    }

    void initializeZero() {
        for(size_t i = 0; i < n; ++i) {
            for(size_t j = 0; j < m; ++j) {
                matrix[i][j] = 0;
            }
        }
    }

    void read() {
        for(size_t i = 0; i < n; ++i) {
            for(size_t j = 0; j < m; ++j) {
                cin >> matrix[i][j];
            }
        }
    }

    void print() {
        for(size_t i = 0; i < n; ++i) {
            for(size_t j = 0; j < m; ++j) {
                if(abs(matrix[i][j]) < ZERO)
                    matrix[i][j] = 0;
                cout << matrix[i][j] << " ";
            }
            cout << "\n";
        }
    }

    friend Matrix<T> operator*(const Matrix<T> &left, const Matrix<T> &right) {
        Matrix<T> res(left.n, right.m);
        if(left.m != right.n)
            return res;
        for (size_t i = 0; i < left.n; ++i) {
            for (size_t j = 0; j < right.m; ++j) {
                for (size_t k = 0; k < left.m; ++k) {
                    res.matrix[i][j] += left.matrix[i][k] * right.matrix[k][j];
                }
            }
        }
        return res;
    }

    friend Matrix<T> operator-(const Matrix<T> &left, const Matrix<T> &right) {
        Matrix<T> res(left.n, left.n);
        for (size_t i = 0; i < left.n; ++i) {
            for (size_t j = 0; j < left.n; ++j) {
                res[i][j] = left[i][j] - right[i][j];
            }
        }
        return res;
    }

    double sign(double v) {
        return (v >= 0) ? 1 : -1;
    }

    bool fullMatrixCheck(Matrix<T> &am, vector<complex<T>> &old_lambdas, int& iter) {
        // under diagonal check
        double sum = 0;
        for(size_t i = 0; i < n - 1; ++i) { // for real numbers
            for (size_t j = i + 1; j < n; ++j) {
                sum += pow(am[i][j], 2);
            }
            if(sqrt(sum) <= EPS) { // exit condition for real numbers
                //cout << "first false!" << endl;
                return false;
            }
            sum = 0;
        }

        // complex check
        vector<complex<T>> new_lambdas(n, 0);
        for(size_t j = 0; j < n; ++j) {
            if(j < n - 1 && abs(am[j + 1][j]) > EPS) {
                T a = 1;
                T b = - am[j][j] - am[j + 1][j + 1];
                T c = - am[j][j + 1] * am[j + 1][j] + am[j][j] * am[j + 1][j + 1];
                T D = pow(b, 2) - 4 * a * c;
                if (D < ZERO) {
                    complex <T> l1 (- b / (2 * a), sqrt(- D) / (2 * a));
                    complex <T> l2 (- b / (2 * a), - sqrt(- D) / (2 * a));
                    new_lambdas[j] = l1;
                    new_lambdas[j + 1] = l2;
                    ++j;
                } else {
                    new_lambdas[j] = am[j][j];
                }
            } else {
                new_lambdas[j] = am[j][j];
            }
        }

        if(old_lambdas.size() != 0) {
            for(size_t i = 0; i < old_lambdas.size(); ++i) {
                complex<T> delta = new_lambdas[i] - old_lambdas[i];
                if(norm(delta) > EPS) {
                    old_lambdas = new_lambdas;
                    return true;
                }
            }
        } else {
            old_lambdas = new_lambdas; // nothing yet to compare to
            return true;
        }
        return false;
    }

    Matrix<T> Householder(vector<T> &v) {
        Matrix<T> vvt(n);
        Matrix<T> E(n);
        E.initializeOne();
        double vtv = 0;

        for (size_t i = 0; i < n; ++i) { // v * vT
            for (size_t j = 0; j < n; ++j) {
                vvt[i][j] = v[i] * v[j];
            }
        }
        for(size_t i = 0; i < n; ++i) { // vT * v
            vtv += pow(v[i], 2);
        }

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                vvt[i][j] = 2 * vvt[i][j] / vtv;
            }
        }
        Matrix<T> res(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                res[i][j] = E[i][j] - vvt[i][j];
            }
        }
        return res;
    }

    vector<complex<T>> QRalgorithm(Matrix<T> &a, Matrix<T> &q, Matrix<T> &r, int &iterations) {
        vector<complex<T>> lambdas;
        while(fullMatrixCheck(a, lambdas, iterations)) {
            iterations++;
            vector<T> v(n, 0);
            r = a;
            q.initializeOne();
            for (size_t i = 0; i < n - 1; ++i) {
                // I
                vector <T> v(n, 0);

                // II
                double sum = 0;
                for (size_t k = i; k < n; ++k) {
                    sum += pow(r[k][i], 2);
                }
                v[i] = r[i][i] + sign(r[i][i]) * sqrt(sum);

                // III
                for (size_t j = i + 1; j < n; ++j) {
                    v[j] = r[j][i];
                }

                Matrix<T> h = Householder(v);
                r = h * r;
                q = q * h;
            }
            a = r * q;
        }
        return lambdas;
    }


private:
    vector<vector<T>> matrix;
};

/*
0.00001
3
 1  2  5
-8  0 -6
 7 -9 -7
λ_1≈-11,932 λ_2≈-1,279 λ_3≈7,211
(-11.932,0)   (7.21055,0)   (-1.27853,0)

0.0001
 4
 0 5 5 2
 1 8 9 4
 1 2 8 2
 0 6 7 4
λ_1=1 λ_2≈-0,675 λ_3≈3,474 λ_4≈16,201
(16.2011,0)   (3.47407,0)   (0.999924,0)   (-0.675075,0)

 0.0001
 4
 2 7 -8 6
 4 4 0 -7
 -1 -3 6 3
 9 -7 -2 -8
 λ_1≈-16,026    λ_2≈5,366   λ_3≈7,330+1,880*i    λ_4≈7,330-1,880*i
(-16.0257,0)   (7.32309,1.85919)   (7.32309,-1.85919)   (5.3795,0)

 https://matrixcalc.org/vectors.html

 */

int main() {
    int n;
    cin >> EPS;
    cin >> n;
    cout.precision(6);
    Matrix<double> A(n);
    A.read();
    Matrix<double> U(n);
    vector<double> lambda(n, 0);
    int k = 0;
    int iterations = 0;
    cout << "EPS =   " << EPS << endl;
    Matrix<double> Q(n);
    Matrix<double> R(n);
    vector<complex<double>> lambdas = A.QRalgorithm(A, Q, R, iterations);
    cout << "Number of iterations:   " << iterations << endl;
    cout << "\nEigen values:   ";
    for(size_t i = 0; i < lambdas.size(); ++i)
        cout << lambdas[i] << "   ";
    cout << endl;


}