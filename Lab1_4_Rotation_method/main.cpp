#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

double EPS = 0.001;
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

    bool tcheck(Matrix<T> &a) {
        double sum = 0;
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i != j) {
                    sum += pow(a[i][j], 2);
                }
            }
        }
        sum = sqrt(sum);
        return sum > EPS;
    }

    Matrix<T> transpose(Matrix<T> &a) {
        Matrix<T> tr(n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                tr[j][i] = a[i][j];
        return tr;
    }

    void maxNotDiagonal(int &i, int &j) {
        double max = abs(matrix[0][1]); i = 0; j = 1;
        for(int m = 0; m < n; ++m) {
            for(int l = 0; l < n; ++l) {
                if(m != l) {
                    if(abs(matrix[m][l]) > max) {
                        i = m;
                        j = l;
                        max = abs(matrix[m][l]);
                    }
                }
            }
        }
    }

    vector<T> RotationJacobiMethod(Matrix<T> &a_0, Matrix<T> &Ufinal, int &k) {
        k = 0;
        Matrix<T> a = a_0;
        Ufinal.initializeOne();
        Matrix<T> U(n);

        while(tcheck(a)) {
            int i_max = 0;
            int j_max = 0;
            a.maxNotDiagonal(i_max, j_max);

            double phi = M_PI / 4;
            if(abs(a[i_max][i_max] - a[j_max][j_max]) > ZERO)
                phi = 0.5 * atan(2 * a[i_max][j_max] / (a[i_max][i_max] - a[j_max][j_max]));

            U.initializeOne();
            U[i_max][i_max] = cos(phi); U[i_max][j_max] = - sin(phi);
            U[j_max][i_max] = sin(phi); U[j_max][j_max] = cos(phi);

            a = transpose(U) * a * U;
            Ufinal = Ufinal * U;
            ++k;
        }

        vector<T> lambda(n, 0);
        for(int i = 0; i < n; ++i)
            lambda[i] = a[i][i];
        return lambda;
    }


private:
    vector<vector<T>> matrix;
};

/*
0.001
3
 8  -3   9
-3   8  -2
 9  -2  -8
 */

int main() {
    int n;
    cin >> EPS;
    cin >> n;
    cout.precision(5);
    Matrix<double> A(n);
    A.read();
    Matrix<double> U(n);
    vector<double> lambda(n, 0);
    int k = 0;

    cout << "EPS =   " << EPS << endl;
    lambda = A.RotationJacobiMethod(A, U, k);
    cout << "Number of iterations:   " << k;
    cout << "\nEigenvalues:   ";
    for(int i = 0; i < n; ++i) {
        cout << lambda[i] << " ";
    }
    cout << "\nEigenvectors: " << endl;
    for(int j = 0; j < U.n; ++j) {
        cout << "#" << j << ":   ";
        for(int i = 0; i < U.n; ++i) {
            cout << U[i][j] << "  ";
        }
        cout << "\n";
    }
}

/*
EPS =   0.001
Number of iterations:   5
Eigenvalues:   14.1 5.95 -12.1
Eigenvectors:
#0:   0.783  -0.504  0.364
#1:   0.471  0.863  0.18
#2:   -0.405  0.0305  0.914
 */