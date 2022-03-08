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
/*
5

18 -9
2  -9  -4
-9  21 -8
-4 -10  5
7   12

-81 71 -39 64 3
*/
    void read() {
        cin >> b[0] >> c[0];
        for (int i = 1; i < n - 1; ++i) {
            cin >> a[i] >> b[i] >> c[i];
        }
        cin >> a[n - 1] >> b[n - 1];
    }

    vector<T> TridiagonalSolver(vector<T> &d) {
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

int main() {
    int n;
    cin >> n;
    cout.precision(3);
    TridiagonalWorker<double> tri(n);
    tri.read();
    vector<double> d(n, 0);
    for(int i = 0; i < n; ++i)
        cin >> d[i];
    vector<double> result = tri.TridiagonalSolver(d);
    for(int i = 0; i < n; ++i)
        cout << result[i] << " ";
    cout << endl;
}