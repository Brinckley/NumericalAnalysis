#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

double EPS = 0.001;

static double f(double x) {
    return x * exp(x) + x * x - 1;
}

static double df(double x) {
    return (x + 1) * exp(x) + 2 * x;
}

static double Newton(double xo, int &iter) {
    double x_0 = 0;
    double x_k1 = xo;
    do {
        x_0 = x_k1;
        iter++;
        if(df(x_0) == 0) break; // f(x) * f''(x) > 0 -> f'(x) != 0
        x_k1 = x_0 - f(x_0) / df(x_0);
    } while (abs(x_k1 - x_0) > EPS);
    return x_k1;
}

static double phi(double x) {
    return - sqrt(- x * exp(x) + 1);
}

static double dphi(double x) {
    return (- 1) * exp(x) * (x + 1) / (2 * sqrt(1 - x * exp(x)));
}

static double Iteration(double l, double r, int &iter) {
    // graphic check -> phi ∈ [l, r]
    double q = min(max(abs(dphi(l)), abs(dphi(r))), 1 - EPS); // monotonous!!!! for my max is right value
    double q_ = q / (1 - q);

    double x_0 = 0;
    double x_k1 = l;
    do {
        x_0 = x_k1;
        iter++;
        x_k1 = phi(x_0);;
    } while (EPS < q_ * abs(x_k1 - x_0));
    return x_k1;
}

int main() {
    cin >> EPS;
    double xo = 0;
    cin >> xo;
    int iter = 0;
    double x1 = Newton(xo, iter);
    cout << "Iter number:  " << iter << endl;
    cout << "Solution:  " << x1 << endl;

    iter = 0;
    double l, r;
    cin >> l >> r;
    double x2 = Iteration(l, r, iter);
    cout << "Iter number:  " << iter << endl;
    cout << "Solution:  " << x2 << endl;

}

/*
0.1
-1.5
-3 0

Iter number:  2
Solution:  -1.16825
Iter number:  2
Solution:  -1.16917
 */
