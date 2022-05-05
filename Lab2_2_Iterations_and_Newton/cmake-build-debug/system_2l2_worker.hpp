#ifndef LAB2_2_ITERATIONS_AND_NEWTON_SYSTEM_2L2_WORKER_HPP
#define LAB2_2_ITERATIONS_AND_NEWTON_SYSTEM_2L2_WORKER_HPP

#include "matrix.hpp"
#include <math.h>
using namespace std;

template <class T>
class system_2l2_worker {
public:
    T EPS = 0.001;
    system_2l2_worker(T EPS) {
        this->EPS = EPS;
    }

    pair<T, T> newtonSolver(T x1_0, T x2_0, int &iter) {
        T x1_k = x1_0;
        T x2_k = x2_0;

        do {
            iter++;
            x1_0 = x1_k; x2_0 = x2_k;
            matrix<T> A1, A2, J;
            updateNewtonMatrices(x1_0, x2_0, A1, A2, J);

            if(J.det2() == 0) break;
            x1_k = x1_0 - A1.det2() / J.det2();
            x2_k = x2_0 - A2.det2() / J.det2();
        } while (vec2Norm({x1_0, x2_0}, {x1_k, x2_k}) > EPS);

        cout << "Func1 result:  " << func1(x1_k, x2_k)  << endl;
        cout << "Func2 result:  " << func2(x1_k, x2_k)  << endl;
        return {x1_k, x2_k};
    }

private:
    T func1(T x1, T x2) {
        return 2 * x1 - cos(x2);
    }

    T d1func1(T x1, T x2) {
        return 2;
    }

    T d2func1(T x1, T x2) {
        return sin(x2);
    }

    T func2(T x1, T x2) {
        return 2 * x2 - exp(x1);
    }

    T d1func2(T x1, T x2) {
        return - exp(x1);
    }

    T d2func2(T x1, T x2) {
        return 2;
    }


    void updateNewtonMatrices(T x1_0, T x2_0, matrix<T> &A1, matrix<T> &A2, matrix<T> &J) {
        A1 = matrix<T>(2, vector<T> {func1(x1_0, x2_0), d2func1(x1_0, x2_0),
                                       func2(x1_0, x2_0), d2func2(x1_0, x2_0)});
        A2 = matrix<T>(2, vector<T> {d1func1(x1_0, x2_0), func1(x1_0, x2_0),
                                       d1func2(x1_0, x2_0), func2(x1_0, x2_0)});
        J = matrix<T>(2, vector<T> {d1func1(x1_0, x2_0), d2func1(x1_0, x2_0),
                                      d1func2(x1_0, x2_0), d2func2(x1_0, x2_0)});
    }

    T vec2Norm(pair<T, T> x1, pair<T, T> x2) {
        T sum = 0;
        sum += abs(x1.first - x2.first);
        sum += abs(x1.second - x2.second);
        return sum;
    }



};


#endif //LAB2_2_ITERATIONS_AND_NEWTON_SYSTEM_2L2_WORKER_HPP
