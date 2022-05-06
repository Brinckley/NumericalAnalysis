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
/*
    pair<T, T> newtonSolverSimple(T x1_0, T x2_0, int &iter) {
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
*/
    vector<T> newtonSolver(vector<T> x_0, int &iter) {
        vector<T> x_k = x_0;

        do {
            x_0 = x_k;

            matrix<T> J = myJacobiTaskBuilder(x_0);
            matrix<T> Ji = J.invertible(J, x_0.size());
            vector<T> f = myFTaskBuilder(x_0);
            vector<T> dx = Ji * f;

            iter++;
            x_k = x_0 - dx;
        } while (vecNorm(x_0, x_k) > EPS);

        cout << "Func1 result:  " << func1(x_k[0], x_k[1])  << endl;
        cout << "Func2 result:  " << func2(x_k[0], x_k[1])  << endl;
        return x_k;
    }

private:
    T func1(T x1, T x2) {
        return 2 * x1 - cos(x2);
    }

    T func2(T x1, T x2) {
        return 2 * x2 - exp(x1);
    }

    matrix<T> myJacobiTaskBuilder(vector<T> x) {
        return matrix<T>(2, vector<T> {2, - exp(x[0]), sin(x[1]), 2});
    }

    vector<T> myFTaskBuilder(vector<T> x) {
        return vector<T>{func1(x[0], x[1]), func2(x[0], x[1])};
    }

    T vecNorm(vector<T> v1, vector<T> v2) {
        T sum = 0;
        for(int i = 0; i < v1.size(); ++i) {
            sum += abs(v1[i] - v2[i]);
        }
        return sum;
    }
/*
    void updateNewtonMatrices(T x1_0, T x2_0, matrix<T> &A1, matrix<T> &A2, matrix<T> &J) {
        A1 = matrix<T>(2, vector<T> {func1(x1_0, x2_0), d2func1(x1_0, x2_0),
                                     func2(x1_0, x2_0), d2func2(x1_0, x2_0)});
        A2 = matrix<T>(2, vector<T> {d1func1(x1_0, x2_0), func1(x1_0, x2_0),
                                     d1func2(x1_0, x2_0), func2(x1_0, x2_0)});
        J = matrix<T>(2, vector<T> {d1func1(x1_0, x2_0), d2func1(x1_0, x2_0),
                                    d1func2(x1_0, x2_0), d2func2(x1_0, x2_0)});
    }*/

    T vec2Norm(pair<T, T> x1, pair<T, T> x2) {
        T sum = 0;
        sum += abs(x1.first - x2.first);
        sum += abs(x1.second - x2.second);
        return sum;
    }


    T phi1(T x1, T x2) { return cos(x2) / 2; }
    T phi2(T x1, T x2) { return exp(x1) / 2; }

    T d1phi1(T x1, T x2) { return 0; }
    T d2phi1(T x1, T x2) { return - sin(x2) / 2; }
    T d1phi2(T x1, T x2) { return exp(x1) / 2; }
    T d2phi2(T x1, T x2) { return 0; }

    pair<T, T> iterationSolver(T x1_0, T x2_0, int &iter) {
        T x1_k = x1_0;
        T x2_k = x2_0;
        matrix<T> dphi(2, {d1phi1(x1_0, x2_0), d2phi1(x1_0, x2_0),
                           d1phi2(x1_0, x2_0), d2phi2(x1_0, x2_0)});

        T q = 0;
        T q_ = q / (1 - q);

        do{
            x1_k = x1_0;
            x2_k = x2_0;


        } while (EPS < q_ * vec2Norm({x1_k, x2_k}, {x1_0, x2_0}));

        return {x1_k, x2_k};
    }


};


#endif //LAB2_2_ITERATIONS_AND_NEWTON_SYSTEM_2L2_WORKER_HPP
