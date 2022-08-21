//
// Created by alex on 20.08.22.
//

#ifndef LAB3_4_DIFF_DIFF_HPP
#define LAB3_4_DIFF_DIFF_HPP

#include <vector>

using namespace std;

class Diff {
public:
    Diff(vector<double> &x, vector<double> &y) {
        this->x = x;
        this->y = y;
        n = x.size();
    }

    double diff_first(double x_) {
        int x_index = 0;
        for (int i = 0; i < n - 2; ++i) {
            if (x_ >= x[i] && x_ <= x[i + 1]) {
                x_index = i;
                break;
            }
        }

        return
        (y[x_index + 1] - y[x_index + 0]) / (x[x_index + 1] - x[x_index + 0])
        +
        ((y[x_index + 2] - y[x_index + 1]) / (x[x_index + 2] - x[x_index + 1])
        - (y[x_index + 1] - y[x_index + 0]) / (x[x_index + 1] - x[x_index + 0]))
        * (2.0 * x_ - x[x_index] - x[x_index + 1])
        / (x[x_index + 2] - x[x_index]);
    }

    double diff_second(double x_) {
        int index = 0;
        for (int i = 0; i < n - 2; ++i) {
            if (x_ >= x[i] && x_ <= x[i + 1]) {
                index = i;
                break;
            }
        }
        return
        2.0
        * ((y[index + 2] - y[index + 1]) / (x[index + 2] - x[index + 1])
        - (y[index + 1] - y[index]) / (x[index + 1] - x[index]))
        / (x[index + 2] - x[index]);
    }
private:
    int n;
    vector<double> x;
    vector<double> y;
};


#endif //LAB3_4_DIFF_DIFF_HPP
