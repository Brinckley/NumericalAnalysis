//
// Created by alex on 14.08.22.
//

#ifndef LAB3_2_SPLINE_POLYNOM_HPP
#define LAB3_2_SPLINE_POLYNOM_HPP

#include <iostream>
#include <vector>
#include <cmath>

using namespace  std;

class polynom {
private:
    constexpr static double EPS = 0.00001;

    vector<double> poly;
    int n;

    void clear() {
        while (n > 1 and abs(poly.back()) < EPS) {
            poly.pop_back();
            n--;
        }
    }

public:
    polynom() {
        n = 1;
        poly = vector<double>(n);
    }

    polynom(int n) {
        this->n=  n;
        poly = vector<double>(n);
    }

    polynom(const vector<double> &v) {
        poly = v;
        n = v.size();
    }

    double& operator[](int i) {
        return poly[i];
    }

    friend polynom operator+(const polynom &lhs, const polynom &rhs) {
        polynom res(max(lhs.n, rhs.n));
        for (int i = 0; i < res.n; ++i) {
            if(i < lhs.n) {
                if(i < rhs.n) {
                    res[i] = lhs.poly[i] + rhs.poly[i];
                } else {
                    res[i] = lhs.poly[i];
                }
            } else {
                res[i] = rhs.poly[i];
            }
        }
        return res;
    }

    friend polynom operator-(const polynom &lhs, const polynom &rhs) {
        polynom res(max(lhs.poly.size(), rhs.poly.size()));
        for (int i = 0; i < res.n; ++i) {
            if(i < lhs.n) {
                if(i < rhs.n) {
                    res[i] = lhs.poly[i] - rhs.poly[i];
                } else {
                    res[i] = lhs.poly[i];
                }
            } else {
                res[i] = -rhs.poly[i];
            }
        }
        return res;
    }

    friend polynom operator*(double lambda, const polynom &p) {
        polynom res(p.poly);
        for (int i = 0; i < res.poly.size(); ++i) {
            res[i] *= lambda;
        }
        return res;
    }

    friend polynom operator/(const polynom &p, double a) {
        polynom res(p);
        for (int i = 0; i < res.poly.size(); ++i) {
            res.poly[i] /= a;
        }
        return res;
    }

    friend polynom operator*(const polynom &lhs, const polynom &rhs) {
        polynom res(lhs.poly.size() + rhs.poly.size());
        for (int i = 0; i < lhs.poly.size(); ++i) {
            for (int j = 0; j < rhs.poly.size(); ++j) {
                res.poly[i + j] += lhs.poly[i] * rhs.poly[j];
            }
        }
        res.clear();
        return res;
    }

    friend ostream& operator<<(ostream &out, const polynom &p) {
        for (int i = 0; i < p.n; ++i) {
            out << p.poly[i] << " * " <<  "x^" << i;
            if (i != p.n - 1) {
                out << " + ";
            } else {
                out << endl;
            }
        }
        return out;
    }

    double calculate(double xi) {
        double res = 0;
        for(int i = 0; i < n; ++i) {
            res += poly[i] * pow(xi, i);
        }
        return res;
    }
};


#endif //LAB3_2_SPLINE_POLYNOM_HPP
