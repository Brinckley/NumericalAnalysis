#include <iostream>
#include <vector>

#include "squares.hpp"

using namespace std;

int main() {
    vector<double> x { 0.1, 0.5, 0.9, 1.3, 1.7, 2.1 };
    vector<double> y { -2.2026, -0.19315, 0.79464, 1.5624, 2.2306, 2.8419 };
    /*
     * ----- 1 -----
     * -1.77444 * x^0 + 2.37582 * x^1
     * PHI1 = 0.985412
     * ----- 2 -----
     * -2.46046 * x^0 + 4.40618 * x^1 + -0.922892 * x^2
     * PHI2 = 0.171386
     */

    // vector<double> x { -1.7, -1.2, -0.7, -0.2, 0.3, 0.8 };
    // vector<double> y { 0.52796, 0.43372, 0.24333, 0.03275, 0.12149, 1.4243};
    /*
     * ----- 1 -----
     * 0.549667 * x^0 + 0.190539 * x^1
     * PHI1 = 1.12033
     * ----- 2 -----
     * 0.244887 * x^0 + 0.711367 * x^1 + 0.578698 * x^2
     * PHI2 = 0.338921
     */

    cout << "----- 1 -----" << endl;
    squares square(x, y, 2);
    square.square_build();
    square.printPoly();
    cout << "PHI1 = " << square.squares_result() << endl;

    cout << "----- 2 -----" << endl;
    squares square1(x, y, 3);
    square1.square_build();
    square1.printPoly();
    cout << "PHI2 = " << square1.squares_result() << endl;
}
