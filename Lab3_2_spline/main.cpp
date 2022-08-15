#include "polynom.hpp"
#include "spline.hpp"

int main() {
    double _x = 0.8;
    vector<double> x {0.1, 0.5, 0.9, 1.3, 1.7};
    vector<double> f {-2.2026, -0.19315, 0.79464, 1.5624, 2.2036};
    spline s(x, f);
    s.builder();
    cout << s.solver(_x) << endl;

}
