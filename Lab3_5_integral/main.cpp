#include <iostream>

#include "integral.hpp"

int main() {
    double h1 = 0.5;
    double h2 = 0.25;
    double x0 = 0;
    double xk = 2;

    integral_worker worker(x0, xk, h1);

    cout << "\nStep = " << h1 << endl;
    double r1 = worker.rectangle_method();
    double t1 = worker.trapezoid_method();
    double s1 = worker.simpson_method();
    cout << "Rectangle method : " << r1 << endl;
    cout << "Trapezoid method : " << t1 << endl;
    cout << "Simpson method   : " << s1 << endl;


    worker = integral_worker(x0, xk, h2);

    cout << "\nStep = " << h2 << endl;
    double r2 = worker.rectangle_method();
    double t2 = worker.trapezoid_method();
    double s2 = worker.simpson_method();
    cout << "Rectangle method : " << r2 << endl;
    cout << "Trapezoid method : " << t2 << endl;
    cout << "Simpson method   : " << s2 << endl;


    double k = 0.001;
    cout << "\nRunge-Romberg-Richardson" << endl;
    cout << "Rectangle method : " << integral_worker::Runge_Romberg_Richardson_method(r1, r2, k, 2) << endl;
    cout << "Rectangle method : " << integral_worker::Runge_Romberg_Richardson_method(t1, t2, k, 2) << endl;
    cout << "Simpson method   : " << integral_worker::Runge_Romberg_Richardson_method(s1, s2, k, 4) << endl;  // parabola

}
