#include <cmath>
#include <cassert>
#include "runge.h"
#include <limits>

double test_function(double x, double y) {
    return x + y;
}

void test_macheps() {
    RungeKuttaSolver solver(test_function);
    double macheps = solver.macheps();
    assert(std::abs(macheps - std::numeric_limits<double>::epsilon()) < 1e-10);
}

void test_calc_h_min() {
    RungeKuttaSolver solver(test_function);
    double a = 1.0, b = 2.0;
    double h_min = solver.calc_h_min(a, b);
    assert(h_min > 0);
}

void test_runge_calc_error() {
    RungeKuttaSolver solver(test_function);
    double y0 = 1.0, y1 = 1.1;
    double error = solver.runge_calc_error(y0, y1);
    assert(error > 0);
}

void test_get_x_y() {
    RungeKuttaSolver solver(test_function);
    double x_start = 0.0;
    double x_end = 1.0;
    double y_0 = 0.0;
    double h = 0.0001;
    double x_l;
    double result = solver.runge_kutt_o4(x_start, x_end, y_0, h, x_l);

    std::vector<double> x_arr = solver.get_x();
    std::vector<double> y_arr = solver.get_y();

    assert(x_arr.size() == y_arr.size());
    assert(x_arr[0] == x_start);
    assert(x_arr.back() == x_end);
}

int main() {
    test_macheps();
    test_calc_h_min();
    test_runge_calc_error();
    test_get_x_y();
    return 0;
}