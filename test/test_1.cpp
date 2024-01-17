#include <cmath>
#include <cassert>
#include "runge.h"

#include <iostream>

const double h = 0.00001; //step for tests (tests may fail if h is not compatible enough with eps)
//I mean, some functions may have bad precision for this h so they cant fall into eps
const double eps = 1e-2; //precision for test results

double test_function(double x, double y) {
    return x + y;
}

void test_runge_kutt_o4() {
    RungeKuttaSolver solver(test_function);
    double x_start = 0.0;
    double x_end = 1.0;
    double y_0 = 0.0;
    double x_l;
    double result = solver.runge_kutt_o4(x_start, x_end, y_0, h, x_l);

    // Check the result with some expected value
    double expected = M_E - 2.0;
    assert(std::abs(result - expected) < eps);
}

double test_function_2(double x, double y) {
    return x * x + y;
}

void test_runge_kutt_o4_2() {
    RungeKuttaSolver solver(test_function_2);
    double x_start = 0.0;
    double x_end = 1.0;
    double y_0 = 1.0;
    double x_l;
   
    double result = solver.runge_kutt_o4(x_start, x_end, y_0, h, x_l);

    // Check the result with some expected value
    double expected = 3 * M_E - 5; // This is an approximation
    assert(std::abs(result - expected) < eps);
}

double test_function_3(double x, double y) {
    return sin(x) + y;
}

void test_runge_kutt_o4_3() {
    RungeKuttaSolver solver(test_function_3);
    double x_start = 0.0;
    double x_end = 1.0;
    double y_0 = 1.0;
    double x_l;
     
    double result = solver.runge_kutt_o4(x_start, x_end, y_0, h, x_l);

    // Check the result with some expected value
    double expected = 0.5 * (3 * M_E - sin(1) - cos(1)); 
    assert(std::abs(result - expected) < eps);
}

double test_function_4(double x, double y)
{
    return x - y;
}

void test_runge_kutt_o4_4()
{
    RungeKuttaSolver solver(test_function_4);
    double x_start = 0.0;
    double x_end = 1.0;
    double y_0 = 2.0;
    double x_l;



    double result = solver.runge_kutt_o4(x_start, x_end, y_0, h, x_l, true);

    
    // Check the result with some expected value
    double expected = 2 * M_E - 1.0;
    assert(std::abs(result - expected) < eps);
}

int main() {
    test_runge_kutt_o4();
    test_runge_kutt_o4_2();
    test_runge_kutt_o4_3();
    test_runge_kutt_o4_4();
    return 0;
}