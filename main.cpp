#include "include/math_function.h"
#include "include/runge.h"
#include "include/file_io.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// рунге-кутта 4 порядка точности (формула 33)

/*
 * @return icod
 * @param a - start x
 * @param b - end x
 * @param c - x0
 * @param yc - y(x0)
 * @param h - step
 * @param eps - error
 * @param out - output struct
 * @return icod
 */
int calc(double a, double b, double c, double yc, double h, double eps, OutStruct &out)
{
    auto r = RungeKuttaSolver(f);

    // init vars
    bool is_back = (c == b);
    int icod;
    double y_last, y_last_next, x_last, x_last_next;
    double error = -1, error_p;

    double h_min = r.calc_h_min(a, b);

    bool is_continue = true;

    // init compare buffer value
    y_last_next = r.runge_kutt_o4(a, b, yc, h, x_last_next, is_back);
    int iter = 0;

    std::cout << "i | y | errorP | error | h" << std::endl;

    /*
     * It does following:
     * Calc for h and h/2
     * untill error < eps
     * or error_p <= error (that means error is increasing)
     */
    while (is_continue)
    {
        iter++;

        y_last = y_last_next;
        //x_last = x_last_next;

        h /= 2;
        y_last_next = r.runge_kutt_o4(a, b, yc, h / 2, x_last_next, is_back);

        error_p = error;
        error = r.runge_calc_error(y_last, y_last_next); // правило рунге

        std::cout << iter << ' ' << y_last_next << ' ' << error_p << ' ' << error << ' ' << h << std::endl;

        if (error < eps)
        {
            is_continue = false;
            icod = 0;
        }
        if (error_p <= error && error_p != -1)
        {
            is_continue = false;
            icod = 1;
        }
        if (h < h_min)
        {
            is_continue = false;
            icod = 2;
        }
    }

    out = OutStruct{error, h, icod, x_last_next, y_last_next};
    return icod;
}

int init_calc(const std::vector<double> &inp, OutStruct &out)
{
    return calc(inp[0], inp[1], inp[2], inp[3], inp[4], inp[5], out);
}

int main()
{
    int icod = 0;

    std::string INP_FILE_PATH = "data.txt";
    std::string O_FILE_PATH = "rez.txt";

    std::ifstream ifs(INP_FILE_PATH);
    if (!ifs.is_open())
    {
        icod = -1; // file cannot be accessed
        std::cout << "No input file" << std::endl;
        return icod;
    }

    std::vector<double> input_v;

    if (!read_file(ifs, input_v))
    {
        icod = -2; // cannot read file
        std::cout << "Cant read inp file" << std::endl;
        return icod;
    }

    OutStruct out;
    icod = init_calc(input_v, out);

    std::cout << "icod: " << icod << std::endl;
    std::cout << "Calc stop with code: ";

    switch (icod)
    {
    case 0:
        std::cout << "error < eps, OK" << std::endl;
        break;
    case 1:
        std::cout << "Error is increasing" << std::endl;
        break;
    case 2:
        std::cout << "h < h_min" << std::endl;
        break;

    default:
        break;
    }

    std::ofstream ofs(O_FILE_PATH);

    if (!ofs.is_open())
    {
        icod = -1; // file cannot be accessed
        std::cout << "no access file out" << std::endl;
        return icod;
    }
    ofs << out;
    return icod;
}
	