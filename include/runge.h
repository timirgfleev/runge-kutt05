#include <cmath>
#include <cfloat>
#include <vector>

class RungeKuttaSolver
{
public:
    explicit RungeKuttaSolver(double (*f)(double, double)) : f(f){};

    /*
     * Runge-Kutta 4th order
     * @param x_start - start x
     * @param x_end - end x
     * @param y_0 - y(x_start)
     * @param h - step
     * @param x_l - last x
     * @param back - if true, then x_start is end and x_end is start
     * @return y(x_end)
     */
    double runge_kutt_o4(double x_start, double x_end, double y_0, double h, double &x_l, bool back = false)
    {
        int step_count = int((x_end - x_start) / h);

        double x = x_start, y, yp = y_0;
        double res;

        int start = back ? step_count - 1 : 1; // start skip first
        int end = back ? 0 : step_count;
        int step = back ? -1 : 1;
        back ? h = -h : h = h;

        std::vector<double> y_arr(step_count + 1, 0);
        back ? y_arr[step_count] = y_0 : y_arr[0] = y_0;

        for (int i = start; back ? i >= end : i <= end; i += step)
        {
            x = x_start + h * i;
            y = calc_y1(x, yp, h);
            yp = y;

            y_arr[i] = y;
        }

        x_l = x;
        res = y;

        return res;
    }

    /*
     * @return machine epsilon
     */
    double macheps()
    {
        double e = 1.0;

        while (1.0 + (e / 2.0) > 1.0)
            e /= 2.0;
        return e;
    }

    /*
     * @return min step value according to machine epsilon
     */
    double calc_h_min(double a, double b)
    {
        return macheps() * std::max(std::max(std::abs(a),
                                             std::abs(b)),
                                    DBL_MIN); // DBL_MIN minimal double value
    }

    /*
     * @param y0 - y0
     * @param y1 - y1
     * @param p - order of method (4 by default)
     * @return error value (double)
     */
    double runge_calc_error(double y0, double y1, int p = 4)
    {
        double u = std::abs(y0 - y1);
        double d = (pow(2, p) - 1);
        return u / d;
    }

protected:
    //reference to function f(x, y)
    double (*f)(double, double);

    double calc_k1(double x0, double y0, double h)
    {
        return h * (f(x0, y0));
    }

    double calc_k2(double x0, double y0, double h,
                   double k1)
    {
        return h * (f(x0 + h / 3, y0 + k1 / 3));
    }

    double calc_k3(double x0, double y0, double h,
                   double k1,
                   double k2)
    {
        return h * (f(x0 + 2 * h / 3, y0 - k1 / 3 + k2));
    }

    double calc_k4(double x0, double y0, double h,
                   double k1,
                   double k2,
                   double k3)
    {
        return h * (f(x0 + h, y0 + k1 - k2 + k3));
    }

    /*
     * Runge-Kutta 4th order formula
     * @param x - x0
     * @param y0 - y0
     * @param h - step
     * @return y1
     */
    double calc_y1(double x, double y0, double h)
    {
        double k1 = calc_k1(x, y0, h);
        double k2 = calc_k2(x, y0, h, k1);
        double k3 = calc_k3(x, y0, h, k1, k2);
        double k4 = calc_k4(x, y0, h, k1, k2, k3);
        double y1 = y0 + (k1 + 3 * k2 + 3 * k3 + k4) / 8;
        return y1;
    }
};
