#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cfloat>
//рунге-кутта 4 порядка точности (формула 33)

using namespace std;

//todo: make f pass into class by construct or something; also split code to diff file as lib
double f(double x, double y) {
    return pow((x-y), 2) + 1;
}

bool read_file(ifstream &ifs, vector<double> &res) {
    const int param_count = 6;
    res = vector<double>(param_count);
    for (int i = 0; i < param_count; ++i) {
        ifs >> res[i];
    }
    return ifs.fail();
}

struct OutStruct {
    //OutStruct() : er(0), h(0), icod(-3), x_l(0), y_l(0) {};
    double er;//абсолютная погрешность
    double h; //шаг
    int icod;
    double x_l; //последний x
    double y_l; // последний y
};

ofstream &operator<<(ofstream &ofs, const OutStruct &e) {
    const char sep = ' ';
    ofs << e.er << sep << e.h << sep << e.icod << endl;
    ofs << e.x_l << sep << e.y_l << endl;
    return ofs;
}

bool out_file(ofstream &ofs, OutStruct &data) {
    ofs << data;
    return ofs.fail();
}


class Runge{
public:
    int calc(double a, double b, double c, double yc, double h, double eps, OutStruct &out) {
        bool is_back = (c == b);
        int icod;
        double y_last, y_last_next, x_last, x_last_next;
        double error = -1, error_p;

        double h_min = calc_h_min(a, b);

        bool is_continue = true;
        y_last_next = runge_kutt_o4(a, b, yc, h, x_last_next, is_back);
        int iter = 0;
        cout << "i | y | errorP | error" << endl;
        while (is_continue) {
            iter++;

            y_last = y_last_next;
            x_last = x_last_next;

            h /= 2;
            y_last_next = runge_kutt_o4(a, b, yc, h / 2, x_last_next, is_back);

            error_p = error;
            error = runge_calc_error(y_last, y_last_next);//правило рунге

            cout << iter << ' ' << y_last_next << ' ' << error_p << ' '  << error << endl;

            if (error < eps) {
                is_continue = false;
                icod = 0;
            }
            if (error_p <= error && error_p != -1) {
                is_continue = false;
                icod = 1;
            }
            if (h < h_min) {
                is_continue = false;
                icod = 2;
            }
        }

        out = OutStruct{error, h, icod, x_last_next, y_last_next};
        return icod;
    }

    double runge_kutt_o4(double x_start, double x_end, double y_0, double h, double &x_l, bool back = false) {
        int step_count = int((x_end - x_start) / h);

        double x = x_start, y, yp = y_0;
        //double res = yp * h;
        double res;
        if (!back) {
            for (int i = 1; i <= step_count; ++i) {
                x = x_start + h * i;
                y = calc_y1(x, yp, h);
                //res += y * h;
                yp = y;
            }
            x_l = x;
            res = y;
        } else {
            for (int i = step_count; i > 0; --i) {
                x = x_start + h * i;
                y = calc_y1(x, yp, h);
                //res += y * h;
                yp = y;
            }
            x_l = x;
            res = y;
        }

        return res;
    }
protected:
    double calc_k1(double x0, double y0, double h) {
        return h * (f(x0, y0));
    }

    double calc_k2(double x0, double y0, double h,
                   double k1) {
        return h * (f(x0 + h / 3, y0 + k1 / 3));
    }

    double calc_k3(double x0, double y0, double h,
                   double k1,
                   double k2) {
        return h * (f(x0 + 2 * h / 3, y0 - k1 / 3 + k2));
    }


    double calc_k4(double x0, double y0, double h,
                   double k1,
                   double k2,
                   double k3) {
        return h * (f(x0 + h, y0 + k1 - k2 + k3));
    }

    double calc_y1(double x, double y0, double h) {
        double k1 = calc_k1(x, y0, h);
        double k2 = calc_k2(x, y0, h, k1);
        double k3 = calc_k3(x, y0, h, k1, k2);
        double k4 = calc_k4(x, y0, h, k1, k2, k3);
        double y1 = y0 + (k1 + 3 * k2 + 3 * k3 + k4) / 8;
        return y1;
    }

    double runge_calc_error(double y0, double y1, int p = 4) {
        double u = abs(y0 - y1);
        double d = (pow(2, p) - 1);
        return u / d;
    }

    double macheps() {
        double e = 1.0;

        while (1.0 + (e / 2.0) > 1.0)
            e /= 2.0;
        return e;
    }

    double calc_h_min(double a, double b) {
        return macheps() * max(max(abs(a), abs(b)), DBL_MIN);//DBL_MIN Минимальное положительное значение.
    }
};

/*
 * res:
 * a, b, c, yc,
 * H start, eps
 */

int init_calc(const vector<double> &inp, OutStruct &out) {
    return Runge().calc(inp[0], inp[1], inp[2], inp[3], inp[4], inp[5], out);
}

int main() {
    int icod = 0;

    string INP_FILE_PATH = "data.txt";
    string O_FILE_PATH = "rez.txt";

    ifstream ifs(INP_FILE_PATH);
    if (!ifs.is_open()) {
        icod = -1; //file cannot be accessed
        return icod;
    }

    vector<double> input_v;

    if (read_file(ifs, input_v)) {
        icod = -2; //cannot read file
        return icod;
    }

    OutStruct out;
    icod = init_calc(input_v, out);
    ofstream ofs(O_FILE_PATH);
    ofs << out;
    return icod;
}
