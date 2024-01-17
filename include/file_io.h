#pragma once

#include <istream>
#include <ostream>

#include <vector>


/*
 * @return false if fail
 * @param i_stream - file stream
 * @param res - vector to write
 */
bool read_file(std::istream &i_stream, std::vector<double> &res)
{
    const int param_count = 6;
    res = std::vector<double>(param_count);
    for (int i = 0; i < param_count; ++i)
    {
        i_stream >> res[i];
    }
    return !i_stream.fail();
}

struct OutStruct
{
    // OutStruct() : er(0), h(0), icod(-3), x_l(0), y_l(0) {};
    double er; // абсолютная погрешность
    double h;  // шаг
    int icod;
    double x_l; // последний x
    double y_l; // последний y
};

std::ostream &operator<<(std::ostream &ofs, const OutStruct &e)
{
    const char sep = ' ';
    ofs << e.er << sep << e.h << sep << e.icod << std::endl;
    ofs << e.x_l << sep << e.y_l << std::endl;
    return ofs;
}

/*
 * @return true if fail
 * @param out_stream - file stream
 * @param data - data to write
 */
bool out_file(std::ostream &out_stream, OutStruct &data)
{
    out_stream << data;
    return !out_stream.fail();
}