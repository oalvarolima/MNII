#pragma once

#include "../utils/Fn.hpp"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>


struct itrInfo
{
    uint32_t itrNum;
    double Ivalue;
    double EPS;
};

// Fórmulas fechadas
double trapezoidal(const Fn& f, double upper, double lower);
double simpson13(const Fn& f, double upper, double lower);
double simpson38(const Fn& f, double upper, double lower);
double boole(const Fn& f, double upper, double lower);

// Fórmulas abertas
double opened2Points(const Fn& f, double upper, double lower);
double milne(const Fn& f, double upper, double lower);
double opened4points(const Fn &f, double upper, double lower);
double opened5points(const Fn &f, double upper, double lower);

std::vector<itrInfo> integrate(const Fn& f, double upper, double lower, double EPS,
                               double(intervalIntegrator)(const Fn& f, double, double));

void printResults(const std::vector<itrInfo> &results, const std::string &methodName);
