#pragma once

#include "../utils/Fn.hpp"

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

// Gauss Legendre
double gauss_legendre_2_points(const Fn &f, double upper, double lower);
double gauss_legendre_3_points(const Fn &f, double upper, double lower);
double gauss_legendre_4_points(const Fn &f, double upper, double lower);

// Gauss Hermite
double gauss_hermite_2_points(std::string f_inicial);
double gauss_hermite_3_points(std::string f_inicial);
double gauss_hermite_4_points(std::string f_inicial);

// Gauss Laguerre
double gauss_laguerre_2_points(std::string f_inicial);
double gauss_laguerre_3_points(std::string f_inicial);
double gauss_laguerre_4_points(std::string f_inicial);

// Gauss Chebyshev
double gauss_chebyshev_2_points(std::string f_inicial);
double gauss_chebyshev_3_points(std::string f_inicial);
double gauss_chebyshev_4_points(std::string f_inicial);

std::vector<itrInfo> integrate(const Fn& f, double upper, double lower, double EPS,
                               double(intervalIntegrator)(const Fn& f, double, double));

void printResults(const std::vector<itrInfo> &results, const std::string &methodName);
