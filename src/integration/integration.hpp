#pragma once

#include <cmath>
#include <iostream>
#include <vector>

#include "../utils/Fn.hpp"

struct itrInfo
{
    uint32_t itrNum;
    double Ivalue;
    double EPS;
};

// Newton-Cotes
namespace NC {
    namespace closed {
        double two(const Fn& f, double upper, double lower);
        double three(const Fn& f, double upper, double lower);
        double four(const Fn& f, double upper, double lower);
        double five(const Fn& f, double upper, double lower);
    }

    namespace opened {
        double two(const Fn& f, double upper, double lower);
        double three(const Fn& f, double upper, double lower);
        double four(const Fn &f, double upper, double lower);
        double five(const Fn &f, double upper, double lower);
    }
}

// Quadraturas de Gauss
namespace Gauss {
    namespace Legendre {
        double two(const Fn &f, double upper, double lower);
        double three(const Fn &f, double upper, double lower);
        double four(const Fn &f, double upper, double lower);
    }

    namespace Hermite {
        double two(std::string f_inicial);
        double three(std::string f_inicial);
        double four(std::string f_inicial);
    }

    namespace Laguerre {
        double two(std::string f_inicial);
        double three(std::string f_inicial);
        double four(std::string f_inicial);
    }

    namespace Chebyshev {
        double two(std::string f_inicial);
        double three(std::string f_inicial);
        double four(std::string f_inicial);
    }
}

std::vector<itrInfo> integrate(const Fn& f, double upper, double lower, double EPS,
                               double(intervalIntegrator)(const Fn& f, double, double));

void printResults(const std::vector<itrInfo> &results, const std::string &methodName);
