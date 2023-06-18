#pragma once

#include <vector>

#include "../utils/Fn.hpp"

class SpecialIntegrator {
public:
    enum Type {
        HERMITE = 0,
        LAGUERRE = 1,
        CHEBYSHEV = 2
    };

    enum Degree {
        TWO = 0,
        THREE = 1,
        FOUR = 2,
    };

    SpecialIntegrator(Type type, Degree degree);

    double integrate(const Fn &f);

    struct TableValues {
        std::vector<double> weights;
        std::vector<double> alfa;
    };

private:
    const TableValues tableValues;
    const std::string weightFunction;
    const std::string functionToAdd;
};