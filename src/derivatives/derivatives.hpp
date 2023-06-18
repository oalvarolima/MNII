#pragma once

#include <vector>

#include "../utils/Fn.hpp"

class Derivator {
public:
    enum Order {
        FIRST = 0,
        SECOND = 1,
        THIRD = 2,
    };

    enum Type {
        FORWARD = 0,
        BACKWARD = 0,
        CENTRAL = 1,
    };

    enum Degree {
        E1 = 0,
        E2 = 1,
        E3 = 2,
        E4 = 3,
    };

    Derivator(Type type, Order order, Degree degree);

    double derivate(const Fn& fn, double x, double interval) const;

    struct Table {
        double multiplier;
        std::vector<double> weights;
        int intervalExponent;
    };

private:
    Table values;
    Type type;
};