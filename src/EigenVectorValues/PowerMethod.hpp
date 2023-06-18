#pragma once

#include <iostream>

#include "../utils/Matrix.hpp"

class PowerMethod {
    struct result {
        Matrix vector;
        double value;

        void print();
    };
public:
    static PowerMethod::result regular(const Matrix &A, double tolerance);
    static PowerMethod::result inverse(const Matrix &A, double tolerance);
    static PowerMethod::result shifted(const Matrix &A, double shiftment, double tolerance);
};
