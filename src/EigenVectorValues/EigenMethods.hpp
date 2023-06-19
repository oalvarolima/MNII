#pragma once

#include <iostream>

#include "SimilarityTransformer.hpp"

struct QRMethodResult {
    Matrix eigenValues;
    Matrix accumQs;
};
QRMethodResult QRMethod(const Matrix& A, double tolerance);

struct JacobiMethodResult {
    Matrix eigenValues;
    Matrix accumJs;
};
JacobiMethodResult JacobiMethod(const Matrix& A, double tolerance);

class PowerMethod {
    struct result {
        Matrix vector;
        double value;

        void print();
    };
public:
    static result regular(const Matrix &A, double tolerance);
    static result inverse(const Matrix &A, double tolerance);
    static result shifted(const Matrix &A, double shiftment, double tolerance);
};