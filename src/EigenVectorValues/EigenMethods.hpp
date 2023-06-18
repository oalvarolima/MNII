#pragma once

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