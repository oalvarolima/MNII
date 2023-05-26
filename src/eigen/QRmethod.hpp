#pragma once

#include "mtxUtils.hpp"

namespace QR {
    struct methodResult {
        Matrix eigenValues;
        Matrix accumQs;
    };
    methodResult method(const Matrix &A, double tolerance);

    struct QandR {
        Matrix Q, R;
    };
    QandR decomp(const Matrix& A);

    Matrix makeOldHHMatrix(const Matrix &A, int col, int row);
}