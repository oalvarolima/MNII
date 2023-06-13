#pragma once

#include "mtxUtils.hpp"

namespace Jacobi {
    struct scanResult {
        Matrix diagMatrix;
        Matrix accumJ;
    };

    struct methodResult {
        Matrix eigenValues;
        Matrix accumJs;
    };

    methodResult makeDiagMtx(const Matrix &A, double tolerance);
    scanResult jacobiScan(const Matrix &A);
    Matrix jacobiMatrixIJ(const Matrix &A, int row, int col);
}
