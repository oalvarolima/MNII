#pragma once

#include "../utils/Matrix.hpp"

struct LUdecomp {
    Matrix L, U;
};
class LU {
public:
    static LUdecomp decomp(const Matrix& A);
    static Matrix inverse(const Matrix& A);
    static Matrix solver(const LUdecomp& lu, const Matrix& y);
};

struct QRdecomp {
    Matrix Q, R;
};
class QR {
public:
    static QRdecomp decomp(const Matrix& A);
private:
    static Matrix makeOldHHMatrix(const Matrix &A, int col, int row);
};

struct HHtransformed {
    Matrix tridiagMatrix;
    Matrix accumHHs;
};
class HouseHolder {
public:
    static HHtransformed transform(const Matrix& A);
private:
    static Matrix HHMatrixToColI(const Matrix &A, uint32_t i);
};

struct JacobiScanResult {
    Matrix diagonalMatrix;
    Matrix accumJ;
};
class Jacobi {
public:
    static JacobiScanResult transform(const Matrix& A);
private:
    static Matrix jacobiMatrixIJ(const Matrix &A, int row, int col);
};

//TODO
struct SVDdecomp {
    Matrix U, S, V;
};
class SVD {
public:
    static SVDdecomp decomp(const Matrix& A);
};

bool isDiagonal(const Matrix &A, double EPS);