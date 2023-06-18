#include "EigenMethods.hpp"

QRMethodResult QRMethod(const Matrix &A, double tolerance) {
    Matrix result = A;
    Matrix accumQs = Matrix::Identity(A.rows(), A.cols());
    do {
        QRdecomp qr = QR::decomp(result);
        result = qr.R * qr.Q;
        accumQs *= qr.Q;
    } while(!isDiagonal(result, tolerance));

    return {result, accumQs};
}

JacobiMethodResult JacobiMethod(const Matrix &A, double tolerance) {
    Matrix result = A;
    Matrix accumJs = Matrix::Identity(A.rows(), A.cols());
    do {
        JacobiScanResult jacobi = Jacobi::transform(result);
        result = jacobi.diagonalMatrix;
        accumJs *= jacobi.accumJ;
    } while(!isDiagonal(result, tolerance));

    return {result, accumJs};
}