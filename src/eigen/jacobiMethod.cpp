#include "jacobiMethod.hpp"

Jacobi::methodResult Jacobi::makeDiagMtx(const Matrix &A, double tolerance) {
    Matrix diagMatrix = A, accumJs = Matrix::Identity(A.rows(), A.cols());
    do {
        scanResult result = jacobiScan(diagMatrix);
        accumJs *= result.accumJ;
        diagMatrix = result.diagMatrix;
    } while(sumOfSquaredTermsUnderDiag(diagMatrix) > tolerance);

    return {diagMatrix.diagonal(), accumJs};
}

Jacobi::scanResult Jacobi::jacobiScan(const Matrix& A) {
    Matrix Anew = A, accumJ = Matrix::Identity(A.rows(), A.cols());
    for (int col = 0; col < A.cols() - 1; col++) {
        for (int row = col + 1; row < A.rows(); row++) {
            Matrix J = Jacobi::jacobiMatrixIJ(Anew, row, col);
            Anew = J.transpose() * Anew * J;
            accumJ *= J;
        }
    }

    return {Anew, accumJ};
}

Matrix Jacobi::jacobiMatrixIJ(const Matrix &A, int row, int col) {
    Matrix J = Matrix::Identity(A.rows(), A.cols());
    if(abs(A(row, col)) <= .000001)
        return J;

    double theta;
    if(abs(A(row, row) - A(col, col)) <= .000001)
        theta = M_PI/4.;
    else
        theta = .5 * (atan((-2.*A(row, col)) / (A(row, row) - A(col, col))));

    J(row, row) = J(col, col) = cos(theta);
    J(row, col) = sin(theta);
    J(col, row) = -J(row, col);

    return J;
}