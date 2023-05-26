#include "QRmethod.hpp"

QR::methodResult QR::method(const Matrix &A, double tolerance) {
    Matrix diagMatrix = A;
    Matrix accumQs = Matrix::Identity(A.rows(), A.cols());

    do {
        QandR qr = decomp(diagMatrix);
        diagMatrix = qr.R * qr.Q;
        accumQs *= qr.Q;
    } while (sumOfSquaredTermsUnderDiag(diagMatrix) > tolerance);

    return {diagMatrix.diagonal(), accumQs};
}

QR::QandR QR::decomp(const Matrix& A) {
    Matrix R = A;
    Matrix Qtransposed = Matrix::Identity(A.rows(), A.cols());

    for(int col = 0; col < A.cols() - 1; col++) {
        for(int row = col + 1; row < A.rows(); row++) {
            Matrix J = makeOldHHMatrix(R, col, row);
            R = J * R;
            Qtransposed = J * Qtransposed;
        }
    }

    return {Qtransposed.transpose(), R};
}

Matrix QR::makeOldHHMatrix(const Matrix &A, int col, int row) {
    Matrix J = Matrix::Identity(A.rows(), A.cols());

    if(abs(A(row, col)) <= .000001)
        return J;

    double theta;
    if(abs(A(col, col)) <= .000001)
        theta = A(row, col) < 0 ? M_PI / 2. : -M_PI / 2.;
    else
        theta = atan(-A(row, col) / A(col, col));

    J(row, row) = J(col, col) = cos(theta);
    J(row, col) = sin(theta);
    J(col, row) = -J(row, col);

    return J;
}