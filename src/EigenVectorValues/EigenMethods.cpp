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

PowerMethod::result PowerMethod::regular(const Matrix &A, double tolerance) {
    Matrix currV = Matrix::Random(A.rows(), 1);
    Matrix prevV;
    double currLambda = 0, prevLambda, error = 1;
    while(error > tolerance) {
        prevV = currV;
        prevLambda = currLambda;
        currV = A * prevV.col(0).normalized();
        currLambda = currV.col(0).dot(prevV.col(0));
        error = std::fabs((currLambda - prevLambda) / currLambda);
    }

    return {prevV, currLambda};
}

PowerMethod::result PowerMethod::inverse(const Matrix &A, double tolerance) {
    result regularResult = regular(A.inverse(), tolerance);
    return {regularResult.vector, 1/regularResult.value};
}

PowerMethod::result PowerMethod::shifted(const Matrix &A, double shiftment, double tolerance) {
    result inverseResult = inverse(A - (shiftment*Matrix::Identity(A.rows(), A.cols())), tolerance);
    return {inverseResult.vector, inverseResult.value + shiftment};
}

void PowerMethod::result::print() {
    std::cout << "\nAuto-valor: " << value << "Auto-vetor:\n" << vector << "\n";
}