#include "SimilarityTransformer.hpp"

LUdecomp LU::decomp(const Matrix &A) {
    Matrix L = Matrix::Identity(A.rows(), A.cols()), U = A;

    for (uint32_t col = 0; col < U.cols(); col++)
        for (uint32_t row = col + 1; row < U.rows(); row++) {
            L(row, col) = U(row, col) /  U(col, col);
            U.row(row) -= L(row, col)*U.row(col);
        }

    return {L, U};
}

Matrix LU::inverse(const Matrix &A) {
    Matrix inverse = Matrix::Zero(A.rows(), A.cols());
    LUdecomp lu = decomp(A);
    Matrix ID = Matrix::Identity(A.rows(), A.cols());
    for(uint32_t i = 0; i < A.cols(); i++)
        inverse.col(i) = solver(lu, ID.col(i));

    return inverse;
}

Matrix LU::solver(const LUdecomp &lu, const Matrix &y) {
    Matrix result = y;
    for(uint32_t row = 0; row < y.rows(); row++) {
        double sum = 0;
        for(uint32_t col = 0; col < row; col++)
            sum += lu.L(row, col) * result(col, 0);

        result(row, 0) = y(row, 0) - sum;
    }
    for(int row = result.rows() - 1; row >= 0; row--) {
        double sum = 0;
        for(int col = result.rows() - 1; col > row; col--)
            sum += lu.U(row, col) * result(col, 0);

        result(row, 0) = (result(row, 0) - sum) / lu.U(row, row);
    }
    return result;
}

QRdecomp QR::decomp(const Matrix &A) {
    Matrix transpQ = Matrix::Identity(A.rows(), A.cols());
    Matrix R = A;

    for(int col = 0; col < A.cols() - 1; col++)
        for(int row = col + 1; row < A.rows(); row++) {
            Matrix J = makeOldHHMatrix(R, col, row);
            R = J * R;
            transpQ = J * transpQ;
        }

    return {transpQ.transpose(), R};
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

HHtransformed HouseHolder::transform(const Matrix &A) {
    Matrix tridiagMatrix = A;
    Matrix accumHH = Matrix::Identity(A.rows(), A.cols());

    for(uint32_t i = 0; i < A.rows() - 2; i++) {
        Matrix HHMatrix = HHMatrixToColI(tridiagMatrix, i);
        tridiagMatrix = HHMatrix.transpose() * tridiagMatrix * HHMatrix;
        accumHH *= HHMatrix;
    }

    return {tridiagMatrix, accumHH};
}

Matrix HouseHolder::HHMatrixToColI(const Matrix &A, uint32_t i) {
    Matrix w = A.col(i);
    Matrix wl = Matrix::Zero(A.rows(), 1);

    wl(i+1, 0) = w.norm();
    Matrix n = (w - wl).col(0).normalized();

    return Matrix::Identity(A.rows(), A.cols()) - (2 * (n * n.transpose()));
}

JacobiScanResult Jacobi::transform(const Matrix &A) {
    Matrix newA = A;
    Matrix accumJ = Matrix::Identity(A.rows(), A.cols());
    for (int col = 0; col < A.cols() - 1; col++)
        for (int row = col + 1; row < A.rows(); row++) {
            Matrix J = jacobiMatrixIJ(newA, row, col);
            newA = J.transpose() * newA * J;
            accumJ *= J;
        }

    return {newA, accumJ};
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

//TODO
SVDdecomp SVD::decomp(const Matrix &A) {
    Matrix DUMMY = Matrix::Random(1, 1);
    return {DUMMY, DUMMY, DUMMY};
}

bool isDiagonal(const Matrix &A, double EPS) {
    double sum = 0;
    for (int row = 0; row < A.rows(); row++)
        for (int col = 0; col < row; col++)
            sum += A(row, col) * A(row, col);

    return abs(sum) < EPS;
}
