#include "mtxUtils.hpp"

mtxUtils::LU mtxUtils::decompLU(const Matrix& M) {
    Matrix U = M, L = Matrix::Identity(M.rows(), M.cols());
 
    for (uint32_t col = 0; col < U.cols(); col++) {
        for (uint32_t row = col + 1; row < U.rows(); row++) {
            L(row, col) = U(row, col) /  U(col, col);
            U.row(row) -= L(row, col)*U.row(col);
        }
    }

    return {L, U};
}
