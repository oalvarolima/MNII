#include "houseHolderMethod.hpp"

HouseHolder::result HouseHolder::makeTridiagMatrix(const Matrix& A) {
    Matrix accumHHs = Matrix::Identity(A.rows(), A.cols());
    Matrix tridiagMatrix = A;
    
    for(uint32_t i = 0; i < A.rows() - 2; i++) {
        Matrix HHMatrix = makeHHcolI(tridiagMatrix, i);
        tridiagMatrix = HHMatrix.transpose() * tridiagMatrix * HHMatrix;
        accumHHs *= HHMatrix;
    }

    return {tridiagMatrix, accumHHs};
}

Matrix HouseHolder::makeHHcolI(const Matrix& A, uint32_t i) {
    Matrix w(A.rows(), 1);
    w.fill(0);

    Matrix wl(A.rows(), 1);
    wl.fill(0);

    for(uint32_t d = i + 1; d < A.rows(); d++) {
        w(d, 0) = A(d, i);
    }

    wl(i+1, 0) = w.norm();
    Matrix n = (w - wl).col(0).normalized();

    return Matrix::Identity(A.rows(), A.cols()) - (2 * (n * n.transpose()));
}
