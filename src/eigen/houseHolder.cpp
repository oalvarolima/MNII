#include "houseHolder.hpp"

Matrix makeHHcolI(const Matrix& A, uint32_t i);

HH::result HH::makeHHMatrix(const Matrix& A) {
    Matrix H(A.rows(), A.cols()), Hi, Ai; 
    Matrix Am1 = A;
    
    for(uint32_t i = 0; i < A.rows() - 2; i++) {
        Hi = makeHHcolI(Am1, i);
        Ai = Hi.transpose() * Am1 * Hi;
        Am1 = Ai;
        H = H*Hi;
    }

    return {Ai, H};
}


Matrix makeHHcolI(const Matrix& A, uint32_t i) {
    Matrix w(A.rows(), 1), wl(A.rows(), 1);
    w.fill(0); wl.fill(0);

    for(uint32_t d = i + 1; d < A.rows(); d++) {
        w(d, 0) = A(d, i);
    }

    wl(i+1, 0) = w.norm();
    Matrix n = w - wl;
    n.col(0).normalize();

    return Matrix::Identity(A.rows(), A.cols()) - (2 * (n * n.transpose()));
}
