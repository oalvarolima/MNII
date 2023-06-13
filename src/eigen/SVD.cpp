#include "SVD.hpp"

SVD::result SVD::decompSVD(const Matrix &A) {
    Matrix  bSl = A * A.transpose(),
            bSr = A.transpose() * A;

    powerMethod::allEigenResults leastSMatrix = A.rows() < A.cols() ?
            powerMethod::getAllEigenResults(bSl, .000001) :
            powerMethod::getAllEigenResults(bSr, .000001);

    LOG(leastSMatrix.eigenValues);
    LOG(leastSMatrix.eigenVectors);

    return {bSl, bSl, bSl};
}