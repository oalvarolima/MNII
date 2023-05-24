#ifndef MTX_UTILS_HPP
#define MTX_UTILS_HPP

#include <Eigen/Dense>

typedef Eigen::MatrixXd Matrix;

struct eigenResult {
    Matrix vector;
    double value;
};

namespace mtxUtils { 

    struct LU {
        Matrix L, U;
    };

    LU decompLU(const Matrix& M);
}

#endif 
