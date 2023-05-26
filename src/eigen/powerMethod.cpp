#include "powerMethod.hpp"

eigenResult powerMethod::regular(const Matrix& A, double tolerance) {
    Matrix vk(A.rows(), 1), vkOld;
    vk.fill(1);

    double oldEigenVal, eigenVal = 0;
    double error = 1;
    while (error > tolerance) {
        oldEigenVal = eigenVal;
        vkOld = vk;
        vkOld.col(0).normalize();
        vk = A * vkOld;
        eigenVal = vkOld.col(0).dot(vk.col(0));
        error = std::fabs((eigenVal - oldEigenVal) / eigenVal);
    }

    return {vkOld, eigenVal};
}

eigenResult powerMethod::inverse(const Matrix& A, double tolerance) {
    eigenResult result = regular(A.inverse(), tolerance);
    return {result.vector, 1. / result.value};
}

eigenResult powerMethod::shifted(const Matrix& A, double shiftment, double tolerance) {
    eigenResult result = inverse(A - (shiftment*Matrix::Identity(A.rows(), A.cols())), tolerance);
    return {result.vector, result.value + shiftment};
}

void powerMethod::print(const eigenResult &result) {
    std::cout << "\nAuto-valor: " << result.value;
    std::cout << "\nAuto-vetor correspondente:\n" << result.vector << '\n';
}
