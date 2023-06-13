#include "powerMethod.hpp"

eigenResult powerMethod::regular(const Matrix& A, double tolerance) {
    Matrix vk(A.rows(), 1), vkOld;
    vk.fill(1);
    vk(0, 0) = 0;

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

    return {vk.col(0).normalized(), eigenVal};
}

eigenResult powerMethod::inverse(const Matrix& A, double tolerance) {
    eigenResult result = regular(A.inverse(), tolerance);
    return {result.vector, 1. / result.value};
}

eigenResult powerMethod::shifted(const Matrix& A, double shiftment, double tolerance) {
    eigenResult result = inverse(A - (shiftment*Matrix::Identity(A.rows(), A.cols())), tolerance);
    return {result.vector, result.value + shiftment};
}

bool equals(double n1, double n2) {
    return std::fabs(n1 - n2) < .1;
}

powerMethod::allEigenResults powerMethod::getAllEigenResults(const Matrix &A, double tolerance) {
    Matrix eigenVectors(A.rows(), A.cols());
    Matrix eigenValues = Matrix::Identity(A.rows(), A.cols());

    eigenResult result = regular(A, tolerance);
    eigenVectors.col(0) = result.vector;
    eigenValues(0, 0) = result.value;

    result = inverse(A, tolerance);
    eigenVectors.col(A.cols() - 1) = result.vector;
    eigenValues(A.rows() - 1, A.cols() - 1) = result.value;

    double currEigenValue, newEigenValue, increment = 1.;
    currEigenValue = newEigenValue = eigenValues(0, 0);
    for (int i = 1; i < A.rows() - 1; ++i) {
        double shiftment = currEigenValue;
        while(equals(newEigenValue, currEigenValue)) {
            result = shifted(A, shiftment, tolerance);
            newEigenValue = result.value;
            shiftment -= increment;
        }

        eigenVectors.col(i) = result.vector;
        eigenValues(i, i) = result.value;
        currEigenValue = newEigenValue;
    }

    return {eigenVectors, eigenValues};
}

void powerMethod::print(const eigenResult &result) {
    std::cout << "\nAuto-valor: " << result.value;
    std::cout << "\nAuto-vetor correspondente:\n" << result.vector << '\n';
}
