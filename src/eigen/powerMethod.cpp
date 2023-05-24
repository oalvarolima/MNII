#include "powerMethod.hpp"

eigenResult powerMethod::regular(const Matrix& A, double tolerance) {
    Matrix vk(A.rows(), 1), vkOld;
    vk.fill(1);

    double eigenValOld, eigenVal = 0;
    double error = tolerance+tolerance;
    while (error > tolerance) {
        eigenValOld = eigenVal;
        vkOld = vk;
        vkOld.col(0).normalize();
        vk = A * vkOld;
        eigenVal = vkOld.col(0).dot(vk.col(0));
        error = std::fabs( (eigenVal - eigenValOld) / eigenVal);
    }

    return {vkOld, eigenVal};
}

#define LOG(x) std::cout << x << std::endl;

eigenResult powerMethod::inverse(const Matrix& A, double tolerance) {
    eigenResult a = regular(A.inverse(), tolerance);
    return {a.vector, 1./a.value};
}

eigenResult powerMethod::shifted(const Matrix& A, double shiftment, double tolerance) {
    eigenResult result = inverse(A - (shiftment*Matrix::Identity(A.rows(), A.cols())), tolerance);
    return {result.vector, result.value + shiftment};
}

void powerMethod::print(eigenResult a) {
    std::cout << "\nAuto-valor: " << a.value;
    std::cout << "\nAuto-vetor correspondente:\n" << a.vector << '\n';
}
