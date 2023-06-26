#include "src/NumMethods.hpp"

int main() {
    Matrix A(3, 3);
    A << 1, 1, -1,
         1, -2, 3,
         2, 3, 1;
    LOG('\n' << A.inverse() << "\n");
    LOG(LU::inverse(A));
}