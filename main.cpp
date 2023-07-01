#include "src/NumMethods.hpp"
#include <iomanip>

int main() {
    std::cout << std::fixed << std::setprecision(9);
    const double EPS = 1e-6;
    LOG("\n=-=-=-=-=-=-=-=-=- Potência -=-=-=-=-=-=-=-=-=\n");
    Matrix A(5, 5);
    A << 40, 8, 4, 2, 1, 8, 30, 12, 6, 2, 4, 12, 20, 1, 2, 2, 6, 1, 25, 4, 1, 2, 2, 4, 5;
    auto terceira_result = PowerMethod::regular(A, EPS);
    auto terceira_shifted_result_I = PowerMethod::shifted(A, EPS, 31);
    auto terceira_shifted_result_II = PowerMethod::shifted(A, EPS, 23);
    auto terceira_shifted_result_III = PowerMethod::shifted(A, EPS, 11);
    auto terceira_inverse_result = PowerMethod::inverse(A, EPS);
    terceira_result.print();
    terceira_shifted_result_I.print();
    terceira_shifted_result_II.print();
    terceira_shifted_result_III.print();
    terceira_inverse_result.print();

    LOG("\n=-=-=-=-=-=-=-=-=-=-=-=-=- metodo QR =-=-=-=-=-=-=-=-=-=-=-=-=-\n");
    auto qr_result = QRMethod(A, EPS);
    LOG(qr_result.accumQs);
    LOG("\n\n\tauto-valores: \n" << qr_result.eigenValues);

    LOG("\n=-=-=-=-=-=-=-=-=-=-=-=-=- metodo Jacobi =-=-=-=-=-=-=-=-=-=-=-=-=-\n");
    auto jacobi_result = JacobiMethod(A, EPS);
    LOG(jacobi_result.accumJs);
    LOG("\n\n\tauto-valores: \n" << jacobi_result.eigenValues);

    LOG("\n=-=-=-=-=-=-=-=-=-=-=-=-=- teste com QR e householder =-=-=-=-=-=-=-=-=-=-=-=-=-\n");
    auto hh_result = HouseHolder::transform(A);
    auto qr_result_II = QRMethod(hh_result.tridiagMatrix, EPS);
    LOG(qr_result_II.accumQs);
    LOG("\n\n\tP = HP\n" << hh_result.accumHHs*qr_result_II.accumQs);
}