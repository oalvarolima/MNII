#include "src/NumMethods.hpp"
#include <iomanip>

int main() {
    std::cout << std::fixed << std::setprecision(9);
    LOG("\n=-=-=-=-=-=-=-=-=- Primeira matriz -=-=-=-=-=-=-=-=-=\n");
    const double EPS = .000001;
    Matrix primeira(3, 3);
    primeira << 5, 2, 1,
                2, 3, 1,
                1, 1, 2;
    auto result = PowerMethod::regular(primeira, EPS);
    auto inverse_result = PowerMethod::inverse(primeira, EPS);
    auto shifted_result = PowerMethod::shifted(primeira, EPS, 2);
    result.print();
    inverse_result.print();
    shifted_result.print();

    LOG("\n=-=-=-=-=-=-=-=-=- Segunda matriz -=-=-=-=-=-=-=-=-=\n");
    Matrix segunda(3, 3);
    segunda << -14,  1,  -2,
                 1, -1,   1,
                -2,  1, -11;
    auto segunda_result = PowerMethod::regular(segunda, EPS);
    auto segunda_inverse_result = PowerMethod::inverse(segunda, EPS);
    auto segunda_shifted_result = PowerMethod::shifted(segunda, EPS, -10);
    segunda_result.print();
    segunda_inverse_result.print();
    segunda_shifted_result.print();

    LOG("\n=-=-=-=-=-=-=-=-=- Terceira matriz -=-=-=-=-=-=-=-=-=\n");
    Matrix terceira(5, 5);
    terceira << 40, 8, 4, 2, 1, 8, 30, 12, 6, 2, 4, 12, 20, 1, 2, 2, 6, 1, 25, 4, 1, 2, 2, 4, 5;
    auto terceira_result = PowerMethod::regular(terceira, EPS);
    auto terceira_inverse_result = PowerMethod::inverse(terceira, EPS);
    auto terceira_shifted_result = PowerMethod::shifted(terceira, EPS, -4);
    auto terceira_shifted_result_I = PowerMethod::shifted(terceira, EPS, 33);
    auto terceira_shifted_result_II = PowerMethod::shifted(terceira, EPS, 22);
    auto terceira_shifted_result_III = PowerMethod::shifted(terceira, EPS, 11);
    terceira_result.print();
    terceira_shifted_result_I.print();
    terceira_shifted_result_II.print();
    terceira_shifted_result_III.print();
    terceira_inverse_result.print();
}