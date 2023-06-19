#include "src/NumMethods.hpp"

int main() {
    LOG("\n\texponencial simples -> " << simpleExponential(Function("1 / cbrt(x^2)"), 1, -1, 4));
    LOG("\texponencial dupla   -> " << doubleExponential(Function("1 / sqrt(4-x^2)"), 0, -2, 2.5));
}