#include "src/NumMethods.hpp"

int main() {
    LOG("\nexponencial simples -> " << simpleExponential(Fn("1 / cbrt(x^2)"), 1, -1, 4));
    LOG("exponencial dupla -> " << doubleExponential(Fn("1 / sqrt(4-x^2)"), 0, -2, 2.5));
}