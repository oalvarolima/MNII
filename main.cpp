#include "src/NumMethods.hpp"

int main() {
    Fn f("(sin(2 * x) + 4*(x^2) + 3*x)^2");
    printResults(integrate(f, 1, 0, .000001, NC::opened::two), "Newton-cotes 5 pontos");
    printResults(integrate(f, 1, 0, .000001, NC::opened::three), "Newton-cotes 3 pontos");
    printResults(integrate(f, 1, 0, .000001, Gauss::Legendre::two), "Gauss-Legendre 2 pontos");
    printResults(integrate(f, 1, 0, .000001, Gauss::Legendre::three), "Gauss-Legendre 3 pontos");
    printResults(integrate(f, 1, 0, .000001, Gauss::Legendre::four), "Gauss-Legendre 4 pontos");
}