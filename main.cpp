#include "src/NumMethods.hpp"

int main()
{
    Fn f("e^(-x^2)");
//    printResults(integrate(f, 1, 0, .000001, boole), "Newton-Cotes fechado grau 4");
//    printResults(integrate(f, 1, 0, .000001, opened5points), "Newton-Cotes aberto grau 4");
    runIntegralCalc();
}
