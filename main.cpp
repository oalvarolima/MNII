#include "src/NumMethods.hpp"

int main()
{
    Fn f("sqrt(e^(3*x) + 4*(x^2))");
    LOG(secondDeriv_ctr_E4(f, 2, .0001));
}