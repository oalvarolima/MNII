#include "src/NumMethods.hpp"

#include <iostream>

#define LOG(x) std::cout << x << std::endl;

int main()
{
    std::string t1("sqrt(e^(3*x) + 4*(x^2))");
    Fn f(t1);
    LOG(secondDeriv_ctr_E4(f, 2, .0001));

}