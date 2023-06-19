#include "src/NumMethods.hpp"

int main() {
    Derivator a(Derivator::FORWARD, Derivator::THIRD, Derivator::E2);
    Derivator b(Derivator::BACKWARD, Derivator::THIRD, Derivator::E2);
    Derivator c(Derivator::CENTRAL, Derivator::THIRD, Derivator::E1);
    LOG("Forward segunda E3: " << a.derivate(Fn("x^4"), 1, 0.00001));
    LOG("Backward segunda E3: " << b.derivate(Fn("x^4"), 1, 0.00001));
    LOG("Central segunda E1: " << c.derivate(Fn("x^4"), 1, 0.00001));
}