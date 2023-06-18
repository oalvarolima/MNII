#include "src/NumMethods.hpp"

int main() {
    Derivator a(Derivator::FORWARD, Derivator::SECOND, Derivator::E3);
    Derivator b(Derivator::BACKWARD, Derivator::SECOND, Derivator::E3);
    Derivator c(Derivator::CENTRAL, Derivator::SECOND, Derivator::E1);
    LOG(a.derivate(Fn("x^2"), 1, 0.00001));
    LOG(b.derivate(Fn("x^2"), 1, 0.00001));
    LOG(c.derivate(Fn("x^2"), 1, 0.00001));
}