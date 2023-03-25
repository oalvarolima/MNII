#pragma once

#include "../utils/Fn.hpp"

// Calcula a segunda de derivada de f com ordem de erro O(h^4) usando a filosfia central
double secondDeriv_ctr_E4(const Fn &f, double x, double interval);