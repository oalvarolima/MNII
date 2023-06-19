#pragma once

#include "tinyexpr.h"
#include <string>
#include <iostream>
#include <cmath>

#define LOG(x) std::cout << x << std::endl

class Fn {
public:
    Fn(const std::string &expr);

    double eval(double x) const;
    void test(double x) const;

    std::string expression;

private:
    te_expr* compiledExpr;
    double* x;
};