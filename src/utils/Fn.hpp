#pragma once

#include "tinyexpr.h"
#include <string>
#include <iostream>

#define LOG(x) std::cout << x << std::endl;

class Fn
{
public:
    Fn(const std::string &expr);

    double eval(double x) const;
private:
    te_expr* compiledExpr;
    double* x;
};
