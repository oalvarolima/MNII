#pragma once

#include "tinyexpr.h"
#include <iostream>

#define LOG(x) std::cout << x << std::endl

class Fn
{
public:
    Fn(const std::string &expr);

    double eval(double x);

private:
    te_expr* compiledExpr;
    double x;
};