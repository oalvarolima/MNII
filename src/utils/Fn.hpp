#pragma once

#include "tinyexpr.h"
#include <string>

class Fn
{
public:
    Fn(const std::string &expr);

    double eval(double x) const;
private:
    te_expr* compiledExpr;
    double* x;
};