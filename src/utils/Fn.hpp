#pragma once

#include "tinyexpr.h"
#include <string>
#include <iostream>
#include <math.h>

#define LOG(x) std::cout << x << std::endl

class Fn
{
public:
    Fn(const std::string &expr);
    ~Fn();

    double eval(double x) const;
    void test(double x) const;
    
private:
    std::string expre;
    te_expr* compiledExpr;
    double* x;
};
