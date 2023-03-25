#include "Fn.hpp"

Fn::Fn(const std::string &expr)
{
    x = new double;
    te_variable variables[] = {"x", x};

    const char *charsExpr = expr.c_str();
    int errorCode;
    compiledExpr = te_compile(charsExpr, variables, 1, &errorCode);

    if( !compiledExpr )
    {
        printf("%s\n", charsExpr);
        printf("\t%*s^\nError near here", errorCode - 1, "");
    }
}

double Fn::eval(double x) const
{
    *this->x = x;
    return te_eval(compiledExpr);
}