#include "Fn.hpp"

double logn(double x, double base) {
    return log(x)/log(base);
}

double cbrt_(double x) {
    return std::cbrt(x);
}

Fn::Fn(const std::string &expr) : expression(expr) {
    x = new double;
    te_variable variables[] = {
            {"logn", (const void *)logn,  TE_FUNCTION2},
            {"cbrt", (const void *) cbrt_, TE_FUNCTION1},
            {"x",    x}
    };

    const char *charsExpr = expr.c_str();
    int errorCode;
    compiledExpr = te_compile(charsExpr, variables, 3, &errorCode);

    if(!compiledExpr) {
        printf("\n\t%s\n", charsExpr);
        printf("\t%*s^\nError near here", errorCode - 1, "");
    }
}

double Fn::eval(double x) const {
    *this->x = x;
    return te_eval(compiledExpr);
}

void Fn::test(double x) const {
    std::cout << "\nf(x) = " << expression << "\n";
    std::cout << "f(" << x << ") = " << eval(x) << std::endl;
}