#include "Fn.hpp"

double logn(double x, double base) {
    return log(x)/log(base);
}

Fn::Fn(const std::string &expr) : expression(expr) {
    x = new double;
    te_variable variables[] = {
            {"logn", (const void *)logn, TE_FUNCTION2},
            {"x",    x}
    };

    const char *charsExpr = expr.c_str();
    int errorCode;
    compiledExpr = te_compile(charsExpr, variables, 2, &errorCode);

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

Fn::~Fn() {
    te_free(compiledExpr);
    delete x;
}
