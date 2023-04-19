#include "integration.hpp"

double integratorFunction(const Fn& f, double upper, double lower, uint32_t intervalsNum, 
                          double(intervalIntegrator)(const Fn&, double, double))
{
    double interval = (upper - lower) / intervalsNum;
    double result = 0;
    for( uint32_t n = 0; n < intervalsNum; n++, lower += interval )
    {
        result += intervalIntegrator(f, lower + interval, lower);
    }

    return result;
}

std::vector<itrInfo> integrate(const Fn& f, double upper, double lower, double EPS,
                               double(intervalIntegrator)(const Fn& f, double, double))
{
    std::vector<itrInfo> results;
    const uint32_t MaxItr = 100000;
    double Error = EPS + 69.696969;
    double result, oldResult = 0;
    for( uint32_t intervalsNum = 1; Error >= EPS && intervalsNum <= MaxItr; intervalsNum++ )
    {
        result = integratorFunction(f, upper, lower, intervalsNum, intervalIntegrator);
        Error = fabs( (result - oldResult) / result );
        
        results.push_back( {intervalsNum, result, Error} );

        oldResult = result;
    }

    return results;
}

double trapezoidal(const Fn& f, double upper, double lower)
{
    return .5*(upper - lower)*(f.eval(upper) + f.eval(lower));
}

double simpson13(const Fn& f, double upper, double lower)
{
    double h = (upper - lower) / 2;
    return (1./3)*h*(f.eval(lower) + 4*f.eval(lower + h) + f.eval(upper));
}

double simpson38(const Fn& f, double upper, double lower)
{
    double h = (upper - lower) / 3;
    return (3./8)*h*(f.eval(lower) + 3*f.eval(lower + h) + 3*f.eval(lower + 2*h) + f.eval(upper));
}

double boole(const Fn& f, double upper, double lower)
{
    double h = (upper - lower) / 4;
    return (2./45)*h*(7*f.eval(lower) + 32*f.eval(lower + h) + 12*f.eval(lower + 2*h) + 32*f.eval(lower + 3*h) + 7*f.eval(upper));
}

double opened2Points(const Fn &f, double upper, double lower)
{
    double interval = (upper - lower) / 3;
    return (3./2)*interval*(f.eval(lower + interval) + f.eval(lower + 2*interval));
}

double milne(const Fn &f, double upper, double lower)
{
    double interval = (upper - lower) / 4;
    return (4./3)*interval*(2*f.eval(lower + interval) - f.eval(lower + 2*interval) + 2*(f.eval(lower + 3*interval)));
}

double opened4points(const Fn &f, double upper, double lower)
{
    double interval = (upper - lower) / 5;
    return (5./24)*interval*(11*f.eval(lower + interval) + f.eval(lower + 2*interval) + (f.eval(lower + 3*interval)) + 11*(f.eval(lower + 4*interval)));
}


double opened5points(const Fn &f, double upper, double lower)
{
    double interval = (upper - lower) / 6;
    return (3./20)*interval*(22*f.eval(lower + interval) - 28*f.eval(lower + 2*interval) + 52*(f.eval(lower + 3*interval)) - 28*(f.eval(lower + 4*interval)) + 22*(f.eval(lower + 5*interval)));
}

void printResults(const std::vector<itrInfo> &results, const std::string &methodName)
{
    printf("\n                            %s\n\n", methodName.c_str());
    printf("N de intervalos               Integral                    EPS\n");
    for(auto [itrNum, Ivalue, EPS] : results)
    {
        printf("%s%u                     %s%.8f%s%.8f\n",
               std::string(8 - std::to_string(itrNum).length(), ' ').c_str(),
               itrNum,
               Ivalue < 0 ? "" : " ",
               Ivalue,
               std::string(26 - std::to_string(Ivalue).length(), ' ').c_str(),
               EPS);
    }
}
