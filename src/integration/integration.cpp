#include "integration.hpp"

double trapezoidalRule(const Fn &f,  double upper, double lower, uint32_t intervalsNum)
{
    double interval = (upper - lower) / intervalsNum;
    double result = f.eval(lower) + f.eval(upper);
    for( uint32_t i = 1; i < intervalsNum; i++ )
    {
        result += 2*f.eval(lower + i*interval);
    }
    return interval*.5*result;
}

double simpson13(const Fn &f,  double upper, double lower, uint32_t intervalsNum)
{
    double interval = (upper - lower) / intervalsNum;
    double result = f.eval(lower) + f.eval(upper) + 4*f.eval(lower + interval);
    for ( uint32_t i = 2 ; i < intervalsNum ; i += 2) 
    {
        result += 2*f.eval(lower + i*interval) + 4*f.eval(lower + (i+1)*interval);
    }

    return interval*(1./3)*result;
}

double simpson38(const Fn &f, double upper, double lower, uint32_t intervalsNum)
{
    double interval = (upper - lower) / intervalsNum;
    double result = f.eval(lower) + f.eval(upper) + 3*f.eval(lower + interval)+ 3*f.eval(lower + 2*interval);
    for ( uint32_t i = 3 ; i < intervalsNum ; i += 3) 
    {
        result += 2*f.eval(lower + i*interval) + 3*f.eval(lower + (i+1)*interval) + 3*f.eval(lower + (i+2)*interval);
    }

    return interval*(3./8)*result;
}

std::vector<itrInfo> integrate(const Fn &f, double upper, double lower, double EPS,
                        double (*integrator)(const Fn &, double, double, uint32_t),
                        uint32_t max_itr, uint32_t offSet)
{
    std::vector<itrInfo> results;
    double currError = EPS + 69;
    double oldResult = 0, result;
    for( uint32_t i = offSet; currError > EPS && (i/offSet) <= max_itr; i += offSet)
    {
        result = integrator(f, upper, lower, i);
        currError = fabs( (oldResult - result) / result );
        oldResult = result;

        results.push_back( {i, result, currError} );
    }

    return results;
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
