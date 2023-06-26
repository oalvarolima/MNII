#include "Integrator.hpp"

IntegratorIterations Integrator::integrate(const Function& f, double upper, double lower, double relativeErrorTolerance) {
    IntegratorIterations iterationsData;

    double currError = 1;
    double currIntegralValue, prevIntegralValue = 0;
    uint32_t currIntervalsNum = 1;
    while(currError >= relativeErrorTolerance && currIntervalsNum <= maxIntervalsNum) {
        currIntegralValue = integrateIntervals(f, currIntervalsNum, upper, lower);

        currError = std::abs((currIntegralValue - prevIntegralValue)/currIntegralValue);
        iterationsData.add(currIntervalsNum, currIntegralValue, currError);

        prevIntegralValue = currIntegralValue;
        currIntervalsNum = currIntervalsNum << 1;
    }

    return iterationsData;
}

double Integrator::integrateIntervals(const Function &f, uint32_t intervalsNum, double upper, double lower) {
    double intervalSize = (upper - lower)/intervalsNum;
    double currLower = lower;
    double currUpper = lower + intervalSize;
    double integralValue = 0;
    for(uint32_t i = 0; i < intervalsNum; i++) {
        integralValue += intervalIntegrator(f, currUpper, currLower);
        currLower = currUpper;
        currUpper += intervalSize;
    }

    return integralValue;
}

const std::vector<std::vector<NewtonCotes::TableValues>> NewtonCotesTable =
{
    {
        {1, {1, 1}, 1./2, NewtonCotes::CLOSED},
        {2, {1, 4, 1}, 1./3, NewtonCotes::CLOSED},
        {3, {1, 3, 3, 1}, 3./8, NewtonCotes::CLOSED},
        {4, {7, 32, 12, 32, 7}, 2./45, NewtonCotes::CLOSED}
    },
    {
        {3, {1, 1}, 3. / 2, NewtonCotes::OPENED},
        {4, {2, -1, 2}, 4./3, NewtonCotes::OPENED},
        {5, {11, 1, 1, 11}, 5./24, NewtonCotes::OPENED},
        {6, {22, -28, 52, -28, 22}, 3./20, NewtonCotes::OPENED}
    }
};
NewtonCotes::NewtonCotes(Type type, Degree degree) : tableValues(NewtonCotesTable[type][degree]) {}

double NewtonCotes::intervalIntegrator(const Function &f, double upper, double lower) {
    double interval = (upper - lower) / tableValues.intervalsNum;
    if(tableValues.type == OPENED)
        lower += interval;

    double sum = 0;
    for(uint32_t i = 0; i < tableValues.weights.size(); i++)
        sum += tableValues.weights[i] * f.eval(lower + i * interval);

    return tableValues.alfa * interval * sum;
}

const std::vector<GaussLegendre::TableValues> GaussLegendreTable =
{
    {{1, 1}, {-1. / sqrt(3), 1. / sqrt(3)}},
    {{.8888888888888888, .5555555555555556, .5555555555555556}, {0, -0.7745966692414834, 0.7745966692414834}},
    {{0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538}, {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526}},
    {{0.01761400713915, 0.0406014298004, 0.06267204833411, 0.083276741576705, 0.10193011981724, 0.1181945319615, 0.13168863844918, 0.14209610931838, 0.1491729864726, 0.15275338713073, 0.1527533871307, 0.1491729864726, 0.1420961093184, 0.13168863844918, 0.11819453196152, 0.1019301198172, 0.0832767415767, 0.06267204833411, 0.04060142980039, 0.0176140071392}, {-0.9931285991851, -0.96397192727791, -0.91223442825133, -0.83911697182222, -0.74633190646015, -0.63605368072652, -0.51086700195083, -0.37370608871542, -0.22778585114165, -0.0765265211335, 0.0765265211335, 0.22778585114165, 0.37370608871542, 0.51086700195083, 0.63605368072652, 0.74633190646015, 0.8391169718222, 0.91223442825133, 0.96397192727791, 0.9931285991851}}
};
GaussLegendre::GaussLegendre(GaussLegendre::Degree degree) : tableValues(GaussLegendreTable[degree]) {}

double GaussLegendre::intervalIntegrator(const Function &f, double upper, double lower) {
    double sum = 0;
    for(uint32_t i = 0; i < tableValues.weights.size(); i++)
        sum += tableValues.weights[i] * f.eval(changeOfVariable(tableValues.alfa[i], upper, lower));

    return (upper - lower) * .5 * sum;
}

double GaussLegendre::changeOfVariable(double alfa, double upper, double lower) {
    return .5 * ((upper + lower) + (upper - lower) * alfa);
}

std::string replaceAll(std::string str, const std::string &from, const std::string &to) {
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
    return str;
}

double doubleExponential(const Function& f, double upper, double lower, double c) {
    GaussLegendre gl20(GaussLegendre::TWENTY);

    std::string halfPi = std::to_string(M_PI / 2.);
    std::string dx_ds = std::to_string((upper-lower)*.5);
    std::string x_s = std::to_string((lower+upper)*.5) + " + " + dx_ds + " * tanh(" + halfPi + " * sinh(x) )";
    dx_ds += " * " + halfPi + " * ( (cosh(x)) / (cosh(" + halfPi + " * sinh(x)))^2 )";

    std::string newFunction = replaceAll(f.expression, "x", "(" + x_s + ")") + " * (" + dx_ds + ")";

    return gl20.integrate(Function(newFunction), c, -c, .0000001).finalResult();
}

double simpleExponential(const Function& f, double upper, double lower, double c) {
    GaussLegendre gl20(GaussLegendre::TWENTY);

    std::string halfPi = std::to_string(M_PI / 2.);
    std::string dx_ds = std::to_string((upper-lower)*.5);
    std::string x_s = std::to_string((lower+upper)*.5) + " + " + dx_ds + " * tanh(x)";
    dx_ds +=  " / (cosh(x)^2)";

    std::string newFunction = replaceAll(f.expression, "x", "(" + x_s + ")") + " * (" + dx_ds + ")";

    return gl20.integrate(Function(newFunction), c, -c, .0000001).finalResult();
}