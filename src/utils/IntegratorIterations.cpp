#include "IntegratorIterations.hpp"

IntegratorIterations::IntegratorIterations() = default;

void IntegratorIterations::add(uint32_t intervalsNum, double integralValue, double error) {
    results.push_back({intervalsNum, integralValue, error});
}

void IntegratorIterations::printResults() {
    printf("N de intervalos               Integral                    Erro atual\n");
    for(auto [intervalsNum, IntegralValue, Error] : results) {
        printf("%s%u                     %s%.8f%s%.8f\n",
               std::string(8 - std::to_string(intervalsNum).length(), ' ').c_str(),
               intervalsNum,
               IntegralValue < 0 ? "" : " ",
               IntegralValue,
               std::string(26 - std::to_string(IntegralValue).length(), ' ').c_str(),
               Error);
    }
}

double IntegratorIterations::finalResult() {
    return results.back().integralValue;
}
