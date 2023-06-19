#pragma once

#include <cstdint>
#include <string>
#include <vector>

class IntegratorIterations {
    struct Iteration {
        uint32_t intervalsNum;
        double integralValue;
        double error;
    };

public:
    IntegratorIterations();

    void add(uint32_t intervalsNum, double integralValue, double error);
    void printResults();
    double finalResult();

private:
    std::vector<Iteration> results;

};
