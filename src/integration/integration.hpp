#pragma once

#include "../utils/Fn.hpp"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>

#define LOG(x) std::cout << x << std::endl;

struct itrInfo
{
    uint32_t itrNum;
    double Ivalue;
    double EPS;
};

std::vector<itrInfo> integrate(const Fn &f, double upper, double lower, double EPS, double (*integrator)(const Fn &, double, double, uint32_t), uint32_t max_itr, uint32_t offSet);
double trapezoidalRule(const Fn &f,  double upper, double lower, uint32_t intervalsNum);
double simpson13(const Fn &f,  double upper, double lower, uint32_t intervalsNum);
double simpson38(const Fn &f, double upper, double lower, uint32_t intervalsNum);

void printResults(const std::vector<itrInfo> &results, const std::string &methodName);
