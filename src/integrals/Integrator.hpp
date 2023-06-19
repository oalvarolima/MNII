#pragma once

#include "../utils/Fn.hpp"
#include "../utils/IntegratorIterations.hpp"

class Integrator {
public:
    IntegratorIterations integrate(const Fn& f, double upper, double lower, double relativeErrorTolerance);

protected:
    double integrateIntervals(const Fn &f, uint32_t intervalsNum, double upper, double lower);
    virtual double intervalIntegrator(const Fn &f, double lower, double upper) = 0;

protected:
    const uint32_t maxIntervalsNum = 1024;
};

class NewtonCotes : public Integrator {
public:
    enum Type {
        CLOSED = 0,
        OPENED = 1
    };

    enum Degree {
        TWO = 0,
        THREE = 1,
        FOUR = 2,
        FIVE = 3
    };
    NewtonCotes(Type type, Degree degree);

    struct TableValues {
        double intervalsNum;
        std::vector<double> weights;
        double alfa;
        Type type;
    };

private:
    double intervalIntegrator(const Fn &f, double lower, double upper) override;

private:
    const TableValues tableValues;
};

class GaussLegendre : public Integrator {
public:
    enum Degree {
        TWO = 0,
        THREE = 1,
        FOUR = 2,
        TWENTY = 3
    };

    GaussLegendre(Degree degree);

    struct TableValues {
        std::vector<double> weights;
        std::vector<double> alfa;
    };

private:
    double intervalIntegrator(const Fn &f, double lower, double upper) override;
    static double changeOfVariable(double alfa, double upper, double lower);

private:
    const TableValues tableValues;
};