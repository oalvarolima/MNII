#include "derivatives.hpp"

const std::vector<std::vector<std::vector<Derivator::Table>>> tableValues = {
        {//forward
                {//first
                        {1      , {-1., 1.}              , 1},//E1
                        {.5     , {-3., 4., -1}          , 1},//E2
                        {1. / 6 , {-11., 18, -9,  2}     , 1},//E3
                        {1. / 12, {-25., 48, -36, 16, -3}, 1},//E4
                },
                {//second
                        {1,       {1, 2,  1}                     , 2},
                        {1,       {2,  -5, 4,   -1}              , 2},
                        {1. / 12, {35, -104, 114, -56, 11}       , 2},
                        {1. / 12, {45, -154, 214, -156, 61,  -10}, 2},
                },
                {//third
                        {1,  {-1, 3, -3, 1},     3},
                        {.5,     {-5, 18, -24, 14, -3},    3},
                        {1. / 4,   {-17, 71, -118, 98,  -41,  7},           3},
                        {1. / 8, {-49, 232, -461, 496, -307, 104, -15}, 3},
                }
        },
        {//backward
                {//first
                        {1      , {1., -1.}             , 1},//E1
                        {.5     , {3., -4., 1}          , 1},//E2
                        {1. / 6 , {11., -18, 9,  -2}    , 1},//E3
                        {1. / 12, {25., -48, 36, -16, 3}, 1},//E4
                },
                {//second
                        {1,       {1, -2,  1}                     , 2},
                        {1,       {2,  -5, 4,   -1}              , 2},
                        {1. / 12, {35, -104, 114, -56, 11}       , 2},
                        {1. / 12, {45, -154, 214, -156, 61, -10}, 2},
                },
                {//third
                        {1,  {1, -3, 3, -1},     3},
                        {.5,     {5, -18, -24, -14, 3},    3},
                        {1. / 4,   {17, -71, 118, -98,  41,  -7},           3},
                        {1. / 8, {49, -232, 461, -496, 307, -104, 15}, 3},
                }
        },
        {//central
                {
                        {.5, {-1,  0, 1}, 1},//E2
                        {1. / 12, {-1,  -8, 0, 8, -1}, 1},//E4
                        {1. / 60, {-1,   9,  -45, 0, 45, -9, 1}, 1},//E6
                        {1. / 840, {3,   -32, 168, -672, 0, 672, -168, 32, -3}, 1},//E8
                },
                {
                        {1., {1, -2, 1}, 2},
                        {1. / 12, {-1, 16, -30, 16, -1}, 2},
                        {1. / 180, {2,  -27,  270, -490, 270, -27, 2}, 2},
                        {1. / 5040, {-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9}, 2},
                },
                {
                        {.5, {-1, 2, 0,  -2, 1}, 3},
                        {1. / 8, {1,  -8, 13,  0,  8, -1}, 3},
                        {1. / 240, {-7,  72, -338, 488, -488, 338, -72, 7}, 3},
                        {1. / 30240, {205, -2522, 14607, -52428, 70098, 0, -70098, 52428, -14607, 2522, -205}, 3},
                }
        },
};

Derivator::Derivator(Derivator::Type type, Derivator::Order order, Derivator::Degree degree) : values(tableValues[type][order][degree]), type(type), order(order) {}

double Derivator::derivate(const Function &fn, double x, double interval) const {
    int direction = type == BACKWARD? -1 : 1;
    x = type == CENTRAL ? x - interval * (values.weights.size() >> 1) : x;

    double result = 0;
    for (int i = 0; i < values.weights.size(); ++i)
        result += values.weights[i] * fn.eval(x + i * direction * interval);

    return result * values.multiplier / pow(interval, values.intervalExponent);
}