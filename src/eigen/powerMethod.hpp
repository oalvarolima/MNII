#pragma once

#include <cmath>
#include <iostream>
#include "../utils/Fn.hpp"

#include "mtxUtils.hpp"

namespace powerMethod {
    eigenResult regular(const Matrix &A, double tolerance);
    eigenResult inverse(const Matrix &A, double tolerance);
    eigenResult shifted(const Matrix &A, double shiftment, double tolerance);

    struct allEigenResults {
        Matrix eigenVectors;
        Matrix eigenValues;
    };
    allEigenResults getAllEigenResults(const Matrix &A, double tolerance);

    void print(const eigenResult &result);
}
