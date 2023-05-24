#pragma once

#include <cmath>
#include <iostream>

#include "mtxUtils.hpp"

namespace powerMethod {
    eigenResult regular(const Matrix& A, double tolerance); 
    eigenResult inverse(const Matrix& A, double tolerance); 
    eigenResult shifted(const Matrix& A, double shiftment, double tolerance); 

    void print(eigenResult a);
}
