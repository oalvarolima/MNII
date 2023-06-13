#pragma once

#include "powerMethod.hpp"

namespace SVD {
    struct result {
        Matrix U, S, V;
    };

    result decompSVD(const Matrix &A);
}