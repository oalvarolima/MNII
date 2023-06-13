#pragma once

#include "mtxUtils.hpp"

namespace HouseHolder {
    struct result {
        Matrix triDiagM;
        Matrix accumHHs;
    };

    result makeTridiagMatrix(const Matrix &A);
    Matrix makeHHcolI(const Matrix &A, uint32_t i);
}
