#pragma once

#include "mtxUtils.hpp"

namespace HouseHolder {
    struct result {
        Matrix triDiagM;
        Matrix accumHHs;
    };

    result makeTridiagMatrix(const Matrix& A);
}
