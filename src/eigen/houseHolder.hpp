#pragma once

#include "mtxUtils.hpp"
#include "../utils/Fn.hpp"

namespace HH {
    struct result {
        Matrix triDiagM;
        Matrix accumH;
    };

    result makeHHMatrix(const Matrix& A); 
}
