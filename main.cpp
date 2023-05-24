#include "src/NumMethods.hpp"

int main() {
    Matrix M1(3, 3);
    M1 << 5, 2, 1,
          2, 3, 1,
          1, 1, 2;
    LOG("\nPara a Matrix M1: ");
    LOG(M1);

    powerMethod::print(powerMethod::regular(M1, .000001));
    powerMethod::print(powerMethod::inverse(M1, .000001));
    powerMethod::print(powerMethod::shifted(M1, 3, .000001));

    Matrix M2(3, 3);
    M2 << -14,  1,  -2,
            1, -1,   1,
           -2,  1, -11;
    LOG("\nPara a Matrix M2: ");
    LOG(M2);

    powerMethod::print(powerMethod::regular(M2, .000001));
    powerMethod::print(powerMethod::inverse(M2, .000001));
    powerMethod::print(powerMethod::shifted(M2, -10, .000001));


    Matrix M3(5, 5);
    M3 << 40,  8,  4,  2, 1,
          8, 30, 12,  6, 2,
          4, 12, 20,  1, 2,
          2,  6,  1, 25, 4,
          1,  2,  2,  4, 5;
    LOG("\nPara a Matrix M3: ");
    LOG(M3);

    powerMethod::print(powerMethod::regular(M3, .000001));
    powerMethod::print(powerMethod::inverse(M3, .000001));
    powerMethod::print(powerMethod::shifted(M3, 11, .000001));
    powerMethod::print(powerMethod::shifted(M3, 23, .000001));
    powerMethod::print(powerMethod::shifted(M3, 32, .000001));

    HH::result HH_M3 = HH::makeHHMatrix(M3);
    eigenResult regular = powerMethod::regular(HH_M3.triDiagM, .000001);
    regular.vector = HH_M3.accumH * regular.vector;
    powerMethod::print(regular);
    eigenResult inverse = powerMethod::inverse(HH_M3.triDiagM, .000001);
    powerMethod::print(inverse);


    eigenResult s1 = powerMethod::shifted(HH_M3.triDiagM, 11, .000001);
    powerMethod::print(s1);
    eigenResult s2 = powerMethod::shifted(HH_M3.triDiagM, 11, .000001);
    powerMethod::print(s2);
    eigenResult s3 = powerMethod::shifted(HH_M3.triDiagM, 11, .000001);
    powerMethod::print(s3);
}
