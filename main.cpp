#include "src/NumMethods.hpp"

int main() {
    const double EPS = .000001;

    Matrix M1(5, 5);
    M1 << 40, 8, 4, 2, 1,
            8, 30, 12, 6, 2,
            4, 12, 20, 1, 2,
            2, 6, 1, 25, 4,
            1, 2, 2, 4, 5;

    LOG("\nPara a matrix M1");
    powerMethod::print(powerMethod::regular(M1, EPS));
    powerMethod::print(powerMethod::shifted(M1, 33, EPS));
    powerMethod::print(powerMethod::shifted(M1, 20, EPS));
    powerMethod::print(powerMethod::shifted(M1, 11, EPS));
    powerMethod::print(powerMethod::inverse(M1, EPS));


    HouseHolder::result result = HouseHolder::makeTridiagMatrix(M1);
    LOG("\n\nMatrix triadiagonal");
    LOG(result.triDiagM);

    LOG("\n\nH1*H2*H3");
    LOG(result.accumHHs);

    LOG("\n\n\n Usando a Matrix tridiagonal");

    eigenResult e1 = powerMethod::regular(result.triDiagM, EPS);
    e1.vector = result.accumHHs * e1.vector;
    powerMethod::print(e1);

    eigenResult e2 = powerMethod::shifted(result.triDiagM, 33, EPS);
    e2.vector = result.accumHHs * e2.vector;
    powerMethod::print(e2);

    eigenResult e3 = powerMethod::shifted(result.triDiagM, 20, EPS);
    e3.vector = result.accumHHs * e3.vector;
    powerMethod::print(e3);

    eigenResult e4 = powerMethod::shifted(result.triDiagM, 11, EPS);
    e4.vector = result.accumHHs * e4.vector;
    powerMethod::print(e4);

    eigenResult e5 = powerMethod::inverse(result.triDiagM, EPS);
    e5.vector = result.accumHHs * e5.vector;
    powerMethod::print(e5);

    QR::methodResult qrresult = QR::method(M1, EPS);
    LOG("\n\nauto-valores por QR");
    LOG(qrresult.eigenValues);

    LOG("\n\nQ1*Q2*Q3");
    LOG(qrresult.accumQs);
}