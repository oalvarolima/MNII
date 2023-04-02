#include "src/NumMethods.hpp"

void tarefa1()
{
    // A formula da derivada central com erro O(h^4) pela
    // diferença das equações de taylor é igual a fórmula
    // pela interpolação de Newton.
    Fn f("sqrt(e^(3*x) + 4*(x^2))");
    
    uint32_t max_itr = 50;
    double EPS = 69;
    double deriv = 69.69, old_deriv = 0;
    double interval = .5;
    printf("\n   intervalo               f(x)                      f''(x)                   EPS\n\n"); 
    while( EPS > .000001 && max_itr-- )
    {
        deriv = secondDeriv_ctr_E4(f, 2, interval);
        EPS = fabs((deriv - old_deriv) / deriv);
        printf("    %.5f               %.5f                  %.5f                 %.5f\n",
                interval,
                20.48,
                deriv,
                EPS);
        interval /= 2;
        old_deriv = deriv;
    }
}

void tarefa1_2()
{
    Matrix imgGrayScaleColorMatrix = grayScaleImgMatrix("ednaldoPereira.jpeg");
    Matrix smoothdImg = applyGaussianFilter(imgGrayScaleColorMatrix);
    writeImgOfMatrix(edgeDetector(smoothdImg), "ednaldoPereiraSobel.jpg");
    writeImgOfMatrix(laplacianEdgeDetector(imgGrayScaleColorMatrix), "ednaldoPereiraLP.jpg");
}

int main()
{
   tarefa1(); 
   tarefa1_2();
   //runIcalculator();
}
