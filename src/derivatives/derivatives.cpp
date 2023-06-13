#include "derivatives.hpp"

double secondDeriv_ctr_E4(const Fn &f, double x, double interval) {
    double fx_m_2i = f.eval(x - 2*interval);
    double fx_m_i = f.eval(x - interval);
    double fx = f.eval(x);
    double fx_p_i = f.eval(x + interval);
    double fx_p_2i = f.eval(x + 2*interval);
    return (-fx_m_2i + 16*fx_p_i - 30*fx + 16*fx_m_i - fx_p_2i) *
           ( 1./(12*interval*interval) );
}
