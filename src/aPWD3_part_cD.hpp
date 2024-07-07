#pragma once
#ifndef aPWD3_cD_HPP
#define aPWD3_cD_HPP

#include <complex>
#include "./numlib/WignerSymbol.hpp"

namespace aPWD3_cD
{
    // spin projection of 3N forces.
    std::complex<double> Gt(int idxbra, int idxket, double ppmag, double qpmag, double pmag, double qmag,
                            double thetaq, double thetapp, double phipp, double thetaqp, double phiqp,
                            util::WignerSymbols &wigner);
} // end namespace aPWD3_cD

#endif // aPWD3_cD_HPP