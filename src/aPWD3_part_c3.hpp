#pragma once
#ifndef aPWD3_c3_HPP
#define aPWD3_c3_HPP

#include <complex>
#include "./numlib/WignerSymbol.hpp"

namespace aPWD3_c3
{
    // spin projection of 3N forces.
    std::complex<double> Gt(int idxbra, int idxket, double ppmag, double qpmag, double pmag, double qmag,
                            double thetaq, double thetapp, double phipp, double thetaqp, double phiqp,
                            util::WignerSymbols &wigner);
} // end namespace aPWD3_c3

#endif // aPWD3_c3_HPP