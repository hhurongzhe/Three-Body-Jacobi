#pragma once
#ifndef aPWD3_c1_HPP
#define aPWD3_c1_HPP

#include <complex>
#include "./numlib/WignerSymbol.hpp"

namespace aPWD3_c1
{
// spin projection of 3N forces.
std::complex<double> Gt(int idxbra, int idxket, double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp, util::WignerSymbols &wigner);
} // end namespace aPWD3_c1

#endif // aPWD3_c1_HPP