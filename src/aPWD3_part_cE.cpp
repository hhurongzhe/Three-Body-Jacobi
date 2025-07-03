#include <iostream>
#include <complex>
#include <tuple>
#include <vector>
#include <cmath>

#include "constants.hpp"
#include "./numlib/spherical_harmonics.hpp"
#include "aPWD3_part_cE.hpp"

namespace aPWD3_cE
{
// to treat "integer \times complex".
std::complex<double> operator*(const int &c1, const std::complex<double> &c2) { return (c1 * 1.0) * c2; }
using constants::fpi;
using constants::gA;
using constants::Mpi;
using constants::PI;
using spherical_harmonics::Ybra;
using spherical_harmonics::Yket;

constexpr double sqrt3 = 1.732050807568877;

double Sin(double x) { return sin(x); }
double Cos(double x) { return cos(x); }
double Power(double x, int n) { return std::pow(x, n); }
double Sqrt(double x) { return std::sqrt(x); }
std::complex<double> Complex(double x, double y)
{
    std::complex<double> value{x, y};
    return value;
}

// spin projection of 3N forces.
std::complex<double> Gt(int idxbra, int idxket, double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp, util::WignerSymbols &wigner)
{
    const double Fct = 1.0 / (fpi * fpi * fpi * fpi);

    std::complex<double> mtx;
    if (idxbra == -1 && idxket == -1)
    {
        return 0.0;
    }
    else if (idxbra == 1 && idxket == 1)
    {
        mtx = Fct * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 1 && idxket == 10)
    {
        mtx = Fct * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 2 && idxket == 2)
    {
        mtx = -3 * Fct * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 2 && idxket == 12)
    {
        mtx = -3 * Fct * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 3 && idxket == 3)
    {
        mtx = (-3 * Fct *
               (Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) + Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) + Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner))) /
              5.;
        return mtx;
    }
    else if (idxbra == 3 && idxket == 15)
    {
        mtx = (-3 * Fct *
               (Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) + Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) + Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner))) /
              5.;
        return mtx;
    }
    else if (idxbra == 3 && idxket == 16)
    {
        mtx = (-3 * Fct *
               (Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) + Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) + Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner))) /
              5.;
        return mtx;
    }
    else if (idxbra == 4 && idxket == 4)
    {
        mtx = -3 * Fct * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 1, 1, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 5 && idxket == 5)
    {
        mtx = -(Fct * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 6 && idxket == 6)
    {
        mtx = Fct * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 1, 1, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 7 && idxket == 7)
    {
        mtx = (Fct * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner))) / 3.;
        return mtx;
    }
    else if (idxbra == 8 && idxket == 8)
    {
        mtx = (Fct * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner))) / 3.;
        return mtx;
    }
    else if (idxbra == 9 && idxket == 9)
    {
        mtx = (Fct * (Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) + Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) + Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) +
                      Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) + Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner))) /
              5.;
        return mtx;
    }
    else if (idxbra == 10 && idxket == 1)
    {
        mtx = Fct * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 10 && idxket == 10)
    {
        mtx = Fct * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 11 && idxket == 11)
    {
        mtx = (Fct * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner))) / 3.;
        return mtx;
    }
    else if (idxbra == 12 && idxket == 2)
    {
        mtx = -3 * Fct * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 12 && idxket == 12)
    {
        mtx = -3 * Fct * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 13 && idxket == 13)
    {
        mtx = -(Fct * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 14)
    {
        mtx = -(Fct * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 3)
    {
        mtx = (-3 * Fct *
               (Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) + Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) + Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner))) /
              5.;
        return mtx;
    }
    else if (idxbra == 15 && idxket == 15)
    {
        mtx = (-3 * Fct *
               (Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) + Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) + Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner))) /
              5.;
        return mtx;
    }
    else if (idxbra == 15 && idxket == 16)
    {
        mtx = (-3 * Fct *
               (Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) + Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) + Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner))) /
              5.;
        return mtx;
    }
    else if (idxbra == 16 && idxket == 3)
    {
        mtx = (-3 * Fct *
               (Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) + Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) + Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner))) /
              5.;
        return mtx;
    }
    else if (idxbra == 16 && idxket == 15)
    {
        mtx = (-3 * Fct *
               (Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) + Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) + Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner))) /
              5.;
        return mtx;
    }
    else if (idxbra == 16 && idxket == 16)
    {
        mtx = (-3 * Fct *
               (Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) + Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) + Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner))) /
              5.;
        return mtx;
    }

    else
    {
        // std::cerr << "unkown channel index: (" << idxbra << "," << idxket << ") !" << std::endl;
        return 0;
    }
}
} // end namespace aPWD3_cE
