#include <iostream>
#include <complex>
#include <tuple>
#include <vector>
#include <cmath>

#include "constants.hpp"
#include "./numlib/spherical_harmonics.hpp"
#include "aPWD3_part_c1.hpp"

namespace aPWD3_c1
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
double get_q2x(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = ppmag * Cos(phipp) * Sin(thetapp) + (qmag * Sin(thetaq)) / 2. - (qpmag * Cos(phiqp) * Sin(thetaqp)) / 2.;
    return value;
}
double get_q2y(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = ppmag * Sin(phipp) * Sin(thetapp) - (qpmag * Sin(phiqp) * Sin(thetaqp)) / 2.;
    return value;
}
double get_q2z(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = -pmag + ppmag * Cos(thetapp) + (qmag * Cos(thetaq)) / 2. - (qpmag * Cos(thetaqp)) / 2.;
    return value;
}
double get_q3x(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = -(ppmag * Cos(phipp) * Sin(thetapp)) + (qmag * Sin(thetaq)) / 2. - (qpmag * Cos(phiqp) * Sin(thetaqp)) / 2.;
    return value;
}
double get_q3y(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = -(ppmag * Sin(phipp) * Sin(thetapp)) - (qpmag * Sin(phiqp) * Sin(thetaqp)) / 2.;
    return value;
}
double get_q3z(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = pmag - ppmag * Cos(thetapp) + (qmag * Cos(thetaq)) / 2. - (qpmag * Cos(thetaqp)) / 2.;
    return value;
}

// spin projection of 3N forces.
std::complex<double> Gt(int idxbra, int idxket, double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp, util::WignerSymbols &wigner)
{
    double q2x = get_q2x(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q2y = get_q2y(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q2z = get_q2z(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q3x = get_q3x(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q3y = get_q3y(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q3z = get_q3z(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q2q2 = q2x * q2x + q2y * q2y + q2z * q2z;
    double q3q3 = q3x * q3x + q3y * q3y + q3z * q3z;
    double Ftpe1 = Power(gA / (2.0 * fpi), 2) / (q2q2 + Mpi * Mpi) / (q3q3 + Mpi * Mpi) * (-4.0 * Mpi * Mpi) / (fpi * fpi);

    std::complex<double> mtx;
    if (idxbra == -1 && idxket == -1)
    {
        return 0.0;
    }
    else if (idxbra == 1 && idxket == 1)
    {
        mtx = -(Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 1 && idxket == 7)
    {
        mtx = (Ftpe1 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 1 && idxket == 8)
    {
        mtx = (Ftpe1 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              3.;
        return mtx;
    }
    else if (idxbra == 1 && idxket == 10)
    {
        mtx = -(Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 2 && idxket == 2)
    {
        mtx = -(Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 2 && idxket == 3)
    {
        mtx = -((Ftpe1 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                 (Sqrt(3) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 0, 2, thetaq, wigner) + Sqrt(3) * (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) - Sqrt(2) * q2x * q3x * Yket(2, 0, 0, 2, thetaq, wigner) -
                  Sqrt(2) * q2y * q3y * Yket(2, 0, 0, 2, thetaq, wigner) + 2 * Sqrt(2) * q2z * q3z * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 1, 0, 2, thetaq, wigner) -
                  Sqrt(3) * q2x * q3z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 1, 0, 2, thetaq, wigner) + Sqrt(3) * q2x * q3x * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3x * Yket(2, 2, 0, 2, thetaq, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Yket(2, 2, 0, 2, thetaq, wigner) - Sqrt(3) * q2y * q3y * Yket(2, 2, 0, 2, thetaq, wigner))) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 2 && idxket == 5)
    {
        mtx = (Ftpe1 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              Sqrt(2);
        return mtx;
    }
    else if (idxbra == 2 && idxket == 12)
    {
        mtx = -(Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 2 && idxket == 15)
    {
        mtx = -((Ftpe1 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                 (Sqrt(3) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 2, 0, thetaq, wigner) + Sqrt(3) * (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) - Sqrt(2) * q2x * q3x * Yket(2, 0, 2, 0, thetaq, wigner) -
                  Sqrt(2) * q2y * q3y * Yket(2, 0, 2, 0, thetaq, wigner) + 2 * Sqrt(2) * q2z * q3z * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 1, 2, 0, thetaq, wigner) -
                  Sqrt(3) * q2x * q3z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 1, 2, 0, thetaq, wigner) + Sqrt(3) * q2x * q3x * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3x * Yket(2, 2, 2, 0, thetaq, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Yket(2, 2, 2, 0, thetaq, wigner) - Sqrt(3) * q2y * q3y * Yket(2, 2, 2, 0, thetaq, wigner))) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 2 && idxket == 16)
    {
        mtx = -((Ftpe1 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                 (Sqrt(3) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 2, 2, thetaq, wigner) + Sqrt(3) * (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) - Sqrt(2) * q2x * q3x * Yket(2, 0, 2, 2, thetaq, wigner) -
                  Sqrt(2) * q2y * q3y * Yket(2, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(2) * q2z * q3z * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 1, 2, 2, thetaq, wigner) -
                  Sqrt(3) * q2x * q3z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 1, 2, 2, thetaq, wigner) + Sqrt(3) * q2x * q3x * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3x * Yket(2, 2, 2, 2, thetaq, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Yket(2, 2, 2, 2, thetaq, wigner) - Sqrt(3) * q2y * q3y * Yket(2, 2, 2, 2, thetaq, wigner))) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 2)
    {
        mtx = -((Ftpe1 *
                 (Sqrt(3) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(2) * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2z * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(3) * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * q2x * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2y * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 0, 0, thetaq, wigner)) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 3)
    {
        mtx =
            -0.05 *
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      2 * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 0, 2, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 0, 2, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) -
                      Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 0, 2, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 0, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 0, 2, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 0, 2, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 0, 2, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 0, 2, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 5)
    {
        mtx = (Ftpe1 * (2 * Sqrt(3) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Sqrt(2) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 4) * q2y * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * q2x * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Sqrt(3) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, 2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + Sqrt(2) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner)) -
                        Sqrt(2) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2y * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2x * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        2 * Sqrt(3) * q2z * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2z * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        2 * Sqrt(3) * q2x * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2y * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              (2. * Sqrt(10));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 12)
    {
        mtx = -((Ftpe1 *
                 (Sqrt(3) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(2) * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2z * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(3) * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * q2x * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2y * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 2, 2, thetaq, wigner)) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 13)
    {
        mtx = (Ftpe1 * (-(Sqrt(3) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner)) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2x * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        2 * q2x * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * q2x * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 2 * q2y * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - 2 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 0, 2, 2, thetaq, wigner)) -
                        Sqrt(3) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        2 * Sqrt(2) * q2z * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2z * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2x * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                     Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              (2. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 14)
    {
        mtx = -0.25 *
              (Ftpe1 * (Sqrt(6) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        2 * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * q2y * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * q2x * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 2 * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        2 * Sqrt(2) * q2x * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * q2y * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(2) * q2x * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(2) * q2y * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        2 * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        Sqrt(6) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        2 * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        4 * q2z * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * q2z * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * q2z * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * q2x * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * q2y * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                     2 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              Sqrt(5);
        return mtx;
    }
    else if (idxbra == 3 && idxket == 15)
    {
        mtx =
            -0.05 *
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      2 * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 2, 0, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 2, 0, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) -
                      Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 0, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 2, 0, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 2, 0, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 2, 0, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 2, 0, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 2, 0, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 16)
    {
        mtx =
            -0.05 *
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      2 * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 2, 2, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 2, 2, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) -
                      Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 2, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 2, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 2, 2, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 2, 2, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 2, 2, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 2, 2, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 4 && idxket == 4)
    {
        mtx = 3 * Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 1, 1, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 4 && idxket == 13)
    {
        mtx = (Ftpe1 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 2, 2, thetaq, wigner) +
                (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner))) /
              Sqrt(2);
        return mtx;
    }
    else if (idxbra == 4 && idxket == 14)
    {
        mtx = Ftpe1 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
              ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 2, 2, thetaq, wigner) +
               (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 5 && idxket == 2)
    {
        mtx = (Ftpe1 *
               ((-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              Sqrt(2);
        return mtx;
    }
    else if (idxbra == 5 && idxket == 3)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(6) * (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * (q2y * q3x - q2x * q3y) * Yket(2, -1, 0, 2, thetaq, wigner) +
                                                                                     (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 0, 0, 2, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(3) * (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(2, 0, 0, 2, thetaq, wigner) +
                                                                                    Sqrt(3) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 1, 0, 2, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 0, 2, thetaq, wigner) +
                                                                                    Complex(0, 1) * Sqrt(6) * ((-(q2y * q3x) + q2x * q3y) * Yket(2, 1, 0, 2, thetaq, wigner) + (Complex(0, 1) * q2z * q3x + q2z * q3y - Complex(0, 1) * q2x * q3z - q2y * q3z) * Yket(2, 2, 0, 2, thetaq, wigner))))) /
              (2. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 5 && idxket == 5)
    {
        mtx = Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) *
              (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 5 && idxket == 12)
    {
        mtx = (Ftpe1 *
               ((-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              Sqrt(2);
        return mtx;
    }
    else if (idxbra == 5 && idxket == 13)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, -1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 2, 2, thetaq, wigner) + (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              Sqrt(6);
        return mtx;
    }
    else if (idxbra == 5 && idxket == 14)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, -1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 2, 2, thetaq, wigner) + (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              (2. * Sqrt(3));
        return mtx;
    }
    else if (idxbra == 5 && idxket == 15)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(6) * (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * (q2y * q3x - q2x * q3y) * Yket(2, -1, 2, 0, thetaq, wigner) +
                                                                                     (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 0, 2, 0, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(3) * (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(2, 0, 2, 0, thetaq, wigner) +
                                                                                    Sqrt(3) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 1, 2, 0, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 2, 0, thetaq, wigner) +
                                                                                    Complex(0, 1) * Sqrt(6) * ((-(q2y * q3x) + q2x * q3y) * Yket(2, 1, 2, 0, thetaq, wigner) + (Complex(0, 1) * q2z * q3x + q2z * q3y - Complex(0, 1) * q2x * q3z - q2y * q3z) * Yket(2, 2, 2, 0, thetaq, wigner))))) /
              (2. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 5 && idxket == 16)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(6) * (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * (q2y * q3x - q2x * q3y) * Yket(2, -1, 2, 2, thetaq, wigner) +
                                                                                     (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 0, 2, 2, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(3) * (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(2, 0, 2, 2, thetaq, wigner) +
                                                                                    Sqrt(3) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 1, 2, 2, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 2, 2, thetaq, wigner) +
                                                                                    Complex(0, 1) * Sqrt(6) * ((-(q2y * q3x) + q2x * q3y) * Yket(2, 1, 2, 2, thetaq, wigner) + (Complex(0, 1) * q2z * q3x + q2z * q3y - Complex(0, 1) * q2x * q3z - q2y * q3z) * Yket(2, 2, 2, 2, thetaq, wigner))))) /
              (2. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 6 && idxket == 6)
    {
        mtx = (Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 1, 1, thetaq, wigner)) / 3.;
        return mtx;
    }
    else if (idxbra == 6 && idxket == 9)
    {
        mtx = (Ftpe1 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
               (Sqrt(3) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 1, 1, thetaq, wigner) + Sqrt(3) * (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 1, 1, thetaq, wigner) - Sqrt(2) * q2x * q3x * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(2) * q2y * q3y * Yket(2, 0, 1, 1, thetaq, wigner) + 2 * Sqrt(2) * q2z * q3z * Yket(2, 0, 1, 1, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 1, 1, 1, thetaq, wigner) - Sqrt(3) * q2x * q3z * Yket(2, 1, 1, 1, thetaq, wigner) +
                Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 1, 1, 1, thetaq, wigner) + Sqrt(3) * q2x * q3x * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3x * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2x * q3y * Yket(2, 2, 1, 1, thetaq, wigner) -
                Sqrt(3) * q2y * q3y * Yket(2, 2, 1, 1, thetaq, wigner))) /
              (3. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 6 && idxket == 11)
    {
        mtx = (Ftpe1 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 2, 2, thetaq, wigner) +
                (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner))) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 7 && idxket == 1)
    {
        mtx = (Ftpe1 *
               ((-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 7 && idxket == 7)
    {
        mtx = (Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) *
               (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              9.;
        return mtx;
    }
    else if (idxbra == 7 && idxket == 8)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + 3 * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) -
                                                                                     3 * Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) - 2 * Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) -
                                                                                    3 * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)) -
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * Sqrt(2) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + 3 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) -
                                                                                    Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              18.;
        return mtx;
    }
    else if (idxbra == 7 && idxket == 9)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (-(Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 1, 1, thetaq, wigner)) + Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, -1, 1, 1, thetaq, wigner) + Sqrt(3) * q2z * q3x * Yket(2, 0, 1, 1, thetaq, wigner) -
                             Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(3) * q2x * q3z * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 0, 1, 1, thetaq, wigner) - Sqrt(2) * q2x * q3x * Yket(2, 1, 1, 1, thetaq, wigner) +
                             Complex(0, 1) * Sqrt(2) * q2y * q3x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2x * q3y * Yket(2, 1, 1, 1, thetaq, wigner) + Sqrt(2) * q2y * q3y * Yket(2, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 1, 1, thetaq, wigner) + (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 1, 1, thetaq, wigner) +
                                                                                    q2z * q3x * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * q2z * q3y * Yket(2, 1, 1, 1, thetaq, wigner) + q2x * q3z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Yket(2, 1, 1, 1, thetaq, wigner) -
                                                                                    2 * q2x * q3x * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 2) * q2x * q3y * Yket(2, 2, 1, 1, thetaq, wigner) + 2 * q2y * q3y * Yket(2, 2, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -1, 1, 1, thetaq, wigner) + Sqrt(3) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 1, 1, thetaq, wigner) -
                                                                                    Sqrt(2) * ((q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, 1, 1, 1, thetaq, wigner) + (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, 2, 1, 1, thetaq, wigner))))) /
              (6. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 7 && idxket == 10)
    {
        mtx = (Ftpe1 *
               ((-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 7 && idxket == 11)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, -1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 2, 2, thetaq, wigner) + (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              (3. * Sqrt(6));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 1)
    {
        mtx = (Ftpe1 *
               (Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 2) * (q2y * q3x - q2x * q3y) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 7)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + 3 * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) -
                                                                                     3 * Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) - 2 * Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) -
                                                                                    3 * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)) -
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * Sqrt(2) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + 3 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) -
                                                                                    Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              18.;
        return mtx;
    }
    else if (idxbra == 8 && idxket == 8)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x * q3x + q2y * q3y + 4 * q2z * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                     6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 1, 1, thetaq, wigner)) -
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + (-8 * q2x * q3x - 8 * q2y * q3y + 4 * q2z * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) -
                                                                                    3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                    2 * (q2x * q3x + q2y * q3y + 4 * q2z * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              36.;
        return mtx;
    }
    else if (idxbra == 8 && idxket == 9)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (-2 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 1, 1, thetaq, wigner) + 2 * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, -1, 1, 1, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 1, 1, thetaq, wigner) -
                             Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 1, 1, thetaq, wigner) - 2 * q2x * q3x * Yket(2, 1, 1, 1, thetaq, wigner) +
                             Complex(0, 2) * q2y * q3x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * q2x * q3y * Yket(2, 1, 1, 1, thetaq, wigner) + 2 * q2y * q3y * Yket(2, 1, 1, 1, thetaq, wigner)) +
                        Sqrt(2) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 1, 1, thetaq, wigner) + (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 1, 1, thetaq, wigner) + q2z * q3x * Yket(2, 1, 1, 1, thetaq, wigner) -
                             Complex(0, 1) * q2z * q3y * Yket(2, 1, 1, 1, thetaq, wigner) + q2x * q3z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Yket(2, 1, 1, 1, thetaq, wigner) - 2 * q2x * q3x * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Yket(2, 2, 1, 1, thetaq, wigner) +
                             Complex(0, 2) * q2x * q3y * Yket(2, 2, 1, 1, thetaq, wigner) + 2 * q2y * q3y * Yket(2, 2, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -1, 1, 1, thetaq, wigner) + Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 1, 1, thetaq, wigner) -
                                                                                    2 * ((q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, 1, 1, 1, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 2, 1, 1, thetaq, wigner))))) /
              (12. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 10)
    {
        mtx = (Ftpe1 *
               (Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 2) * (q2y * q3x - q2x * q3y) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 11)
    {
        mtx = (Ftpe1 *
               (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, -2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(2) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * (q2y * q3x - q2x * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                Sqrt(2) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              (6. * Sqrt(6));
        return mtx;
    }
    else if (idxbra == 9 && idxket == 6)
    {
        mtx = (Ftpe1 *
               (Sqrt(3) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Sqrt(2) * q2x * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2z * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Sqrt(3) * q2z * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2x * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * q2x * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(3) * q2x * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2y * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              (3. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 9 && idxket == 7)
    {
        mtx = (Ftpe1 * (Sqrt(3) * q2z * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Sqrt(3) * q2x * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Sqrt(2) * q2x * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2y * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2x * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Sqrt(2) * q2y * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        q2z * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * q2z * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        q2x * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * q2y * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        2 * q2x * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 2) * q2y * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * q2x * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 2 * q2y * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) - 2 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                        Sqrt(3) * q2z * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Sqrt(3) * q2x * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Sqrt(2) * q2x * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        2 * Sqrt(2) * q2z * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(2) * q2z * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(2) * q2x * q3z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                     Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (6. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 9 && idxket == 8)
    {
        mtx = (Ftpe1 * (Sqrt(6) * q2z * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        2 * q2x * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 2) * q2y * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * q2x * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * q2y * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Sqrt(2) * q2z * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        2 * Sqrt(2) * q2x * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * q2y * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(2) * q2x * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 2 * Sqrt(2) * q2y * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        2 * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                        Sqrt(6) * q2z * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        2 * q2x * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 2 * q2y * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        4 * q2z * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 2 * q2z * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * q2z * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 2 * q2x * q3z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * q2y * q3z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                     2 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (12. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 9 && idxket == 9)
    {
        mtx =
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) +
                      2 * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 1, 1, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 1, 1, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 1, 1, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                      Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 1, 1, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 1, 1, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 1, 1, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 1, 1, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 1, 1, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 1, 1, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner))) /
            60.;
        return mtx;
    }
    else if (idxbra == 9 && idxket == 11)
    {
        mtx = (Ftpe1 * (2 * Sqrt(3) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2z * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Sqrt(2) * q2x * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 4) * q2y * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 4) * q2x * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2z * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Sqrt(6) * q2x * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(3) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, -2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        Sqrt(2) * q2z * q3x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Sqrt(2) * q2x * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(3) * q2y * q3x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2x * q3y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        2 * Sqrt(3) * q2z * q3x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2z * q3y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        2 * Sqrt(3) * q2x * q3z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2y * q3z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner))) /
              (6. * Sqrt(10));
        return mtx;
    }
    else if (idxbra == 10 && idxket == 1)
    {
        mtx = -(Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 10 && idxket == 7)
    {
        mtx = (Ftpe1 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 10 && idxket == 8)
    {
        mtx = (Ftpe1 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              3.;
        return mtx;
    }
    else if (idxbra == 10 && idxket == 10)
    {
        mtx = -(Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 11 && idxket == 6)
    {
        mtx = (Ftpe1 *
               ((q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                (q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 11 && idxket == 7)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (3. * Sqrt(6));
        return mtx;
    }
    else if (idxbra == 11 && idxket == 8)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (6. * Sqrt(3));
        return mtx;
    }
    else if (idxbra == 11 && idxket == 9)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * (q2y * q3x - q2x * q3y) * Yket(2, -1, 1, 1, thetaq, wigner) +
                                                                                     (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, 0, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(3) * (q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(2, 0, 1, 1, thetaq, wigner) +
                                                                                    Sqrt(3) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 1, 1, thetaq, wigner) +
                                                                                    Complex(0, 1) * Sqrt(6) * ((q2y * q3x - q2x * q3y) * Yket(2, 1, 1, 1, thetaq, wigner) + (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * Yket(2, 2, 1, 1, thetaq, wigner))))) /
              (6. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 11 && idxket == 11)
    {
        mtx =
            -0.3333333333333333 * (Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) *
                                   (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 12 && idxket == 2)
    {
        mtx = -(Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 12 && idxket == 3)
    {
        mtx = -((Ftpe1 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                 (Sqrt(3) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 0, 2, thetaq, wigner) + Sqrt(3) * (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) - Sqrt(2) * q2x * q3x * Yket(2, 0, 0, 2, thetaq, wigner) -
                  Sqrt(2) * q2y * q3y * Yket(2, 0, 0, 2, thetaq, wigner) + 2 * Sqrt(2) * q2z * q3z * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 1, 0, 2, thetaq, wigner) -
                  Sqrt(3) * q2x * q3z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 1, 0, 2, thetaq, wigner) + Sqrt(3) * q2x * q3x * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3x * Yket(2, 2, 0, 2, thetaq, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Yket(2, 2, 0, 2, thetaq, wigner) - Sqrt(3) * q2y * q3y * Yket(2, 2, 0, 2, thetaq, wigner))) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 12 && idxket == 5)
    {
        mtx = (Ftpe1 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              Sqrt(2);
        return mtx;
    }
    else if (idxbra == 12 && idxket == 12)
    {
        mtx = -(Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 12 && idxket == 15)
    {
        mtx = -((Ftpe1 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                 (Sqrt(3) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 2, 0, thetaq, wigner) + Sqrt(3) * (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) - Sqrt(2) * q2x * q3x * Yket(2, 0, 2, 0, thetaq, wigner) -
                  Sqrt(2) * q2y * q3y * Yket(2, 0, 2, 0, thetaq, wigner) + 2 * Sqrt(2) * q2z * q3z * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 1, 2, 0, thetaq, wigner) -
                  Sqrt(3) * q2x * q3z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 1, 2, 0, thetaq, wigner) + Sqrt(3) * q2x * q3x * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3x * Yket(2, 2, 2, 0, thetaq, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Yket(2, 2, 2, 0, thetaq, wigner) - Sqrt(3) * q2y * q3y * Yket(2, 2, 2, 0, thetaq, wigner))) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 12 && idxket == 16)
    {
        mtx = -((Ftpe1 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                 (Sqrt(3) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 2, 2, thetaq, wigner) + Sqrt(3) * (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) - Sqrt(2) * q2x * q3x * Yket(2, 0, 2, 2, thetaq, wigner) -
                  Sqrt(2) * q2y * q3y * Yket(2, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(2) * q2z * q3z * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 1, 2, 2, thetaq, wigner) -
                  Sqrt(3) * q2x * q3z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 1, 2, 2, thetaq, wigner) + Sqrt(3) * q2x * q3x * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3x * Yket(2, 2, 2, 2, thetaq, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Yket(2, 2, 2, 2, thetaq, wigner) - Sqrt(3) * q2y * q3y * Yket(2, 2, 2, 2, thetaq, wigner))) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 13 && idxket == 3)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 0, 2, thetaq, wigner) - Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 0, 0, 2, thetaq, wigner) +
                             Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(3) * q2x * q3z * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(2) * q2x * q3x * Yket(2, 1, 0, 2, thetaq, wigner) -
                             Complex(0, 1) * Sqrt(2) * q2y * q3x * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2x * q3y * Yket(2, 1, 0, 2, thetaq, wigner) - Sqrt(2) * q2y * q3y * Yket(2, 1, 0, 2, thetaq, wigner)) -
                        Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 0, 2, thetaq, wigner) + (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) +
                                                                                    q2z * q3x * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * q2z * q3y * Yket(2, 1, 0, 2, thetaq, wigner) + q2x * q3z * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Yket(2, 1, 0, 2, thetaq, wigner) -
                                                                                    2 * q2x * q3x * Yket(2, 2, 0, 2, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Yket(2, 2, 0, 2, thetaq, wigner) + Complex(0, 2) * q2x * q3y * Yket(2, 2, 0, 2, thetaq, wigner) + 2 * q2y * q3y * Yket(2, 2, 0, 2, thetaq, wigner)) -
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -1, 0, 2, thetaq, wigner) + Sqrt(3) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 0, 2, thetaq, wigner) -
                                                                                    Sqrt(2) * ((q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, 1, 0, 2, thetaq, wigner) + (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, 2, 0, 2, thetaq, wigner))))) /
              (2. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 13 && idxket == 4)
    {
        mtx = (Ftpe1 *
               ((q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                (q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              Sqrt(2);
        return mtx;
    }
    else if (idxbra == 13 && idxket == 5)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              Sqrt(6);
        return mtx;
    }
    else if (idxbra == 13 && idxket == 13)
    {
        mtx =
            -0.3333333333333333 * (Ftpe1 * (q2x * q3x + q2y * q3y + q2z * q3z) *
                                   (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 13 && idxket == 14)
    {
        mtx = -0.16666666666666666 * (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + 3 * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) -
                                                                                                            3 * Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                                               Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - 2 * Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) -
                                                                                                           3 * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner)) -
                                               Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * Sqrt(2) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(1, -1, 2, 2, thetaq, wigner) +
                                                                                                           3 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) - Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, 1, 2, 2, thetaq, wigner))));
        return mtx;
    }
    else if (idxbra == 13 && idxket == 15)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 0, thetaq, wigner) - Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 0, 2, 0, thetaq, wigner) +
                             Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(3) * q2x * q3z * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(2) * q2x * q3x * Yket(2, 1, 2, 0, thetaq, wigner) -
                             Complex(0, 1) * Sqrt(2) * q2y * q3x * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2x * q3y * Yket(2, 1, 2, 0, thetaq, wigner) - Sqrt(2) * q2y * q3y * Yket(2, 1, 2, 0, thetaq, wigner)) -
                        Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 2, 0, thetaq, wigner) + (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) +
                                                                                    q2z * q3x * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * q2z * q3y * Yket(2, 1, 2, 0, thetaq, wigner) + q2x * q3z * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Yket(2, 1, 2, 0, thetaq, wigner) -
                                                                                    2 * q2x * q3x * Yket(2, 2, 2, 0, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Yket(2, 2, 2, 0, thetaq, wigner) + Complex(0, 2) * q2x * q3y * Yket(2, 2, 2, 0, thetaq, wigner) + 2 * q2y * q3y * Yket(2, 2, 2, 0, thetaq, wigner)) -
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -1, 2, 0, thetaq, wigner) + Sqrt(3) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 2, 0, thetaq, wigner) -
                                                                                    Sqrt(2) * ((q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, 1, 2, 0, thetaq, wigner) + (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, 2, 2, 0, thetaq, wigner))))) /
              (2. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 13 && idxket == 16)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 2, thetaq, wigner) - Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) - Sqrt(3) * q2z * q3x * Yket(2, 0, 2, 2, thetaq, wigner) +
                             Complex(0, 1) * Sqrt(3) * q2z * q3y * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(3) * q2x * q3z * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(2) * q2x * q3x * Yket(2, 1, 2, 2, thetaq, wigner) -
                             Complex(0, 1) * Sqrt(2) * q2y * q3x * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2x * q3y * Yket(2, 1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3y * Yket(2, 1, 2, 2, thetaq, wigner)) -
                        Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 2, 2, thetaq, wigner) + (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) +
                                                                                    q2z * q3x * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * q2z * q3y * Yket(2, 1, 2, 2, thetaq, wigner) + q2x * q3z * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Yket(2, 1, 2, 2, thetaq, wigner) -
                                                                                    2 * q2x * q3x * Yket(2, 2, 2, 2, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Yket(2, 2, 2, 2, thetaq, wigner) + Complex(0, 2) * q2x * q3y * Yket(2, 2, 2, 2, thetaq, wigner) + 2 * q2y * q3y * Yket(2, 2, 2, 2, thetaq, wigner)) -
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -1, 2, 2, thetaq, wigner) + Sqrt(3) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 2, 2, thetaq, wigner) -
                                                                                    Sqrt(2) * ((q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, 1, 2, 2, thetaq, wigner) + (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, 2, 2, 2, thetaq, wigner))))) /
              (2. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 3)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 0, 2, thetaq, wigner) - 2 * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) - Sqrt(6) * q2z * q3x * Yket(2, 0, 0, 2, thetaq, wigner) +
                             Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(6) * q2x * q3z * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 0, 2, thetaq, wigner) + 2 * q2x * q3x * Yket(2, 1, 0, 2, thetaq, wigner) -
                             Complex(0, 2) * q2y * q3x * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 2) * q2x * q3y * Yket(2, 1, 0, 2, thetaq, wigner) - 2 * q2y * q3y * Yket(2, 1, 0, 2, thetaq, wigner)) -
                        Sqrt(2) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 0, 2, thetaq, wigner) + (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) + q2z * q3x * Yket(2, 1, 0, 2, thetaq, wigner) -
                             Complex(0, 1) * q2z * q3y * Yket(2, 1, 0, 2, thetaq, wigner) + q2x * q3z * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Yket(2, 1, 0, 2, thetaq, wigner) - 2 * q2x * q3x * Yket(2, 2, 0, 2, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Yket(2, 2, 0, 2, thetaq, wigner) +
                             Complex(0, 2) * q2x * q3y * Yket(2, 2, 0, 2, thetaq, wigner) + 2 * q2y * q3y * Yket(2, 2, 0, 2, thetaq, wigner)) -
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -1, 0, 2, thetaq, wigner) + Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 0, 2, thetaq, wigner) -
                                                                                    2 * ((q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, 1, 0, 2, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 2, 0, 2, thetaq, wigner))))) /
              (4. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 4)
    {
        mtx = -((Ftpe1 *
                 (Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 2) * (q2y * q3x - q2x * q3y) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) - (q2x + Complex(0, 1) * q2y) * q3z) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 1, 1, thetaq, wigner)) /
                Sqrt(2));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 5)
    {
        mtx = (Ftpe1 *
               (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, 2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * (-(q2y * q3x) + q2x * q3y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                Sqrt(2) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((-(q2z * (q3x + Complex(0, 1) * q3y)) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 1, 1, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (2. * Sqrt(6));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 13)
    {
        mtx = -0.16666666666666666 * (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + 3 * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) -
                                                                                                            3 * Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                                               Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - 2 * Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) -
                                                                                                           3 * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner)) -
                                               Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * Sqrt(2) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(1, -1, 2, 2, thetaq, wigner) +
                                                                                                           3 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) - Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, 1, 2, 2, thetaq, wigner))));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 14)
    {
        mtx = -0.08333333333333333 * (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x * q3x + q2y * q3y + 4 * q2z * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - 3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                                            6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)) -
                                               Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + (-8 * q2x * q3x - 8 * q2y * q3y + 4 * q2z * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) -
                                                                                                           3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                                               Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(1, -1, 2, 2, thetaq, wigner) +
                                                                                                           3 * Sqrt(2) * (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) + 2 * (q2x * q3x + q2y * q3y + 4 * q2z * q3z) * Yket(1, 1, 2, 2, thetaq, wigner))));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 15)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 0, thetaq, wigner) - 2 * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) - Sqrt(6) * q2z * q3x * Yket(2, 0, 2, 0, thetaq, wigner) +
                             Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(6) * q2x * q3z * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 2, 0, thetaq, wigner) + 2 * q2x * q3x * Yket(2, 1, 2, 0, thetaq, wigner) -
                             Complex(0, 2) * q2y * q3x * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 2) * q2x * q3y * Yket(2, 1, 2, 0, thetaq, wigner) - 2 * q2y * q3y * Yket(2, 1, 2, 0, thetaq, wigner)) -
                        Sqrt(2) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 2, 0, thetaq, wigner) + (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) + q2z * q3x * Yket(2, 1, 2, 0, thetaq, wigner) -
                             Complex(0, 1) * q2z * q3y * Yket(2, 1, 2, 0, thetaq, wigner) + q2x * q3z * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Yket(2, 1, 2, 0, thetaq, wigner) - 2 * q2x * q3x * Yket(2, 2, 2, 0, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Yket(2, 2, 2, 0, thetaq, wigner) +
                             Complex(0, 2) * q2x * q3y * Yket(2, 2, 2, 0, thetaq, wigner) + 2 * q2y * q3y * Yket(2, 2, 2, 0, thetaq, wigner)) -
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -1, 2, 0, thetaq, wigner) + Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 2, 0, thetaq, wigner) -
                                                                                    2 * ((q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, 1, 2, 0, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 2, 2, 0, thetaq, wigner))))) /
              (4. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 16)
    {
        mtx = (Ftpe1 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 2, thetaq, wigner) - 2 * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) - Sqrt(6) * q2z * q3x * Yket(2, 0, 2, 2, thetaq, wigner) +
                             Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3z * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 2, 2, thetaq, wigner) + 2 * q2x * q3x * Yket(2, 1, 2, 2, thetaq, wigner) -
                             Complex(0, 2) * q2y * q3x * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 2) * q2x * q3y * Yket(2, 1, 2, 2, thetaq, wigner) - 2 * q2y * q3y * Yket(2, 1, 2, 2, thetaq, wigner)) -
                        Sqrt(2) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -2, 2, 2, thetaq, wigner) + (q2z * q3x + Complex(0, 1) * q2z * q3y + q2x * q3z + Complex(0, 1) * q2y * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) + q2z * q3x * Yket(2, 1, 2, 2, thetaq, wigner) -
                             Complex(0, 1) * q2z * q3y * Yket(2, 1, 2, 2, thetaq, wigner) + q2x * q3z * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Yket(2, 1, 2, 2, thetaq, wigner) - 2 * q2x * q3x * Yket(2, 2, 2, 2, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Yket(2, 2, 2, 2, thetaq, wigner) +
                             Complex(0, 2) * q2x * q3y * Yket(2, 2, 2, 2, thetaq, wigner) + 2 * q2y * q3y * Yket(2, 2, 2, 2, thetaq, wigner)) -
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * Yket(2, -1, 2, 2, thetaq, wigner) + Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, 0, 2, 2, thetaq, wigner) -
                                                                                    2 * ((q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(2, 1, 2, 2, thetaq, wigner) + (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(2, 2, 2, 2, thetaq, wigner))))) /
              (4. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 2)
    {
        mtx = -((Ftpe1 *
                 (Sqrt(3) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(2) * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2z * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(3) * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * q2x * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2y * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 0, 0, thetaq, wigner)) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 3)
    {
        mtx =
            -0.05 *
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      2 * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 0, 2, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 0, 2, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) -
                      Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 0, 2, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 0, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 0, 2, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 0, 2, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 0, 2, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 0, 2, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 5)
    {
        mtx = (Ftpe1 * (2 * Sqrt(3) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Sqrt(2) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 4) * q2y * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * q2x * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Sqrt(3) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, 2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + Sqrt(2) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner)) -
                        Sqrt(2) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2y * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2x * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        2 * Sqrt(3) * q2z * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2z * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        2 * Sqrt(3) * q2x * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2y * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              (2. * Sqrt(10));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 12)
    {
        mtx = -((Ftpe1 *
                 (Sqrt(3) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(2) * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2z * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(3) * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * q2x * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2y * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 2, 2, thetaq, wigner)) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 13)
    {
        mtx = (Ftpe1 * (-(Sqrt(3) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner)) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2x * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        2 * q2x * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * q2x * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 2 * q2y * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - 2 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 0, 2, 2, thetaq, wigner)) -
                        Sqrt(3) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        2 * Sqrt(2) * q2z * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2z * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2x * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                     Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              (2. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 14)
    {
        mtx = -0.25 *
              (Ftpe1 * (Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        2 * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * q2y * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * q2x * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 2 * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        2 * Sqrt(2) * q2x * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * q2y * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(2) * q2x * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(2) * q2y * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        2 * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        2 * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        4 * q2z * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * q2z * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * q2z * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * q2x * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * q2y * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                     2 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              Sqrt(5);
        return mtx;
    }
    else if (idxbra == 15 && idxket == 15)
    {
        mtx =
            -0.05 *
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      2 * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 2, 0, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 2, 0, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) -
                      Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 0, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 2, 0, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 2, 0, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 2, 0, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 2, 0, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 2, 0, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 16)
    {
        mtx =
            -0.05 *
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      2 * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 2, 2, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 2, 2, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) -
                      Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 2, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 2, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 2, 2, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 2, 2, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 2, 2, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 2, 2, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 16 && idxket == 2)
    {
        mtx = -((Ftpe1 *
                 (Sqrt(3) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(2) * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2z * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(3) * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * q2x * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2y * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 0, 0, thetaq, wigner)) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 16 && idxket == 3)
    {
        mtx =
            -0.05 *
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) +
                      2 * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 0, 2, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 0, 2, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 0, 2, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) -
                      Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 0, 2, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 0, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 0, 2, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 0, 2, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 0, 2, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 0, 2, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 16 && idxket == 5)
    {
        mtx = (Ftpe1 * (2 * Sqrt(3) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Sqrt(2) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 4) * q2y * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * q2x * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Sqrt(3) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, 2) * (q2y * q3x - q2x * q3y) * Yket(1, -1, 1, 1, thetaq, wigner) + Sqrt(2) * (-(q2z * q3x) + Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 1, 1, thetaq, wigner)) -
                        Sqrt(2) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2y * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2x * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        2 * Sqrt(3) * q2z * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2z * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        2 * Sqrt(3) * q2x * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2y * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              (2. * Sqrt(10));
        return mtx;
    }
    else if (idxbra == 16 && idxket == 12)
    {
        mtx = -((Ftpe1 *
                 (Sqrt(3) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(2) * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2z * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Sqrt(3) * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                  Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(3) * q2x * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 1) * Sqrt(3) * q2x * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(3) * q2y * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 2, 2, thetaq, wigner)) /
                Sqrt(5));
        return mtx;
    }
    else if (idxbra == 16 && idxket == 13)
    {
        mtx = (Ftpe1 * (-(Sqrt(3) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner)) - Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2x * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        2 * q2x * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * q2y * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * q2x * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 2 * q2y * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - 2 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 0, 2, 2, thetaq, wigner)) -
                        Sqrt(3) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        2 * Sqrt(2) * q2z * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2z * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2x * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                     Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              (2. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 16 && idxket == 14)
    {
        mtx = -0.25 *
              (Ftpe1 * (Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        2 * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * q2y * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * q2x * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 2 * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(2) * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        2 * Sqrt(2) * q2x * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * q2y * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(2) * q2x * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(2) * q2y * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        2 * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        2 * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        4 * q2z * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * q2z * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * q2z * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * q2x * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * q2y * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x * q3x + q2y * q3y - 2 * q2z * q3z) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y + q2x * q3z - Complex(0, 1) * q2y * q3z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                     2 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              Sqrt(5);
        return mtx;
    }
    else if (idxbra == 16 && idxket == 15)
    {
        mtx =
            -0.05 *
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) +
                      2 * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 2, 0, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 2, 0, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 2, 0, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) -
                      Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 0, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 2, 0, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 2, 0, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 2, 0, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 2, 0, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 2, 0, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 16 && idxket == 16)
    {
        mtx =
            -0.05 *
            (Ftpe1 * (2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) -
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) -
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                      Complex(0, 6) * q2x * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) - 6 * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                      8 * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + 8 * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                      4 * q2z * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                      2 * Sqrt(6) * q2y * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) +
                      2 * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * q2z * q3z * Yket(2, -2, 2, 2, thetaq, wigner) - 3 * (q2z * (q3x - Complex(0, 1) * q3y) + (q2x - Complex(0, 1) * q2y) * q3z) * Yket(2, -1, 2, 2, thetaq, wigner) + Sqrt(6) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * Yket(2, 0, 2, 2, thetaq, wigner)) +
                      Sqrt(6) * q2z * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      Sqrt(6) * q2x * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      6 * q2x * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + 6 * q2y * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) -
                      Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (6 * (q2z * (q3x + Complex(0, 1) * q3y) + (q2x + Complex(0, 1) * q2y) * q3z) * Yket(2, -2, 2, 2, thetaq, wigner) - 6 * (q2x * q3x + q2y * q3y) * Yket(2, -1, 2, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * Yket(2, 0, 2, 2, thetaq, wigner) -
                                                                                   Complex(0, 1) * Sqrt(6) * q2z * q3y * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2x * q3z * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * Yket(2, 0, 2, 2, thetaq, wigner) -
                                                                                   6 * q2x * q3x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3y * Yket(2, 1, 2, 2, thetaq, wigner) + 6 * q2y * q3y * Yket(2, 1, 2, 2, thetaq, wigner)) +
                      2 * Sqrt(6) * q2x * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) -
                      Complex(0, 2) * Sqrt(6) * q2x * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) +
                      6 * q2z * q3x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * q2z * q3y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) +
                      6 * q2x * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) +
                      12 * q2z * q3z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner)));
        return mtx;
    }

    else
    {
        // std::cerr << "unkown channel index: (" << idxbra << "," << idxket << ") !" << std::endl;
        return 0;
    }
}
} // end namespace aPWD3_c1
