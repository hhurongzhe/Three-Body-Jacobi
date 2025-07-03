#include <iostream>
#include <complex>
#include <tuple>
#include <vector>
#include <cmath>

#include "constants.hpp"
#include "./numlib/spherical_harmonics.hpp"
#include "aPWD3_part_c4.hpp"

namespace aPWD3_c4
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
double get_q1x(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = -(qmag * Sin(thetaq)) + qpmag * Cos(phiqp) * Sin(thetaqp);
    return value;
}
double get_q1y(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = qpmag * Sin(phiqp) * Sin(thetaqp);
    return value;
}
double get_q1z(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = -(qmag * Cos(thetaq)) + qpmag * Cos(thetaqp);
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
double get_q4x(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = ppmag * qmag * Cos(thetaq) * Sin(phipp) * Sin(thetapp) - ppmag * qpmag * Cos(thetaqp) * Sin(phipp) * Sin(thetapp) - pmag * qpmag * Sin(phiqp) * Sin(thetaqp) + ppmag * qpmag * Cos(thetapp) * Sin(phiqp) * Sin(thetaqp);
    return value;
}
double get_q4y(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = -(ppmag * qmag * Cos(phipp) * Cos(thetaq) * Sin(thetapp)) + ppmag * qpmag * Cos(phipp) * Cos(thetaqp) * Sin(thetapp) - pmag * qmag * Sin(thetaq) + ppmag * qmag * Cos(thetapp) * Sin(thetaq) + pmag * qpmag * Cos(phiqp) * Sin(thetaqp) - ppmag * qpmag * Cos(phiqp) * Cos(thetapp) * Sin(thetaqp);
    return value;
}
double get_q4z(double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
{
    double value = -(ppmag * qmag * Sin(phipp) * Sin(thetapp) * Sin(thetaq)) + ppmag * qpmag * Cos(phiqp) * Sin(phipp) * Sin(thetapp) * Sin(thetaqp) - ppmag * qpmag * Cos(phipp) * Sin(phiqp) * Sin(thetapp) * Sin(thetaqp);
    return value;
}

// spin projection of 3N forces.
std::complex<double> Gt(int idxbra, int idxket, double ppmag, double qpmag, double pmag, double qmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp, util::WignerSymbols &wigner)
{
    double q1x = get_q1x(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q1y = get_q1y(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q1z = get_q1z(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q2x = get_q2x(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q2y = get_q2y(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q2z = get_q2z(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q3x = get_q3x(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q3y = get_q3y(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q3z = get_q3z(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q4x = get_q4x(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q4y = get_q4y(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q4z = get_q4z(ppmag, qpmag, pmag, qmag, thetaq, thetapp, phipp, thetaqp, phiqp);
    double q2q2 = q2x * q2x + q2y * q2y + q2z * q2z;
    double q3q3 = q3x * q3x + q3y * q3y + q3z * q3z;
    double Ftpe3 = Power(gA / (2.0 * fpi), 2) / (q2q2 + Mpi * Mpi) / (q3q3 + Mpi * Mpi) / (fpi * fpi);

    std::complex<double> mtx;
    if (idxbra == -1 && idxket == -1)
    {
        return 0.0;
    }
    else if (idxbra == 1 && idxket == 2)
    {
        mtx = 2 * Ftpe3 * (-(q2z * q3y * q4x) + q2y * q3z * q4x + q2z * q3x * q4y - q2x * q3z * q4y - q2y * q3x * q4z + q2x * q3y * q4z) * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 1 && idxket == 3)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               (3 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 0, 2, thetaq, wigner) +
                3 * (q2y * q3x * (q4x + Complex(0, 1) * q4y) - q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(2, -1, 0, 2, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) +
                Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + 2 * Sqrt(6) * q2y * q3x * q4z * Yket(2, 0, 0, 2, thetaq, wigner) -
                2 * Sqrt(6) * q2x * q3y * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - 3 * q2y * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + 3 * q2x * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) -
                Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 3) * q2z * q3x * q4z * Yket(2, 1, 0, 2, thetaq, wigner) - 3 * q2z * q3y * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * q2x * q3z * q4z * Yket(2, 1, 0, 2, thetaq, wigner) +
                3 * q2y * q3z * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 3 * q2z * q3y * q4x * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) -
                3 * q2y * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 3 * q2z * q3x * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - 3 * q2x * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) +
                Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner))) /
              Sqrt(15);
        return mtx;
    }
    else if (idxbra == 1 && idxket == 5)
    {
        mtx = Complex(0, 1) * Ftpe3 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
              (Sqrt(2) * (q4x + Complex(0, 1) * q4y) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * q4z * Yket(1, 0, 1, 1, thetaq, wigner) - Sqrt(2) * (q4x - Complex(0, 1) * q4y) * Yket(1, 1, 1, 1, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 1 && idxket == 12)
    {
        mtx = 2 * Ftpe3 * (-(q2z * q3y * q4x) + q2y * q3z * q4x + q2z * q3x * q4y - q2x * q3z * q4y - q2y * q3x * q4z + q2x * q3y * q4z) * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 1 && idxket == 13)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                Complex(0, 2) * (q3z * (q2x * q4x + q2y * q4y) - q2z * (q3x * q4x + q3y * q4y)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 1 && idxket == 14)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                Complex(0, 1) * Sqrt(2) * (q3z * (q2x * q4x + q2y * q4y) - q2z * (q3x * q4x + q3y * q4y)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 1 && idxket == 15)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               (3 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 2, 0, thetaq, wigner) +
                3 * (q2y * q3x * (q4x + Complex(0, 1) * q4y) - q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(2, -1, 2, 0, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) +
                Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + 2 * Sqrt(6) * q2y * q3x * q4z * Yket(2, 0, 2, 0, thetaq, wigner) -
                2 * Sqrt(6) * q2x * q3y * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - 3 * q2y * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + 3 * q2x * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) -
                Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 3) * q2z * q3x * q4z * Yket(2, 1, 2, 0, thetaq, wigner) - 3 * q2z * q3y * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * q2x * q3z * q4z * Yket(2, 1, 2, 0, thetaq, wigner) +
                3 * q2y * q3z * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 3 * q2z * q3y * q4x * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) -
                3 * q2y * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 3 * q2z * q3x * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - 3 * q2x * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) +
                Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner))) /
              Sqrt(15);
        return mtx;
    }
    else if (idxbra == 1 && idxket == 16)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               (3 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 2, 2, thetaq, wigner) +
                3 * (q2y * q3x * (q4x + Complex(0, 1) * q4y) - q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(2, -1, 2, 2, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) +
                Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(6) * q2y * q3x * q4z * Yket(2, 0, 2, 2, thetaq, wigner) -
                2 * Sqrt(6) * q2x * q3y * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - 3 * q2y * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + 3 * q2x * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) -
                Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 3) * q2z * q3x * q4z * Yket(2, 1, 2, 2, thetaq, wigner) - 3 * q2z * q3y * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2x * q3z * q4z * Yket(2, 1, 2, 2, thetaq, wigner) +
                3 * q2y * q3z * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 3 * q2z * q3y * q4x * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) -
                3 * q2y * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 3 * q2z * q3x * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - 3 * q2x * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) +
                Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner))) /
              Sqrt(15);
        return mtx;
    }
    else if (idxbra == 2 && idxket == 1)
    {
        mtx = -2 * Ftpe3 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 2 && idxket == 7)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               (Sqrt(2) *
                    (q2z * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, -3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y - Complex(0, 2) * q3y * q4y - Complex(0, 2) * q3z * q4z) +
                     q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Yket(1, -1, 1, 1, thetaq, wigner) -
                Complex(0, 2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                Sqrt(2) *
                    (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                     q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Yket(1, 1, 1, 1, thetaq, wigner))) /
              3.;
        return mtx;
    }
    else if (idxbra == 2 && idxket == 8)
    {
        mtx =
            -0.3333333333333333 *
            (Ftpe3 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
             ((q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) -
              Complex(0, 1) * Sqrt(2) * (q2z * q3x * q4x + q2x * q3z * q4x + q2z * q3y * q4y + q2y * q3z * q4y - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, 0, 1, 1, thetaq, wigner) +
              (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                  Yket(1, 1, 1, 1, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 2 && idxket == 9)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Yket(2, -2, 1, 1, thetaq, wigner) +
                (q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (-(q3y * q4x) + 2 * q3x * q4y + Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (-(q3x * q4x) - Complex(0, 2) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3z * q4z)) *
                    Yket(2, -1, 1, 1, thetaq, wigner) -
                Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + q2y * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) +
                q2x * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 2) * q2y * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * q2z * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - 2 * q2x * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) +
                Complex(0, 1) * q2y * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 1) * q2x * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + 2 * q2z * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * q2z * q3x * q4z * Yket(2, 1, 1, 1, thetaq, wigner) -
                q2z * q3y * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * q2x * q3z * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - q2y * q3z * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * q2z * q3x * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - q2z * q3y * q4x * Yket(2, 2, 1, 1, thetaq, wigner) -
                Complex(0, 1) * q2x * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - q2y * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - q2z * q3x * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 1) * q2z * q3y * q4y * Yket(2, 2, 1, 1, thetaq, wigner) - q2x * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) +
                Complex(0, 1) * q2y * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 2) * q2x * q3x * q4z * Yket(2, 2, 1, 1, thetaq, wigner) + 2 * q2y * q3x * q4z * Yket(2, 2, 1, 1, thetaq, wigner) + 2 * q2x * q3y * q4z * Yket(2, 2, 1, 1, thetaq, wigner) -
                Complex(0, 2) * q2y * q3y * q4z * Yket(2, 2, 1, 1, thetaq, wigner))) /
              Sqrt(5);
        return mtx;
    }
    else if (idxbra == 2 && idxket == 10)
    {
        mtx = -2 * Ftpe3 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 2 && idxket == 11)
    {
        mtx = -((Ftpe3 * Ybra(0, 0, 0, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                 (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                  Complex(0, 2) * (q3z * (q2x * q4x + q2y * q4y) - q2z * (q3x * q4x + q3y * q4y)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                  Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner))) /
                Sqrt(3));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 1)
    {
        mtx = (Ftpe3 *
               (3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * (Complex(0, 1) * q4x + q4y) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z - q2y * q3z * q4z - q2x * (q3y * (q4x - Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 4 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 4 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 3 && idxket == 6)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (Complex(0, 2) * q3z * q4x + 2 * q3z * q4y - Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) - q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z)) *
                    Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 2) * Sqrt(2) * q2z * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2z * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              Sqrt(10);
        return mtx;
    }
    else if (idxbra == 3 && idxket == 7)
    {
        mtx = (Ftpe3 * (2 * Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 4) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        4 * Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2z * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 3 * q2z * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2x * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 3 * q2y * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        3 * q2z * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        3 * q2x * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2y * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 6 * q2y * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        6 * q2x * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 4) * Sqrt(3) * q2x * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2y * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * ((Complex(0, 1) * q2x + q2y) * q3z * q4z + q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                        2 * Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 4) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        4 * Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2z * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 9 * q2z * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2x * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 9 * q2y * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        9 * q2z * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 3) * q2z * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        9 * q2x * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 3) * q2y * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2z * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 12 * q2z * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3x * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2z * q3y * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2y * q3z * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        3 * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, 1) * (q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (Complex(0, -2) * q2y * q3y * q4x + q2y * q3x * (q4x + Complex(0, 1) * q4y) + q2x * q3y * (q4x + Complex(0, 1) * q4y) - 2 * q2x * q3x * q4y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                             (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (3. * Sqrt(30));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 8)
    {
        mtx = (Ftpe3 * (Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        3 * Sqrt(6) * q2y * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 5 * Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2z * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 12 * q2z * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2x * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 12 * q2y * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        12 * q2z * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2z * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        12 * q2x * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2y * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2x * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 12 * q2y * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        12 * q2x * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2y * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 18) * q2x * q3x * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 18 * q2y * q3x * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        18 * q2x * q3y * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2y * q3y * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        18 * q2x * q3x * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2y * q3x * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 18) * q2x * q3y * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 18 * q2y * q3y * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 8) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 8) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 8) * Sqrt(3) * q2x * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2y * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 12) * Sqrt(3) * q2z * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2x * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2z * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 9 * Sqrt(2) * q2y * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        9 * Sqrt(2) * q2z * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        9 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        9 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 3 * Sqrt(6) * q2y * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        5 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        5 * Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        5 * Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2y * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2z * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2z * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3x * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2z * q3y * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2y * q3z * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        6 * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((q2z * q3z * (Complex(0, 1) * q4x + q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z + (Complex(0, 1) * q2x + q2y) * q3z * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) -
                             Complex(0, 1) * Sqrt(2) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                             3 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (Complex(0, 1) * q4x + q4y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        3 * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, -2) * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y - 3 * q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) *
                                 (Complex(0, -3) * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2y * (q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + 3 * q3y * q4y - 3 * q3z * q4z) +
                                  q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y + q3x * (Complex(0, 3) * q4x + q4y) - Complex(0, 3) * q3z * q4z)) *
                                 Yket(1, 0, 1, 1, thetaq, wigner) +
                             4 * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 3 && idxket == 9)
    {
        mtx = (Complex(0, -0.1) * Ftpe3 *
               (12 * q2z * q3z * q4z * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                3 * q2z * q3z * q4z * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                3 * q2z * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                12 * q2z * q3z * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 3 && idxket == 10)
    {
        mtx = (Ftpe3 *
               (3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * (Complex(0, 1) * q4x + q4y) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z - q2y * q3z * q4z - q2x * (q3y * (q4x - Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 4 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 4 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 3 && idxket == 11)
    {
        mtx = (Ftpe3 * (2 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - 2 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(2) * q2z * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(3) * q2z * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(3) * q2y * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2z * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(3) * q2x * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * q2z * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * q2x * q3z * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * q2z * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * q2y * q3z * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2y * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2y * q3z * q4x * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Sqrt(6) * q2z * q3x * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * q4y * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(3) * (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * Ybra(2, -2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * q4z * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * (q4x - Complex(0, 1) * q4y) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        2 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2z * q3y * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3z * q4z * Ybra(2, 0, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Sqrt(3) * q2z * q3y * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(3) * q2y * q3z * q4x * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(3) * q2z * q3x * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        2 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * Sqrt(3) * q2z * q3y * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 2 * Sqrt(3) * q2y * q3z * q4z * Ybra(2, 2, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(3) * Ybra(2, -1, 0, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((q2z * (Complex(0, 1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) - q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x + q3z * q4y + 2 * q3y * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                             (q4x - Complex(0, 1) * q4y) * (Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 2, 2, thetaq, wigner) + (q2z * (Complex(0, 1) * q3x + q3y) - Complex(0, 1) * q2x * q3z - q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner))))) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 4 && idxket == 6)
    {
        mtx = 2 * Ftpe3 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 1, 1, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 4 && idxket == 7)
    {
        mtx = -((Ftpe3 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                 (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                  Complex(0, 2) * (q3z * (q2x * q4x + q2y * q4y) - q2z * (q3x * q4x + q3y * q4y)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                  Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner))) /
                Sqrt(3));
        return mtx;
    }
    else if (idxbra == 4 && idxket == 8)
    {
        mtx = -((Ftpe3 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                 ((q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                  Complex(0, 1) * Sqrt(2) * (q3z * (q2x * q4x + q2y * q4y) - q2z * (q3x * q4x + q3y * q4y)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                  (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner))) /
                Sqrt(3));
        return mtx;
    }
    else if (idxbra == 4 && idxket == 9)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
               (3 * (Complex(0, 1) * q2z * q3x - q2z * q3y - Complex(0, 1) * q2x * q3z + q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 1, 1, thetaq, wigner) +
                3 * (-(q2y * q3x * (q4x + Complex(0, 1) * q4y)) + q2x * q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z - Complex(0, 1) * q2x * q3z * q4z + q2y * q3z * q4z) * Yket(2, -1, 1, 1, thetaq, wigner) + Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 1, 1, thetaq, wigner) - Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 1, 1, thetaq, wigner) - 2 * Sqrt(6) * q2y * q3x * q4z * Yket(2, 0, 1, 1, thetaq, wigner) +
                2 * Sqrt(6) * q2x * q3y * q4z * Yket(2, 0, 1, 1, thetaq, wigner) + 3 * q2y * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - 3 * q2x * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) +
                Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 3) * q2z * q3x * q4z * Yket(2, 1, 1, 1, thetaq, wigner) + 3 * q2z * q3y * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4z * Yket(2, 1, 1, 1, thetaq, wigner) -
                3 * q2y * q3z * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - 3 * q2z * q3y * q4x * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) +
                3 * q2y * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - 3 * q2z * q3x * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + 3 * q2x * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) -
                Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner))) /
              Sqrt(15);
        return mtx;
    }
    else if (idxbra == 4 && idxket == 11)
    {
        mtx = Ftpe3 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
              (Sqrt(2) * (Complex(0, -1) * q4x + q4y) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * q4z * Yket(1, 0, 2, 2, thetaq, wigner) + Sqrt(2) * (Complex(0, 1) * q4x + q4y) * Yket(1, 1, 2, 2, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 5 && idxket == 1)
    {
        mtx = Complex(0, -1) * Ftpe3 * (q2x * q3x + q2y * q3y + q2z * q3z) *
              (Sqrt(2) * (q4x - Complex(0, 1) * q4y) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * q4z * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * (q4x + Complex(0, 1) * q4y) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) * Yket(0, 0, 0, 0, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 5 && idxket == 6)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (-(q3z * (q2x * q4x + q2y * q4y)) + q2z * (q3x * q4x + q3y * q4y)) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 5 && idxket == 7)
    {
        mtx = (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + q2y * (-(q3z * q4x) + Complex(0, 1) * q3z * q4y + q3x * q4z) + q2x * (Complex(0, 1) * q3z * q4x + q3z * q4y - q3y * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                                                                                     Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * q4z + q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) - q2y * (q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                    2 * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) - q2y * q3z * (q4x + Complex(0, 1) * q4y) + q2x * q3z * (Complex(0, -1) * q4x + q4y) + q2y * q3x * q4z - q2x * q3y * q4z) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * q4z + q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) - q2y * (q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                                                                                    2 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                    Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              3.;
        return mtx;
    }
    else if (idxbra == 5 && idxket == 8)
    {
        mtx = (Ftpe3 *
               (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) * (Complex(0, 1) * q4x + q4y) + q2y * (q3z * q4x - Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x - q3z * q4y - 2 * q3y * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                                                                             2 * (q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2x * (2 * q3y * q4x - Complex(0, 2) * q3y * q4y + Complex(0, 1) * q3z * q4z) + q2y * (-2 * q3x * q4x + Complex(0, 2) * q3x * q4y + q3z * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                             3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * (Complex(0, 1) * q4x + q4y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * Sqrt(2) * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(1, -1, 1, 1, thetaq, wigner) +
                                                                            2 * (2 * q2y * q3x * (q4x + Complex(0, 1) * q4y) - 2 * q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                            Sqrt(2) * (q2z * (q3x - Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) + q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, 1) * q3z * q4x - q3z * q4y - 2 * q3y * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)) -
                2 * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                    ((q2y * q3x * (q4x + Complex(0, 1) * q4y) - q2x * q3y * (q4x + Complex(0, 1) * q4y) + 2 * q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 2) * q2x * q3z * q4z - 2 * q2y * q3z * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) +
                     Sqrt(2) * (-(q2z * q3y * q4x) + q2y * q3z * q4x + q2z * q3x * q4y - q2x * q3z * q4y + 2 * q2y * q3x * q4z - 2 * q2x * q3y * q4z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                     (q2x * q3y * q4x - Complex(0, 1) * q2x * q3y * q4y - Complex(0, 2) * q2z * q3x * q4z - 2 * q2z * q3y * q4z + Complex(0, 2) * q2x * q3z * q4z + q2y * (-(q3x * q4x) + Complex(0, 1) * q3x * q4y + 2 * q3z * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              6.;
        return mtx;
    }
    else if (idxbra == 5 && idxket == 9)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * Sqrt(2) * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * q4z * Yket(2, -2, 1, 1, thetaq, wigner) +
                           3 * Sqrt(2) * (Complex(0, -1) * q2x * q3z * q4x + q2y * q3z * q4x - q2x * q3z * q4y - Complex(0, 1) * q2y * q3z * q4y + q2z * (q3x + Complex(0, 1) * q3y) * (Complex(0, 1) * q4x + q4y) + 2 * q2y * q3x * q4z - 2 * q2x * q3y * q4z) * Yket(2, -1, 1, 1, thetaq, wigner) -
                           4 * Sqrt(3) * q2y * q3x * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + 4 * Sqrt(3) * q2x * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2y * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2x * q3y * q4y * Yket(2, 0, 1, 1, thetaq, wigner) -
                           Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Yket(2, 0, 1, 1, thetaq, wigner) - 2 * Sqrt(3) * q2z * q3y * q4z * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Yket(2, 0, 1, 1, thetaq, wigner) + 2 * Sqrt(3) * q2y * q3z * q4z * Yket(2, 0, 1, 1, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) +
                           3 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner)) +
                      Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * (Complex(0, 1) * q2z * q3x - q2z * q3y - Complex(0, 1) * q2x * q3z + q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 1, 1, thetaq, wigner) - 6 * (q2y * q3x - q2x * q3y) * (q4x + Complex(0, 1) * q4y) * Yket(2, -1, 1, 1, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(6) * q2z * q3x * q4x * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2x * q3z * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2z * q3y * q4y * Yket(2, 0, 1, 1, thetaq, wigner) -
                           Complex(0, 2) * Sqrt(6) * q2y * q3z * q4y * Yket(2, 0, 1, 1, thetaq, wigner) - 6 * q2y * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + 6 * q2x * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2y * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) -
                           Complex(0, 6) * q2x * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2z * q3x * q4x * Yket(2, 2, 1, 1, thetaq, wigner) + 6 * q2z * q3y * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 6) * q2x * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) -
                           6 * q2y * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) + 6 * q2z * q3x * q4y * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3y * q4y * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * q2x * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) +
                           Complex(0, 6) * q2y * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner)) +
                      Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * Sqrt(2) * (Complex(0, 1) * q2z * q3x - q2z * q3y - Complex(0, 1) * q2x * q3z + q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -1, 1, 1, thetaq, wigner) +
                           2 * Sqrt(3) * (-2 * q2y * q3x * (q4x + Complex(0, 1) * q4y) + 2 * q2x * q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z - Complex(0, 1) * q2x * q3z * q4z + q2y * q3z * q4z) * Yket(2, 0, 1, 1, thetaq, wigner) +
                           3 * Sqrt(2) *
                               ((q2z * (Complex(0, 1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) - q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x + q3z * q4y + 2 * q3y * q4z)) * Yket(2, 1, 1, 1, thetaq, wigner) +
                                2 * (q2z * (Complex(0, 1) * q3x + q3y) - Complex(0, 1) * q2x * q3z - q2y * q3z) * q4z * Yket(2, 2, 1, 1, thetaq, wigner))))) /
            (6. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 5 && idxket == 10)
    {
        mtx = Complex(0, -1) * Ftpe3 * (q2x * q3x + q2y * q3y + q2z * q3z) *
              (Sqrt(2) * (q4x - Complex(0, 1) * q4y) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * q4z * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * (q4x + Complex(0, 1) * q4y) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) * Yket(0, 0, 2, 2, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 5 && idxket == 11)
    {
        mtx = (Ftpe3 * (q2x * q3x + q2y * q3y + q2z * q3z) *
               (Complex(0, 1) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * q4z * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * (q4x - Complex(0, 1) * q4y) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                Sqrt(2) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * ((Complex(0, -1) * q4x + q4y) * Yket(1, -1, 2, 2, thetaq, wigner) + (Complex(0, -1) * q4x - q4y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (Complex(0, -1) * q4x + q4y) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * q4z * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 6 && idxket == 3)
    {
        mtx =
            (Ftpe3 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
             ((q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 0, 2, thetaq, wigner) +
              (q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) - q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(2, -1, 0, 2, thetaq, wigner) +
              Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - q2y * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) -
              q2x * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 2) * q2y * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 2) * q2z * q3z * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + 2 * q2x * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) -
              Complex(0, 1) * q2y * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 1) * q2x * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - 2 * q2z * q3z * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 1) * q2z * q3x * q4z * Yket(2, 1, 0, 2, thetaq, wigner) +
              q2z * q3y * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 1) * q2x * q3z * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + q2y * q3z * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 1) * q2z * q3x * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + q2z * q3y * q4x * Yket(2, 2, 0, 2, thetaq, wigner) +
              Complex(0, 1) * q2x * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + q2y * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + q2z * q3x * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 1) * q2z * q3y * q4y * Yket(2, 2, 0, 2, thetaq, wigner) + q2x * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) -
              Complex(0, 1) * q2y * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 2) * q2x * q3x * q4z * Yket(2, 2, 0, 2, thetaq, wigner) - 2 * q2y * q3x * q4z * Yket(2, 2, 0, 2, thetaq, wigner) - 2 * q2x * q3y * q4z * Yket(2, 2, 0, 2, thetaq, wigner) +
              Complex(0, 2) * q2y * q3y * q4z * Yket(2, 2, 0, 2, thetaq, wigner))) /
            Sqrt(5);
        return mtx;
    }
    else if (idxbra == 6 && idxket == 4)
    {
        mtx = 2 * Ftpe3 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 1, 1, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 6 && idxket == 5)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
               (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                Complex(0, 2) * (q3z * (q2x * q4x + q2y * q4y) - q2z * (q3x * q4x + q3y * q4y)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 6 && idxket == 13)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
               (Sqrt(2) *
                    (q2z * (Complex(0, -1) * q3z * q4x + q3z * q4y + Complex(0, 2) * q3x * q4z - 2 * q3y * q4z) + q2x * (-2 * q3y * q4x + Complex(0, 2) * q3y * q4y + q3x * (Complex(0, 3) * q4x + q4y) + Complex(0, 2) * q3z * q4z) -
                     q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Yket(1, -1, 2, 2, thetaq, wigner) +
                Complex(0, 2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner) -
                Sqrt(2) *
                    (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                     q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Yket(1, 1, 2, 2, thetaq, wigner))) /
              3.;
        return mtx;
    }
    else if (idxbra == 6 && idxket == 14)
    {
        mtx =
            (Ftpe3 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
             ((q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) -
              Complex(0, 1) * Sqrt(2) * (q2z * q3x * q4x + q2x * q3z * q4x + q2z * q3y * q4y + q2y * q3z * q4y - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, 0, 2, 2, thetaq, wigner) +
              (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                  Yket(1, 1, 2, 2, thetaq, wigner))) /
            3.;
        return mtx;
    }
    else if (idxbra == 6 && idxket == 15)
    {
        mtx =
            (Ftpe3 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
             ((q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 2, 0, thetaq, wigner) +
              (q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) - q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(2, -1, 2, 0, thetaq, wigner) +
              Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - q2y * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) -
              q2x * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 2) * q2y * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 2) * q2z * q3z * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + 2 * q2x * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) -
              Complex(0, 1) * q2y * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 1) * q2x * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - 2 * q2z * q3z * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 1) * q2z * q3x * q4z * Yket(2, 1, 2, 0, thetaq, wigner) +
              q2z * q3y * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 1) * q2x * q3z * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + q2y * q3z * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 1) * q2z * q3x * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + q2z * q3y * q4x * Yket(2, 2, 2, 0, thetaq, wigner) +
              Complex(0, 1) * q2x * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + q2y * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + q2z * q3x * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 1) * q2z * q3y * q4y * Yket(2, 2, 2, 0, thetaq, wigner) + q2x * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) -
              Complex(0, 1) * q2y * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 2) * q2x * q3x * q4z * Yket(2, 2, 2, 0, thetaq, wigner) - 2 * q2y * q3x * q4z * Yket(2, 2, 2, 0, thetaq, wigner) - 2 * q2x * q3y * q4z * Yket(2, 2, 2, 0, thetaq, wigner) +
              Complex(0, 2) * q2y * q3y * q4z * Yket(2, 2, 2, 0, thetaq, wigner))) /
            Sqrt(5);
        return mtx;
    }
    else if (idxbra == 6 && idxket == 16)
    {
        mtx =
            (Ftpe3 * Ybra(0, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
             ((q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 2, 2, thetaq, wigner) +
              (q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) - q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(2, -1, 2, 2, thetaq, wigner) +
              Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - q2y * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) -
              q2x * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 2) * q2y * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 2) * q2z * q3z * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + 2 * q2x * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) -
              Complex(0, 1) * q2y * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * q2x * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - 2 * q2z * q3z * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * q2z * q3x * q4z * Yket(2, 1, 2, 2, thetaq, wigner) +
              q2z * q3y * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * q2x * q3z * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + q2y * q3z * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * q2z * q3x * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + q2z * q3y * q4x * Yket(2, 2, 2, 2, thetaq, wigner) +
              Complex(0, 1) * q2x * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + q2y * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + q2z * q3x * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 1) * q2z * q3y * q4y * Yket(2, 2, 2, 2, thetaq, wigner) + q2x * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) -
              Complex(0, 1) * q2y * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 2) * q2x * q3x * q4z * Yket(2, 2, 2, 2, thetaq, wigner) - 2 * q2y * q3x * q4z * Yket(2, 2, 2, 2, thetaq, wigner) - 2 * q2x * q3y * q4z * Yket(2, 2, 2, 2, thetaq, wigner) +
              Complex(0, 2) * q2y * q3y * q4z * Yket(2, 2, 2, 2, thetaq, wigner))) /
            Sqrt(5);
        return mtx;
    }
    else if (idxbra == 7 && idxket == 2)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) *
                    (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                     q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) *
                    (q2z * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, -3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y - Complex(0, 2) * q3y * q4y - Complex(0, 2) * q3z * q4z) +
                     q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              3.;
        return mtx;
    }
    else if (idxbra == 7 && idxket == 3)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * Sqrt(2) * ((Complex(0, -1) * q2x + q2y) * q3z * q4z + q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z)) * Yket(2, -2, 0, 2, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * (q2z * q3x * (q4x + Complex(0, 3) * q4y) + q2x * q3z * (q4x + Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, -3) * q4x + q4y) + q2y * q3z * (Complex(0, -3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(2, -1, 0, 2, thetaq, wigner) +
                           4 * Sqrt(3) * q2y * q3x * q4x * Yket(2, 0, 0, 2, thetaq, wigner) + 4 * Sqrt(3) * q2x * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2y * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2z * q3z * q4x * Yket(2, 0, 0, 2, thetaq, wigner) -
                           8 * Sqrt(3) * q2x * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2y * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2x * q3y * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + 4 * Sqrt(3) * q2z * q3z * q4y * Yket(2, 0, 0, 2, thetaq, wigner) -
                           Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - 2 * Sqrt(3) * q2z * q3y * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - 2 * Sqrt(3) * q2y * q3z * q4z * Yket(2, 0, 0, 2, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 1, 0, 2, thetaq, wigner) -
                           3 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 1, 0, 2, thetaq, wigner) +
                           Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3y * q4z * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Yket(2, 1, 0, 2, thetaq, wigner)) +
                      2 * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 0, 2, thetaq, wigner) +
                           3 * (Complex(0, 2) * q2y * q3y * q4x + q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2x * q3y * (q4x - Complex(0, 1) * q4y) - 2 * q2x * q3x * q4y) * Yket(2, -1, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3x * q4x * Yket(2, 0, 0, 2, thetaq, wigner) -
                           Complex(0, 1) * Sqrt(6) * q2x * q3z * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * q4y * Yket(2, 0, 0, 2, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(6) * q2x * q3x * q4z * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3y * q4z * Yket(2, 0, 0, 2, thetaq, wigner) + 3 * q2y * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + 3 * q2x * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) -
                           Complex(0, 6) * q2y * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - 6 * q2x * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) -
                           Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 0, 2, thetaq, wigner) - 3 * q2z * q3y * q4x * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) - 3 * q2y * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) -
                           3 * q2z * q3x * q4y * Yket(2, 2, 0, 2, thetaq, wigner) + Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - 3 * q2x * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) + Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) +
                           Complex(0, 6) * q2x * q3x * q4z * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * q2y * q3x * q4z * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * q2x * q3y * q4z * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3y * q4z * Yket(2, 2, 0, 2, thetaq, wigner)) +
                      Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -1, 0, 2, thetaq, wigner) +
                           2 * Sqrt(3) * (q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 4) * q3y * q4x + 2 * q3x * (q4x - Complex(0, 1) * q4y) - q3z * q4z) + q2x * (2 * q3y * q4x - 4 * q3x * q4y - Complex(0, 2) * q3y * q4y + Complex(0, 1) * q3z * q4z)) *
                               Yket(2, 0, 0, 2, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) *
                               ((q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(2, 1, 0, 2, thetaq, wigner) -
                                2 * ((q2x - Complex(0, 1) * q2y) * q3z * q4z + q2z * (-2 * q3z * q4x + Complex(0, 2) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z)) * Yket(2, 2, 0, 2, thetaq, wigner))))) /
            (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 7 && idxket == 4)
    {
        mtx = -((Ftpe3 *
                 (Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 2) * (-(q3z * (q2x * q4x + q2y * q4y)) + q2z * (q3x * q4x + q3y * q4y)) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 1, 1, thetaq, wigner)) /
                Sqrt(3));
        return mtx;
    }
    else if (idxbra == 7 && idxket == 5)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) - q2y * q3z * (q4x + Complex(0, 1) * q4y) + q2x * q3z * (Complex(0, -1) * q4x + q4y) + q2y * q3x * q4z - q2x * q3y * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) +
                                                                                   Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * q4z + q2x * (q3y * q4x - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) - q2y * (q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                      Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                  2 * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + q2y * (-(q3z * q4x) + Complex(0, 1) * q3z * q4y + q3x * q4z) + q2x * (Complex(0, 1) * q3z * q4x + q3z * q4y - q3y * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                      Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                                                                                  2 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                  Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * q4z + q2x * (q3y * q4x - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) - q2y * (q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
            3.;
        return mtx;
    }
    else if (idxbra == 7 && idxket == 12)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) *
                    (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                     q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) *
                    (q2z * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, -3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y - Complex(0, 2) * q3y * q4y - Complex(0, 2) * q3z * q4z) +
                     q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              3.;
        return mtx;
    }
    else if (idxbra == 7 && idxket == 13)
    {
        mtx = (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, -2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                                                                                     Sqrt(2) *
                                                                                         (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                                                                                          q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                                                                                         Yket(1, 0, 2, 2, thetaq, wigner)) +
                        Sqrt(2) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((q2z * (Complex(0, -1) * q3z * q4x + q3z * q4y + Complex(0, 2) * q3x * q4z - 2 * q3y * q4z) + q2x * (-2 * q3y * q4x + Complex(0, 2) * q3y * q4y + q3x * (Complex(0, 3) * q4x + q4y) + Complex(0, 2) * q3z * q4z) -
                              q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                                 Yket(1, -1, 2, 2, thetaq, wigner) +
                             (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                              q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                                 Yket(1, 1, 2, 2, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) *
                                                                                        (q2z * (Complex(0, -1) * q3z * q4x + q3z * q4y + Complex(0, 2) * q3x * q4z - 2 * q3y * q4z) + q2x * (-2 * q3y * q4x + Complex(0, 2) * q3y * q4y + q3x * (Complex(0, 3) * q4x + q4y) + Complex(0, 2) * q3z * q4z) -
                                                                                         q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                                                                                        Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                    Complex(0, 2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              (3. * Sqrt(3));
        return mtx;
    }
    else if (idxbra == 7 && idxket == 14)
    {
        mtx = (Ftpe3 * (Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (3 * Sqrt(2) * (Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                             2 * (q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (-2 * q3y * q4x + 4 * q3x * q4y + Complex(0, 2) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (-2 * q3x * q4x - Complex(0, 4) * q3y * q4x + Complex(0, 2) * q3x * q4y + q3z * q4z)) *
                                 Yket(1, 0, 2, 2, thetaq, wigner) +
                             Complex(0, 1) * Sqrt(2) * (q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                        Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, -1) * Sqrt(2) * (q2z * q3x * (q4x + Complex(0, 3) * q4y) + q2x * q3z * (q4x + Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, -3) * q4x + q4y) + q2y * q3z * (Complex(0, -3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, -1, 2, 2, thetaq, wigner) +
                             (2 * q2z * (Complex(0, 2) * q3z * q4x + 2 * q3z * q4y - Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, -8) * q3y * q4x + 4 * q3x * (q4x + Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (4 * q3y * q4x - 8 * q3x * q4y + Complex(0, 4) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                                 Yket(1, 0, 2, 2, thetaq, wigner) +
                             3 * Sqrt(2) * (Complex(0, -1) * q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                        2 * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((2 * q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                                 Yket(1, -1, 2, 2, thetaq, wigner) +
                             3 * Sqrt(2) * (q2z * q3y * q4x + q2y * q3z * q4x - q2z * q3x * q4y - q2x * q3z * q4y) * Yket(1, 0, 2, 2, thetaq, wigner) -
                             (2 * q2z * (Complex(0, 2) * q3z * q4x + 2 * q3z * q4y - Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                                 Yket(1, 1, 2, 2, thetaq, wigner)))) /
              (6. * Sqrt(3));
        return mtx;
    }
    else if (idxbra == 7 && idxket == 15)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * Sqrt(2) * ((Complex(0, -1) * q2x + q2y) * q3z * q4z + q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z)) * Yket(2, -2, 2, 0, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * (q2z * q3x * (q4x + Complex(0, 3) * q4y) + q2x * q3z * (q4x + Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, -3) * q4x + q4y) + q2y * q3z * (Complex(0, -3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(2, -1, 2, 0, thetaq, wigner) +
                           4 * Sqrt(3) * q2y * q3x * q4x * Yket(2, 0, 2, 0, thetaq, wigner) + 4 * Sqrt(3) * q2x * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2y * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2z * q3z * q4x * Yket(2, 0, 2, 0, thetaq, wigner) -
                           8 * Sqrt(3) * q2x * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2y * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2x * q3y * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + 4 * Sqrt(3) * q2z * q3z * q4y * Yket(2, 0, 2, 0, thetaq, wigner) -
                           Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - 2 * Sqrt(3) * q2z * q3y * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - 2 * Sqrt(3) * q2y * q3z * q4z * Yket(2, 0, 2, 0, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 1, 2, 0, thetaq, wigner) -
                           3 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 1, 2, 0, thetaq, wigner) +
                           Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3y * q4z * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Yket(2, 1, 2, 0, thetaq, wigner)) +
                      2 * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 2, 0, thetaq, wigner) +
                           3 * (Complex(0, 2) * q2y * q3y * q4x + q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2x * q3y * (q4x - Complex(0, 1) * q4y) - 2 * q2x * q3x * q4y) * Yket(2, -1, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3x * q4x * Yket(2, 0, 2, 0, thetaq, wigner) -
                           Complex(0, 1) * Sqrt(6) * q2x * q3z * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * q4y * Yket(2, 0, 2, 0, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(6) * q2x * q3x * q4z * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3y * q4z * Yket(2, 0, 2, 0, thetaq, wigner) + 3 * q2y * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + 3 * q2x * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) -
                           Complex(0, 6) * q2y * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - 6 * q2x * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) -
                           Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 2, 0, thetaq, wigner) - 3 * q2z * q3y * q4x * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) - 3 * q2y * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) -
                           3 * q2z * q3x * q4y * Yket(2, 2, 2, 0, thetaq, wigner) + Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - 3 * q2x * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) + Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) +
                           Complex(0, 6) * q2x * q3x * q4z * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * q2y * q3x * q4z * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * q2x * q3y * q4z * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * q2y * q3y * q4z * Yket(2, 2, 2, 0, thetaq, wigner)) +
                      Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -1, 2, 0, thetaq, wigner) +
                           2 * Sqrt(3) * (q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 4) * q3y * q4x + 2 * q3x * (q4x - Complex(0, 1) * q4y) - q3z * q4z) + q2x * (2 * q3y * q4x - 4 * q3x * q4y - Complex(0, 2) * q3y * q4y + Complex(0, 1) * q3z * q4z)) *
                               Yket(2, 0, 2, 0, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) *
                               ((q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(2, 1, 2, 0, thetaq, wigner) -
                                2 * ((q2x - Complex(0, 1) * q2y) * q3z * q4z + q2z * (-2 * q3z * q4x + Complex(0, 2) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z)) * Yket(2, 2, 2, 0, thetaq, wigner))))) /
            (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 7 && idxket == 16)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * Sqrt(2) * ((Complex(0, -1) * q2x + q2y) * q3z * q4z + q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z)) * Yket(2, -2, 2, 2, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * (q2z * q3x * (q4x + Complex(0, 3) * q4y) + q2x * q3z * (q4x + Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, -3) * q4x + q4y) + q2y * q3z * (Complex(0, -3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(2, -1, 2, 2, thetaq, wigner) +
                           4 * Sqrt(3) * q2y * q3x * q4x * Yket(2, 0, 2, 2, thetaq, wigner) + 4 * Sqrt(3) * q2x * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2y * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2z * q3z * q4x * Yket(2, 0, 2, 2, thetaq, wigner) -
                           8 * Sqrt(3) * q2x * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2y * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2x * q3y * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + 4 * Sqrt(3) * q2z * q3z * q4y * Yket(2, 0, 2, 2, thetaq, wigner) -
                           Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - 2 * Sqrt(3) * q2z * q3y * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - 2 * Sqrt(3) * q2y * q3z * q4z * Yket(2, 0, 2, 2, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 1, 2, 2, thetaq, wigner) -
                           3 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 1, 2, 2, thetaq, wigner) +
                           Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3y * q4z * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Yket(2, 1, 2, 2, thetaq, wigner)) +
                      2 * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 2, 2, thetaq, wigner) +
                           3 * (Complex(0, 2) * q2y * q3y * q4x + q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2x * q3y * (q4x - Complex(0, 1) * q4y) - 2 * q2x * q3x * q4y) * Yket(2, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3x * q4x * Yket(2, 0, 2, 2, thetaq, wigner) -
                           Complex(0, 1) * Sqrt(6) * q2x * q3z * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * q4y * Yket(2, 0, 2, 2, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(6) * q2x * q3x * q4z * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3y * q4z * Yket(2, 0, 2, 2, thetaq, wigner) + 3 * q2y * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + 3 * q2x * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) -
                           Complex(0, 6) * q2y * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - 6 * q2x * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) -
                           Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 2, 2, thetaq, wigner) - 3 * q2z * q3y * q4x * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) - 3 * q2y * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) -
                           3 * q2z * q3x * q4y * Yket(2, 2, 2, 2, thetaq, wigner) + Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - 3 * q2x * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) + Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) +
                           Complex(0, 6) * q2x * q3x * q4z * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * q2y * q3x * q4z * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * q2x * q3y * q4z * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3y * q4z * Yket(2, 2, 2, 2, thetaq, wigner)) +
                      Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -1, 2, 2, thetaq, wigner) +
                           2 * Sqrt(3) * (q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 4) * q3y * q4x + 2 * q3x * (q4x - Complex(0, 1) * q4y) - q3z * q4z) + q2x * (2 * q3y * q4x - 4 * q3x * q4y - Complex(0, 2) * q3y * q4y + Complex(0, 1) * q3z * q4z)) *
                               Yket(2, 0, 2, 2, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) *
                               ((q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(2, 1, 2, 2, thetaq, wigner) -
                                2 * ((q2x - Complex(0, 1) * q2y) * q3z * q4z + q2z * (-2 * q3z * q4x + Complex(0, 2) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z)) * Yket(2, 2, 2, 2, thetaq, wigner))))) /
            (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 2)
    {
        mtx = -0.3333333333333333 *
              (Ftpe3 *
               (Sqrt(2) * (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                    Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (q2z * q3x * q4x + q2x * q3z * q4x + q2z * q3y * q4y + q2y * q3z * q4y - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) *
                    Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              Sqrt(2);
        return mtx;
    }
    else if (idxbra == 8 && idxket == 3)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * (q2z * q3z * (Complex(0, -1) * q4x + q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + (Complex(0, -1) * q2x + q2y) * q3z * q4z) * Yket(2, -2, 0, 2, thetaq, wigner) +
                           Complex(0, 6) * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y - 3 * q3z * q4z)) * Yket(2, -1, 0, 2, thetaq, wigner) - Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Yket(2, 0, 0, 2, thetaq, wigner) -
                           Sqrt(6) * q2y * q3x * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Yket(2, 0, 0, 2, thetaq, wigner) -
                           Sqrt(6) * q2x * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - 3 * Sqrt(6) * q2y * q3y * q4y * Yket(2, 0, 0, 2, thetaq, wigner) +
                           5 * Sqrt(6) * q2z * q3z * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Yket(2, 0, 0, 2, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3y * q4z * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Yket(2, 0, 0, 2, thetaq, wigner) +
                           5 * Sqrt(6) * q2y * q3z * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 12) * q2z * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - 12 * q2z * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 12) * q2x * q3z * q4x * Yket(2, 1, 0, 2, thetaq, wigner) -
                           12 * q2y * q3z * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - 12 * q2z * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 12) * q2z * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - 12 * q2x * q3z * q4y * Yket(2, 1, 0, 2, thetaq, wigner) +
                           Complex(0, 12) * q2y * q3z * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 12) * q2x * q3x * q4z * Yket(2, 1, 0, 2, thetaq, wigner) - 12 * q2y * q3x * q4z * Yket(2, 1, 0, 2, thetaq, wigner) - 12 * q2x * q3y * q4z * Yket(2, 1, 0, 2, thetaq, wigner) +
                           Complex(0, 12) * q2y * q3y * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 18) * q2x * q3x * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 18 * q2y * q3x * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 18 * q2x * q3y * q4x * Yket(2, 2, 0, 2, thetaq, wigner) -
                           Complex(0, 18) * q2y * q3y * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 18 * q2x * q3x * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 18) * q2y * q3x * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 18) * q2x * q3y * q4y * Yket(2, 2, 0, 2, thetaq, wigner) -
                           18 * q2y * q3y * q4y * Yket(2, 2, 0, 2, thetaq, wigner)) +
                      Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (Complex(0, 6) * Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 0, 2, thetaq, wigner) +
                           3 * Sqrt(2) *
                               (Complex(0, 3) * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2y * (q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + 3 * q3y * q4y - 3 * q3z * q4z) +
                                q2x * (q3y * q4x - Complex(0, 1) * q3y * q4y + q3x * (Complex(0, -3) * q4x + q4y) + Complex(0, 3) * q3z * q4z)) *
                               Yket(2, -1, 0, 2, thetaq, wigner) -
                           Complex(0, 8) * Sqrt(3) * q2z * q3x * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2x * q3z * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2z * q3y * q4y * Yket(2, 0, 0, 2, thetaq, wigner) -
                           Complex(0, 8) * Sqrt(3) * q2y * q3z * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2x * q3x * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2y * q3y * q4z * Yket(2, 0, 0, 2, thetaq, wigner) +
                           Complex(0, 12) * Sqrt(3) * q2z * q3z * q4z * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2x * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + 3 * Sqrt(2) * q2y * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(2) * q2y * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2z * q3z * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + 9 * Sqrt(2) * q2y * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - 9 * Sqrt(2) * q2z * q3z * q4y * Yket(2, 1, 0, 2, thetaq, wigner) -
                           Complex(0, 9) * Sqrt(2) * q2z * q3x * q4z * Yket(2, 1, 0, 2, thetaq, wigner) - 9 * Sqrt(2) * q2z * q3y * q4z * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2x * q3z * q4z * Yket(2, 1, 0, 2, thetaq, wigner) - 9 * Sqrt(2) * q2y * q3z * q4z * Yket(2, 1, 0, 2, thetaq, wigner) +
                           Complex(0, 6) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) +
                           6 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) +
                           Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3y * q4z * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Yket(2, 2, 0, 2, thetaq, wigner)) +
                      Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (18 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) * Yket(2, -2, 0, 2, thetaq, wigner) +
                           12 * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z)) * Yket(2, -1, 0, 2, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(6) * q2y * q3x * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) -
                           Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(6) * q2x * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) +
                           Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - 3 * Sqrt(6) * q2y * q3y * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3z * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Yket(2, 0, 0, 2, thetaq, wigner) +
                           5 * Sqrt(6) * q2z * q3y * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Yket(2, 0, 0, 2, thetaq, wigner) + 5 * Sqrt(6) * q2y * q3z * q4z * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) +
                           Complex(0, 6) * q2x * q3z * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3z * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3x * q4z * Yket(2, 1, 0, 2, thetaq, wigner) +
                           Complex(0, 6) * q2y * q3y * q4z * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 18) * q2z * q3z * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * q2z * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) +
                           Complex(0, 6) * q2z * q3x * q4z * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * q2z * q3y * q4z * Yket(2, 2, 0, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3z * q4z * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * q2y * q3z * q4z * Yket(2, 2, 0, 2, thetaq, wigner)))) /
            (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 4)
    {
        mtx = -((Ftpe3 *
                 (Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 2) * (-(q3z * (q2x * q4x + q2y * q4y)) + q2z * (q3x * q4x + q3y * q4y)) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 1, 1, thetaq, wigner)) /
                Sqrt(6));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 5)
    {
        mtx = (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((q2z * (q3x - Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) + q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, 1) * q3z * q4x - q3z * q4y - 2 * q3y * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (q2x * q3y * q4x - Complex(0, 1) * q2x * q3y * q4y - Complex(0, 2) * q2z * q3x * q4z - 2 * q2z * q3y * q4z + Complex(0, 2) * q2x * q3z * q4z + q2y * (-(q3x * q4x) + Complex(0, 1) * q3x * q4y + 2 * q3z * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                             3 * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * (Complex(0, 1) * q4x + q4y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (3 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(1, -1, 1, 1, thetaq, wigner) +
                                                                                    Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) - q2x * q3y * (q4x + Complex(0, 1) * q4y) + 2 * q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 2) * q2x * q3z * q4z - 2 * q2y * q3z * q4z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                    (q2z * (q3x + Complex(0, 1) * q3y) * (Complex(0, 1) * q4x + q4y) + q2y * (q3z * q4x - Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x - q3z * q4y - 2 * q3y * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Sqrt(2) * (-2 * q2y * q3x * (q4x + Complex(0, 1) * q4y) + 2 * q2x * q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z - Complex(0, 1) * q2x * q3z * q4z + q2y * q3z * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             2 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y - 2 * q2y * q3x * q4z + 2 * q2x * q3y * q4z) * Yket(1, 0, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (2 * q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z - q2y * q3z * q4z + q2x * (-2 * q3y * q4x + Complex(0, 2) * q3y * q4y - Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 12)
    {
        mtx = -0.3333333333333333 *
              (Ftpe3 *
               (Sqrt(2) * (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                    Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (q2z * q3x * q4x + q2x * q3z * q4x + q2z * q3y * q4y + q2y * q3z * q4y - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) *
                    Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              Sqrt(2);
        return mtx;
    }
    else if (idxbra == 8 && idxket == 13)
    {
        mtx =
            (Ftpe3 * (Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                           Sqrt(2) * (2 * q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                               Yket(1, 0, 2, 2, thetaq, wigner) +
                           Complex(0, 1) * (q2z * q3x * (q4x + Complex(0, 3) * q4y) + q2x * q3z * (q4x + Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, -3) * q4x + q4y) + q2y * q3z * (Complex(0, -3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                      Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (Complex(0, -1) * (q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, -1, 2, 2, thetaq, wigner) -
                           Sqrt(2) * (2 * q2z * (Complex(0, 2) * q3z * q4x + 2 * q3z * q4y - Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                               Yket(1, 0, 2, 2, thetaq, wigner) +
                           3 * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                      Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (Sqrt(2) * (q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (-2 * q3y * q4x + 4 * q3x * q4y + Complex(0, 2) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (-2 * q3x * q4x - Complex(0, 4) * q3y * q4x + Complex(0, 2) * q3x * q4y + q3z * q4z)) *
                               Yket(1, -1, 2, 2, thetaq, wigner) -
                           6 * (q2z * q3y * q4x + q2y * q3z * q4x - q2z * q3x * q4y - q2x * q3z * q4y) * Yket(1, 0, 2, 2, thetaq, wigner) +
                           Sqrt(2) * (q2z * (Complex(0, 2) * q3z * q4x + 2 * q3z * q4y - Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, -4) * q3y * q4x + 2 * q3x * (q4x + Complex(0, 1) * q4y) - q3z * q4z) + q2x * (2 * q3y * q4x - 4 * q3x * q4y + Complex(0, 2) * q3y * q4y - Complex(0, 1) * q3z * q4z)) *
                               Yket(1, 1, 2, 2, thetaq, wigner)))) /
            (3. * Sqrt(6));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 14)
    {
        mtx = (Complex(0, 0.16666666666666666) * Ftpe3 *
               (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y + 3 * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) -
                                                                             Sqrt(2) *
                                                                                 (q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y + q3z * q4z) -
                                                                                  Complex(0, 1) * (q2z * (Complex(0, 1) * q3z * q4x + q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2y * (q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + 3 * q3y * q4y + q3z * q4z))) *
                                                                                 Yket(1, 0, 2, 2, thetaq, wigner)) -
                Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) *
                                                                                (q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y + Complex(0, 1) * q3z * q4z) +
                                                                                 q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y + q3z * q4z)) *
                                                                                Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                            2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y + 3 * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)) -
                Sqrt(2) * Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                    ((q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y + Complex(0, 1) * q3z * q4z) +
                      q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y + q3z * q4z)) *
                         Yket(1, -1, 2, 2, thetaq, wigner) +
                     (q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y + q3z * q4z) -
                      Complex(0, 1) * (q2z * (Complex(0, 1) * q3z * q4x + q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2y * (q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + 3 * q3y * q4y + q3z * q4z))) *
                         Yket(1, 1, 2, 2, thetaq, wigner)))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 8 && idxket == 15)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * (q2z * q3z * (Complex(0, -1) * q4x + q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + (Complex(0, -1) * q2x + q2y) * q3z * q4z) * Yket(2, -2, 2, 0, thetaq, wigner) +
                           Complex(0, 6) * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y - 3 * q3z * q4z)) * Yket(2, -1, 2, 0, thetaq, wigner) - Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Yket(2, 0, 2, 0, thetaq, wigner) -
                           Sqrt(6) * q2y * q3x * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Yket(2, 0, 2, 0, thetaq, wigner) -
                           Sqrt(6) * q2x * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - 3 * Sqrt(6) * q2y * q3y * q4y * Yket(2, 0, 2, 0, thetaq, wigner) +
                           5 * Sqrt(6) * q2z * q3z * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Yket(2, 0, 2, 0, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3y * q4z * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Yket(2, 0, 2, 0, thetaq, wigner) +
                           5 * Sqrt(6) * q2y * q3z * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 12) * q2z * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - 12 * q2z * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 12) * q2x * q3z * q4x * Yket(2, 1, 2, 0, thetaq, wigner) -
                           12 * q2y * q3z * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - 12 * q2z * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 12) * q2z * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - 12 * q2x * q3z * q4y * Yket(2, 1, 2, 0, thetaq, wigner) +
                           Complex(0, 12) * q2y * q3z * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 12) * q2x * q3x * q4z * Yket(2, 1, 2, 0, thetaq, wigner) - 12 * q2y * q3x * q4z * Yket(2, 1, 2, 0, thetaq, wigner) - 12 * q2x * q3y * q4z * Yket(2, 1, 2, 0, thetaq, wigner) +
                           Complex(0, 12) * q2y * q3y * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 18) * q2x * q3x * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 18 * q2y * q3x * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 18 * q2x * q3y * q4x * Yket(2, 2, 2, 0, thetaq, wigner) -
                           Complex(0, 18) * q2y * q3y * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 18 * q2x * q3x * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 18) * q2y * q3x * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 18) * q2x * q3y * q4y * Yket(2, 2, 2, 0, thetaq, wigner) -
                           18 * q2y * q3y * q4y * Yket(2, 2, 2, 0, thetaq, wigner)) +
                      Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (Complex(0, 6) * Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 2, 0, thetaq, wigner) +
                           3 * Sqrt(2) *
                               (Complex(0, 3) * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2y * (q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + 3 * q3y * q4y - 3 * q3z * q4z) +
                                q2x * (q3y * q4x - Complex(0, 1) * q3y * q4y + q3x * (Complex(0, -3) * q4x + q4y) + Complex(0, 3) * q3z * q4z)) *
                               Yket(2, -1, 2, 0, thetaq, wigner) -
                           Complex(0, 8) * Sqrt(3) * q2z * q3x * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2x * q3z * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2z * q3y * q4y * Yket(2, 0, 2, 0, thetaq, wigner) -
                           Complex(0, 8) * Sqrt(3) * q2y * q3z * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2x * q3x * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2y * q3y * q4z * Yket(2, 0, 2, 0, thetaq, wigner) +
                           Complex(0, 12) * Sqrt(3) * q2z * q3z * q4z * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2x * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + 3 * Sqrt(2) * q2y * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(2) * q2y * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2z * q3z * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + 9 * Sqrt(2) * q2y * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - 9 * Sqrt(2) * q2z * q3z * q4y * Yket(2, 1, 2, 0, thetaq, wigner) -
                           Complex(0, 9) * Sqrt(2) * q2z * q3x * q4z * Yket(2, 1, 2, 0, thetaq, wigner) - 9 * Sqrt(2) * q2z * q3y * q4z * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2x * q3z * q4z * Yket(2, 1, 2, 0, thetaq, wigner) - 9 * Sqrt(2) * q2y * q3z * q4z * Yket(2, 1, 2, 0, thetaq, wigner) +
                           Complex(0, 6) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) +
                           6 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) +
                           Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3y * q4z * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Yket(2, 2, 2, 0, thetaq, wigner)) +
                      Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (18 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) * Yket(2, -2, 2, 0, thetaq, wigner) +
                           12 * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z)) * Yket(2, -1, 2, 0, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(6) * q2y * q3x * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) -
                           Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(6) * q2x * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) +
                           Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - 3 * Sqrt(6) * q2y * q3y * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3z * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Yket(2, 0, 2, 0, thetaq, wigner) +
                           5 * Sqrt(6) * q2z * q3y * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Yket(2, 0, 2, 0, thetaq, wigner) + 5 * Sqrt(6) * q2y * q3z * q4z * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 6) * q2z * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) +
                           Complex(0, 6) * q2x * q3z * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2z * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2y * q3z * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2x * q3x * q4z * Yket(2, 1, 2, 0, thetaq, wigner) +
                           Complex(0, 6) * q2y * q3y * q4z * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 18) * q2z * q3z * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 6) * q2z * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * q2z * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) +
                           Complex(0, 6) * q2z * q3x * q4z * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * q2z * q3y * q4z * Yket(2, 2, 2, 0, thetaq, wigner) + Complex(0, 6) * q2x * q3z * q4z * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * q2y * q3z * q4z * Yket(2, 2, 2, 0, thetaq, wigner)))) /
            (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 8 && idxket == 16)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * (q2z * q3z * (Complex(0, -1) * q4x + q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + (Complex(0, -1) * q2x + q2y) * q3z * q4z) * Yket(2, -2, 2, 2, thetaq, wigner) +
                           Complex(0, 6) * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y - 3 * q3z * q4z)) * Yket(2, -1, 2, 2, thetaq, wigner) - Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Yket(2, 0, 2, 2, thetaq, wigner) -
                           Sqrt(6) * q2y * q3x * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Yket(2, 0, 2, 2, thetaq, wigner) -
                           Sqrt(6) * q2x * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - 3 * Sqrt(6) * q2y * q3y * q4y * Yket(2, 0, 2, 2, thetaq, wigner) +
                           5 * Sqrt(6) * q2z * q3z * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Yket(2, 0, 2, 2, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3y * q4z * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Yket(2, 0, 2, 2, thetaq, wigner) +
                           5 * Sqrt(6) * q2y * q3z * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 12) * q2z * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - 12 * q2z * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 12) * q2x * q3z * q4x * Yket(2, 1, 2, 2, thetaq, wigner) -
                           12 * q2y * q3z * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - 12 * q2z * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 12) * q2z * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - 12 * q2x * q3z * q4y * Yket(2, 1, 2, 2, thetaq, wigner) +
                           Complex(0, 12) * q2y * q3z * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 12) * q2x * q3x * q4z * Yket(2, 1, 2, 2, thetaq, wigner) - 12 * q2y * q3x * q4z * Yket(2, 1, 2, 2, thetaq, wigner) - 12 * q2x * q3y * q4z * Yket(2, 1, 2, 2, thetaq, wigner) +
                           Complex(0, 12) * q2y * q3y * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 18) * q2x * q3x * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 18 * q2y * q3x * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 18 * q2x * q3y * q4x * Yket(2, 2, 2, 2, thetaq, wigner) -
                           Complex(0, 18) * q2y * q3y * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 18 * q2x * q3x * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 18) * q2y * q3x * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 18) * q2x * q3y * q4y * Yket(2, 2, 2, 2, thetaq, wigner) -
                           18 * q2y * q3y * q4y * Yket(2, 2, 2, 2, thetaq, wigner)) +
                      Ybra(1, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (Complex(0, 6) * Sqrt(2) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 2, 2, thetaq, wigner) +
                           3 * Sqrt(2) *
                               (Complex(0, 3) * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2y * (q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + 3 * q3y * q4y - 3 * q3z * q4z) +
                                q2x * (q3y * q4x - Complex(0, 1) * q3y * q4y + q3x * (Complex(0, -3) * q4x + q4y) + Complex(0, 3) * q3z * q4z)) *
                               Yket(2, -1, 2, 2, thetaq, wigner) -
                           Complex(0, 8) * Sqrt(3) * q2z * q3x * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2x * q3z * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2z * q3y * q4y * Yket(2, 0, 2, 2, thetaq, wigner) -
                           Complex(0, 8) * Sqrt(3) * q2y * q3z * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2x * q3x * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2y * q3y * q4z * Yket(2, 0, 2, 2, thetaq, wigner) +
                           Complex(0, 12) * Sqrt(3) * q2z * q3z * q4z * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2x * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + 3 * Sqrt(2) * q2y * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(2) * q2y * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2z * q3z * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + 9 * Sqrt(2) * q2y * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - 9 * Sqrt(2) * q2z * q3z * q4y * Yket(2, 1, 2, 2, thetaq, wigner) -
                           Complex(0, 9) * Sqrt(2) * q2z * q3x * q4z * Yket(2, 1, 2, 2, thetaq, wigner) - 9 * Sqrt(2) * q2z * q3y * q4z * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2x * q3z * q4z * Yket(2, 1, 2, 2, thetaq, wigner) - 9 * Sqrt(2) * q2y * q3z * q4z * Yket(2, 1, 2, 2, thetaq, wigner) +
                           Complex(0, 6) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) +
                           6 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) +
                           Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3y * q4z * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Yket(2, 2, 2, 2, thetaq, wigner)) +
                      Ybra(1, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (18 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) * Yket(2, -2, 2, 2, thetaq, wigner) +
                           12 * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z)) * Yket(2, -1, 2, 2, thetaq, wigner) +
                           Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2y * q3x * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) -
                           Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) +
                           Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - 3 * Sqrt(6) * q2y * q3y * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3z * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Yket(2, 0, 2, 2, thetaq, wigner) +
                           5 * Sqrt(6) * q2z * q3y * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Yket(2, 0, 2, 2, thetaq, wigner) + 5 * Sqrt(6) * q2y * q3z * q4z * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) +
                           Complex(0, 6) * q2x * q3z * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3z * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3x * q4z * Yket(2, 1, 2, 2, thetaq, wigner) +
                           Complex(0, 6) * q2y * q3y * q4z * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 18) * q2z * q3z * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * q2z * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) +
                           Complex(0, 6) * q2z * q3x * q4z * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * q2z * q3y * q4z * Yket(2, 2, 2, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3z * q4z * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * q2y * q3z * q4z * Yket(2, 2, 2, 2, thetaq, wigner)))) /
            (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 9 && idxket == 2)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (Complex(0, -1) * q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (-(q3y * (q4x + Complex(0, 1) * q4y)) + 2 * q3x * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 2) * q3y * q4x - q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                    Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 2) * Sqrt(2) * q2z * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2z * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 2) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              Sqrt(10);
        return mtx;
    }
    else if (idxbra == 9 && idxket == 3)
    {
        mtx = (Complex(0, 0.1) * Ftpe3 *
               (12 * q2z * q3z * q4z * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 0, 2, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                3 * q2z * q3z * q4z * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 0, 2, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 0, 2, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) -
                3 * q2z * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) -
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 0, 2, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner) -
                12 * q2z * q3z * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 0, 2, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 9 && idxket == 4)
    {
        mtx = (Ftpe3 *
               (3 * Sqrt(2) * (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * (q2x * q3y * q4x - Complex(0, 1) * q2x * q3y * q4y - Complex(0, 1) * q2z * q3x * q4z - q2z * q3y * q4z + Complex(0, 1) * q2x * q3z * q4z + q2y * (-(q3x * q4x) + Complex(0, 1) * q3x * q4y + q3z * q4z)) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - 4 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 4 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 3) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 9 && idxket == 5)
    {
        mtx = (Ftpe3 * (-2 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(2) * q2z * q3y * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Sqrt(2) * q2y * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Sqrt(3) * q2z * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(3) * q2y * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Sqrt(3) * q2z * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * q2z * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * q2x * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * q2z * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * q2y * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2y * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2z * q3y * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Sqrt(6) * q2z * q3x * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2z * q3y * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2x * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Sqrt(3) * (q2z * (Complex(0, 1) * q3x + q3y) - Complex(0, 1) * q2x * q3z - q2y * q3z) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * q4z * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(2) * (q4x - Complex(0, 1) * q4y) * Yket(1, 0, 1, 1, thetaq, wigner)) -
                        2 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 2 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(2) * q2z * q3y * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Sqrt(2) * q2y * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Sqrt(3) * q2z * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(3) * q2y * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Sqrt(3) * q2z * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Sqrt(3) * q2x * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        2 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 2 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 2 * Sqrt(3) * q2z * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 2 * Sqrt(3) * q2y * q3z * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Sqrt(3) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((Complex(0, 1) * q2x * q3z * q4x + q2y * q3z * q4x - q2x * q3z * q4y + Complex(0, 1) * q2y * q3z * q4y + q2z * (q3x - Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) + 2 * q2y * q3x * q4z - 2 * q2x * q3y * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) -
                             (q4x - Complex(0, 1) * q4y) * (Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 1, 1, thetaq, wigner) + (Complex(0, 1) * q2z * q3x + q2z * q3y - Complex(0, 1) * q2x * q3z - q2y * q3z) * Yket(1, 1, 1, 1, thetaq, wigner))))) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 9 && idxket == 12)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (Complex(0, -1) * q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (-(q3y * (q4x + Complex(0, 1) * q4y)) + 2 * q3x * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 2) * q3y * q4x - q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                    Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 2) * Sqrt(2) * q2z * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2z * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 2) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              Sqrt(10);
        return mtx;
    }
    else if (idxbra == 9 && idxket == 13)
    {
        mtx = -0.3333333333333333 *
              (Ftpe3 * (2 * Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 4) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        4 * Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 2 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 3) * q2z * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - 3 * q2z * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 3) * q2x * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - 3 * q2y * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        3 * q2z * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        3 * q2x * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 3) * q2y * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 6 * q2y * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        6 * q2x * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 4) * Sqrt(3) * q2x * q3x * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2y * q3y * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 6 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        6 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        3 * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * ((Complex(0, 1) * q2x + q2y) * q3z * q4z + q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                             Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        2 * Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 4) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        4 * Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 2 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 3) * q2z * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 9 * q2z * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 3) * q2x * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 9 * q2y * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        9 * q2z * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2z * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        9 * q2x * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2y * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 12) * q2z * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 12 * q2z * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 6 * q2z * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 6 * q2y * q3z * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        3 * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, 1) * (q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, -1, 2, 2, thetaq, wigner) +
                             Sqrt(2) * (Complex(0, -2) * q2y * q3y * q4x + q2y * q3x * (q4x + Complex(0, 1) * q4y) + q2x * q3y * (q4x + Complex(0, 1) * q4y) - 2 * q2x * q3x * q4y) * Yket(1, 0, 2, 2, thetaq, wigner) +
                             (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 9 && idxket == 14)
    {
        mtx = (Ftpe3 * (Complex(0, -3) * Sqrt(6) * q2x * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        3 * Sqrt(6) * q2y * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - 5 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - 5 * Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - 5 * Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 12) * q2z * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 12 * q2z * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 12) * q2x * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 12 * q2y * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        12 * q2z * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 12) * q2z * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        12 * q2x * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 12) * q2y * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 12) * q2x * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 12 * q2y * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        12 * q2x * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 12) * q2y * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 18) * q2x * q3x * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - 18 * q2y * q3x * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        18 * q2x * q3y * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 18) * q2y * q3y * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        18 * q2x * q3x * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 18) * q2y * q3x * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 18) * q2x * q3y * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + 18 * q2y * q3y * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 8) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 8) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 8) * Sqrt(3) * q2x * q3x * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 8) * Sqrt(3) * q2y * q3y * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 12) * Sqrt(3) * q2z * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2x * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2z * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        3 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 9 * Sqrt(2) * q2y * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        9 * Sqrt(2) * q2z * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        9 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        9 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        6 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        6 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 6 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 6) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 6 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 6) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        6 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - 6 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 3 * Sqrt(6) * q2y * q3y * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        5 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        5 * Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        5 * Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3x * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 6) * q2x * q3z * q4x * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 6) * q2y * q3z * q4y * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 18) * q2z * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 6) * q2z * q3z * q4x * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 6 * q2z * q3z * q4y * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 6) * q2z * q3x * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 6 * q2z * q3y * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 6) * q2x * q3z * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 6 * q2y * q3z * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        6 * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, -1) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Yket(1, -1, 2, 2, thetaq, wigner) +
                             Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner) -
                             Complex(0, 3) * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                        3 * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, 2) * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y - 3 * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) -
                             Sqrt(2) *
                                 (Complex(0, -3) * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2y * (q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + 3 * q3y * q4y - 3 * q3z * q4z) +
                                  q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y + q3x * (Complex(0, 3) * q4x + q4y) - Complex(0, 3) * q3z * q4z)) *
                                 Yket(1, 0, 2, 2, thetaq, wigner) -
                             Complex(0, 4) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 9 && idxket == 15)
    {
        mtx = (Complex(0, 0.1) * Ftpe3 *
               (12 * q2z * q3z * q4z * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 0, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                3 * q2z * q3z * q4z * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 0, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 0, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) -
                3 * q2z * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) -
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 0, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner) -
                12 * q2z * q3z * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 0, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 9 && idxket == 16)
    {
        mtx = (Complex(0, 0.1) * Ftpe3 *
               (12 * q2z * q3z * q4z * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 2, 2, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                3 * q2z * q3z * q4z * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 2, 2, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 2, 2, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) -
                3 * q2z * q3z * q4z * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) -
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 2, 2, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner) -
                12 * q2z * q3z * q4z * Ybra(2, 2, 1, 1, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 2, 2, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 10 && idxket == 2)
    {
        mtx = 2 * Ftpe3 * (-(q2z * q3y * q4x) + q2y * q3z * q4x + q2z * q3x * q4y - q2x * q3z * q4y - q2y * q3x * q4z + q2x * q3y * q4z) * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 10 && idxket == 3)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               (3 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 0, 2, thetaq, wigner) +
                3 * (q2y * q3x * (q4x + Complex(0, 1) * q4y) - q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(2, -1, 0, 2, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) +
                Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 0, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) - Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + 2 * Sqrt(6) * q2y * q3x * q4z * Yket(2, 0, 0, 2, thetaq, wigner) -
                2 * Sqrt(6) * q2x * q3y * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - 3 * q2y * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + 3 * q2x * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) -
                Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 3) * q2z * q3x * q4z * Yket(2, 1, 0, 2, thetaq, wigner) - 3 * q2z * q3y * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * q2x * q3z * q4z * Yket(2, 1, 0, 2, thetaq, wigner) +
                3 * q2y * q3z * q4z * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 3 * q2z * q3y * q4x * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) -
                3 * q2y * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + 3 * q2z * q3x * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 0, 2, thetaq, wigner) - 3 * q2x * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) +
                Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner))) /
              Sqrt(15);
        return mtx;
    }
    else if (idxbra == 10 && idxket == 5)
    {
        mtx = Complex(0, 1) * Ftpe3 * (q2x * q3x + q2y * q3y + q2z * q3z) * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
              (Sqrt(2) * (q4x + Complex(0, 1) * q4y) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * q4z * Yket(1, 0, 1, 1, thetaq, wigner) - Sqrt(2) * (q4x - Complex(0, 1) * q4y) * Yket(1, 1, 1, 1, thetaq, wigner));
        return mtx;
    }
    else if (idxbra == 10 && idxket == 12)
    {
        mtx = 2 * Ftpe3 * (-(q2z * q3y * q4x) + q2y * q3z * q4x + q2z * q3x * q4y - q2x * q3z * q4y - q2y * q3x * q4z + q2x * q3y * q4z) * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 10 && idxket == 13)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                Complex(0, 2) * (q3z * (q2x * q4x + q2y * q4y) - q2z * (q3x * q4x + q3y * q4y)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 10 && idxket == 14)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                Complex(0, 1) * Sqrt(2) * (q3z * (q2x * q4x + q2y * q4y) - q2z * (q3x * q4x + q3y * q4y)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 10 && idxket == 15)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               (3 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 2, 0, thetaq, wigner) +
                3 * (q2y * q3x * (q4x + Complex(0, 1) * q4y) - q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(2, -1, 2, 0, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) +
                Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 2, 0, thetaq, wigner) + Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) - Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + 2 * Sqrt(6) * q2y * q3x * q4z * Yket(2, 0, 2, 0, thetaq, wigner) -
                2 * Sqrt(6) * q2x * q3y * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - 3 * q2y * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + 3 * q2x * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) -
                Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 3) * q2z * q3x * q4z * Yket(2, 1, 2, 0, thetaq, wigner) - 3 * q2z * q3y * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * q2x * q3z * q4z * Yket(2, 1, 2, 0, thetaq, wigner) +
                3 * q2y * q3z * q4z * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 3 * q2z * q3y * q4x * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) -
                3 * q2y * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + 3 * q2z * q3x * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 2, 0, thetaq, wigner) - 3 * q2x * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) +
                Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner))) /
              Sqrt(15);
        return mtx;
    }
    else if (idxbra == 10 && idxket == 16)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               (3 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 2, 2, thetaq, wigner) +
                3 * (q2y * q3x * (q4x + Complex(0, 1) * q4y) - q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(2, -1, 2, 2, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) +
                Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(6) * q2y * q3x * q4z * Yket(2, 0, 2, 2, thetaq, wigner) -
                2 * Sqrt(6) * q2x * q3y * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - 3 * q2y * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + 3 * q2x * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) -
                Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 3) * q2z * q3x * q4z * Yket(2, 1, 2, 2, thetaq, wigner) - 3 * q2z * q3y * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2x * q3z * q4z * Yket(2, 1, 2, 2, thetaq, wigner) +
                3 * q2y * q3z * q4z * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 3 * q2z * q3y * q4x * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) -
                3 * q2y * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + 3 * q2z * q3x * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 2, 2, thetaq, wigner) - 3 * q2x * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) +
                Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner))) /
              Sqrt(15);
        return mtx;
    }
    else if (idxbra == 11 && idxket == 2)
    {
        mtx = -((Ftpe3 *
                 (Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 2) * (-(q3z * (q2x * q4x + q2y * q4y)) + q2z * (q3x * q4x + q3y * q4y)) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 0, 0, thetaq, wigner)) /
                Sqrt(3));
        return mtx;
    }
    else if (idxbra == 11 && idxket == 3)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * Sqrt(2) * (Complex(0, 1) * q2z * q3x - q2z * q3y - Complex(0, 1) * q2x * q3z + q2y * q3z) * q4z * Yket(2, -2, 0, 2, thetaq, wigner) +
                           3 * Sqrt(2) * (Complex(0, 1) * q2x * q3z * q4x - q2y * q3z * q4x + q2z * (Complex(0, -1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + q2x * q3z * q4y + Complex(0, 1) * q2y * q3z * q4y - 2 * q2y * q3x * q4z + 2 * q2x * q3y * q4z) * Yket(2, -1, 0, 2, thetaq, wigner) +
                           4 * Sqrt(3) * q2y * q3x * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - 4 * Sqrt(3) * q2x * q3y * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2y * q3x * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2x * q3y * q4y * Yket(2, 0, 0, 2, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Yket(2, 0, 0, 2, thetaq, wigner) + 2 * Sqrt(3) * q2z * q3y * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Yket(2, 0, 0, 2, thetaq, wigner) - 2 * Sqrt(3) * q2y * q3z * q4z * Yket(2, 0, 0, 2, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 1, 0, 2, thetaq, wigner) + 3 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 1, 0, 2, thetaq, wigner) -
                           3 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 1, 0, 2, thetaq, wigner)) +
                      Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 0, 2, thetaq, wigner) + 6 * (q2y * q3x - q2x * q3y) * (q4x + Complex(0, 1) * q4y) * Yket(2, -1, 0, 2, thetaq, wigner) -
                           Complex(0, 2) * Sqrt(6) * q2z * q3x * q4x * Yket(2, 0, 0, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3z * q4x * Yket(2, 0, 0, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2z * q3y * q4y * Yket(2, 0, 0, 2, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(6) * q2y * q3z * q4y * Yket(2, 0, 0, 2, thetaq, wigner) + 6 * q2y * q3x * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - 6 * q2x * q3y * q4x * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3x * q4y * Yket(2, 1, 0, 2, thetaq, wigner) +
                           Complex(0, 6) * q2x * q3y * q4y * Yket(2, 1, 0, 2, thetaq, wigner) - Complex(0, 6) * q2z * q3x * q4x * Yket(2, 2, 0, 2, thetaq, wigner) - 6 * q2z * q3y * q4x * Yket(2, 2, 0, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) +
                           6 * q2y * q3z * q4x * Yket(2, 2, 0, 2, thetaq, wigner) - 6 * q2z * q3x * q4y * Yket(2, 2, 0, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * q4y * Yket(2, 2, 0, 2, thetaq, wigner) + 6 * q2x * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner) -
                           Complex(0, 6) * q2y * q3z * q4y * Yket(2, 2, 0, 2, thetaq, wigner)) +
                      Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * Sqrt(2) * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -1, 0, 2, thetaq, wigner) +
                           2 * Sqrt(3) * (2 * q2y * q3x * (q4x + Complex(0, 1) * q4y) - 2 * q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(2, 0, 0, 2, thetaq, wigner) +
                           3 * Sqrt(2) *
                               ((Complex(0, 1) * q2x * q3z * q4x + q2y * q3z * q4x - q2x * q3z * q4y + Complex(0, 1) * q2y * q3z * q4y + q2z * (q3x - Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) + 2 * q2y * q3x * q4z - 2 * q2x * q3y * q4z) * Yket(2, 1, 0, 2, thetaq, wigner) +
                                2 * (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * q4z * Yket(2, 2, 0, 2, thetaq, wigner))))) /
            (6. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 11 && idxket == 4)
    {
        mtx = Ftpe3 * (q2x * q3x + q2y * q3y + q2z * q3z) *
              (Sqrt(2) * (Complex(0, 1) * q4x + q4y) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 2) * q4z * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * (Complex(0, -1) * q4x + q4y) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
              Yket(0, 0, 1, 1, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 11 && idxket == 5)
    {
        mtx = (Ftpe3 * (q2x * q3x + q2y * q3y + q2z * q3z) *
               (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, -2) * q4z * Yket(1, -1, 1, 1, thetaq, wigner) + Sqrt(2) * (Complex(0, 1) * q4x + q4y) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                Sqrt(2) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, 1) * (q4x + Complex(0, 1) * q4y) * Yket(1, -1, 1, 1, thetaq, wigner) + (Complex(0, 1) * q4x + q4y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                Complex(0, 1) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q4x + Complex(0, 1) * q4y) * Yket(1, 0, 1, 1, thetaq, wigner) + 2 * q4z * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 11 && idxket == 12)
    {
        mtx = -((Ftpe3 *
                 (Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Complex(0, 2) * (-(q3z * (q2x * q4x + q2y * q4y)) + q2z * (q3x * q4x + q3y * q4y)) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                  Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
                 Yket(0, 0, 2, 2, thetaq, wigner)) /
                Sqrt(3));
        return mtx;
    }
    else if (idxbra == 11 && idxket == 13)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2z * (q3x + Complex(0, 1) * q3y) * (Complex(0, 1) * q4x + q4y) + q2y * (q3z * q4x - Complex(0, 1) * q3z * q4y - q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x - q3z * q4y + q3y * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                                                                                   Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * q4z + q2x * (q3y * q4x - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) - q2y * (q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                      Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                  2 * (q2z * (q3x - Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) + q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y - q3x * q4z) + q2x * (Complex(0, 1) * q3z * q4x - q3z * q4y + q3y * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                      Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) -
                                                                                  2 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                  Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * q4z + q2x * (q3y * q4x - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) - q2y * (q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
            3.;
        return mtx;
    }
    else if (idxbra == 11 && idxket == 14)
    {
        mtx = (Ftpe3 *
               (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) - q2y * (q3z * q4x - Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, 1) * q3z * q4x + q3z * q4y + 2 * q3y * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                                                                             2 * (2 * q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z - q2y * q3z * q4z + q2x * (-2 * q3y * q4x + Complex(0, 2) * q3y * q4y - Complex(0, 1) * q3z * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                             3 * Sqrt(2) * (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * (q4x - Complex(0, 1) * q4y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                    (3 * Sqrt(2) * (Complex(0, 1) * q2z * q3x - q2z * q3y - Complex(0, 1) * q2x * q3z + q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(1, -1, 2, 2, thetaq, wigner) +
                     (4 * q2x * q3y * q4x - 4 * q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 4) * q2x * q3y * q4y + Complex(0, 2) * q2z * q3x * q4z - 2 * q2z * q3y * q4z - Complex(0, 2) * q2x * q3z * q4z + 2 * q2y * q3z * q4z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                     Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) - q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x + q3z * q4y + 2 * q3y * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                2 * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                    ((q2y * q3x * (q4x + Complex(0, 1) * q4y) - q2x * q3y * (q4x + Complex(0, 1) * q4y) + 2 * q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 2) * q2x * q3z * q4z - 2 * q2y * q3z * q4z) * Yket(1, -1, 2, 2, thetaq, wigner) +
                     Sqrt(2) * (-(q2z * q3y * q4x) + q2y * q3z * q4x + q2z * q3x * q4y - q2x * q3z * q4y + 2 * q2y * q3x * q4z - 2 * q2x * q3y * q4z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                     (q2x * q3y * q4x - Complex(0, 1) * q2x * q3y * q4y - Complex(0, 2) * q2z * q3x * q4z - 2 * q2z * q3y * q4z + Complex(0, 2) * q2x * q3z * q4z + q2y * (-(q3x * q4x) + Complex(0, 1) * q3x * q4y + 2 * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              6.;
        return mtx;
    }
    else if (idxbra == 11 && idxket == 15)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * Sqrt(2) * (Complex(0, 1) * q2z * q3x - q2z * q3y - Complex(0, 1) * q2x * q3z + q2y * q3z) * q4z * Yket(2, -2, 2, 0, thetaq, wigner) +
                           3 * Sqrt(2) * (Complex(0, 1) * q2x * q3z * q4x - q2y * q3z * q4x + q2z * (Complex(0, -1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + q2x * q3z * q4y + Complex(0, 1) * q2y * q3z * q4y - 2 * q2y * q3x * q4z + 2 * q2x * q3y * q4z) * Yket(2, -1, 2, 0, thetaq, wigner) +
                           4 * Sqrt(3) * q2y * q3x * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - 4 * Sqrt(3) * q2x * q3y * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2y * q3x * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2x * q3y * q4y * Yket(2, 0, 2, 0, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Yket(2, 0, 2, 0, thetaq, wigner) + 2 * Sqrt(3) * q2z * q3y * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Yket(2, 0, 2, 0, thetaq, wigner) - 2 * Sqrt(3) * q2y * q3z * q4z * Yket(2, 0, 2, 0, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 1, 2, 0, thetaq, wigner) + 3 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 1, 2, 0, thetaq, wigner) -
                           3 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 1, 2, 0, thetaq, wigner)) +
                      Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 2, 0, thetaq, wigner) + 6 * (q2y * q3x - q2x * q3y) * (q4x + Complex(0, 1) * q4y) * Yket(2, -1, 2, 0, thetaq, wigner) -
                           Complex(0, 2) * Sqrt(6) * q2z * q3x * q4x * Yket(2, 0, 2, 0, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3z * q4x * Yket(2, 0, 2, 0, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2z * q3y * q4y * Yket(2, 0, 2, 0, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(6) * q2y * q3z * q4y * Yket(2, 0, 2, 0, thetaq, wigner) + 6 * q2y * q3x * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - 6 * q2x * q3y * q4x * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 6) * q2y * q3x * q4y * Yket(2, 1, 2, 0, thetaq, wigner) +
                           Complex(0, 6) * q2x * q3y * q4y * Yket(2, 1, 2, 0, thetaq, wigner) - Complex(0, 6) * q2z * q3x * q4x * Yket(2, 2, 2, 0, thetaq, wigner) - 6 * q2z * q3y * q4x * Yket(2, 2, 2, 0, thetaq, wigner) + Complex(0, 6) * q2x * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) +
                           6 * q2y * q3z * q4x * Yket(2, 2, 2, 0, thetaq, wigner) - 6 * q2z * q3x * q4y * Yket(2, 2, 2, 0, thetaq, wigner) + Complex(0, 6) * q2z * q3y * q4y * Yket(2, 2, 2, 0, thetaq, wigner) + 6 * q2x * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner) -
                           Complex(0, 6) * q2y * q3z * q4y * Yket(2, 2, 2, 0, thetaq, wigner)) +
                      Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * Sqrt(2) * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -1, 2, 0, thetaq, wigner) +
                           2 * Sqrt(3) * (2 * q2y * q3x * (q4x + Complex(0, 1) * q4y) - 2 * q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(2, 0, 2, 0, thetaq, wigner) +
                           3 * Sqrt(2) *
                               ((Complex(0, 1) * q2x * q3z * q4x + q2y * q3z * q4x - q2x * q3z * q4y + Complex(0, 1) * q2y * q3z * q4y + q2z * (q3x - Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) + 2 * q2y * q3x * q4z - 2 * q2x * q3y * q4z) * Yket(2, 1, 2, 0, thetaq, wigner) +
                                2 * (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * q4z * Yket(2, 2, 2, 0, thetaq, wigner))))) /
            (6. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 11 && idxket == 16)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * Sqrt(2) * (Complex(0, 1) * q2z * q3x - q2z * q3y - Complex(0, 1) * q2x * q3z + q2y * q3z) * q4z * Yket(2, -2, 2, 2, thetaq, wigner) +
                           3 * Sqrt(2) * (Complex(0, 1) * q2x * q3z * q4x - q2y * q3z * q4x + q2z * (Complex(0, -1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + q2x * q3z * q4y + Complex(0, 1) * q2y * q3z * q4y - 2 * q2y * q3x * q4z + 2 * q2x * q3y * q4z) * Yket(2, -1, 2, 2, thetaq, wigner) +
                           4 * Sqrt(3) * q2y * q3x * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - 4 * Sqrt(3) * q2x * q3y * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2y * q3x * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 4) * Sqrt(3) * q2x * q3y * q4y * Yket(2, 0, 2, 2, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Yket(2, 0, 2, 2, thetaq, wigner) + 2 * Sqrt(3) * q2z * q3y * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Yket(2, 0, 2, 2, thetaq, wigner) - 2 * Sqrt(3) * q2y * q3z * q4z * Yket(2, 0, 2, 2, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 1, 2, 2, thetaq, wigner) + 3 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 1, 2, 2, thetaq, wigner) -
                           3 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 1, 2, 2, thetaq, wigner)) +
                      Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 2, 2, thetaq, wigner) + 6 * (q2y * q3x - q2x * q3y) * (q4x + Complex(0, 1) * q4y) * Yket(2, -1, 2, 2, thetaq, wigner) -
                           Complex(0, 2) * Sqrt(6) * q2z * q3x * q4x * Yket(2, 0, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2x * q3z * q4x * Yket(2, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2z * q3y * q4y * Yket(2, 0, 2, 2, thetaq, wigner) +
                           Complex(0, 2) * Sqrt(6) * q2y * q3z * q4y * Yket(2, 0, 2, 2, thetaq, wigner) + 6 * q2y * q3x * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - 6 * q2x * q3y * q4x * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 6) * q2y * q3x * q4y * Yket(2, 1, 2, 2, thetaq, wigner) +
                           Complex(0, 6) * q2x * q3y * q4y * Yket(2, 1, 2, 2, thetaq, wigner) - Complex(0, 6) * q2z * q3x * q4x * Yket(2, 2, 2, 2, thetaq, wigner) - 6 * q2z * q3y * q4x * Yket(2, 2, 2, 2, thetaq, wigner) + Complex(0, 6) * q2x * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) +
                           6 * q2y * q3z * q4x * Yket(2, 2, 2, 2, thetaq, wigner) - 6 * q2z * q3x * q4y * Yket(2, 2, 2, 2, thetaq, wigner) + Complex(0, 6) * q2z * q3y * q4y * Yket(2, 2, 2, 2, thetaq, wigner) + 6 * q2x * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner) -
                           Complex(0, 6) * q2y * q3z * q4y * Yket(2, 2, 2, 2, thetaq, wigner)) +
                      Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * Sqrt(2) * (Complex(0, -1) * q2z * q3x + q2z * q3y + Complex(0, 1) * q2x * q3z - q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(2, -1, 2, 2, thetaq, wigner) +
                           2 * Sqrt(3) * (2 * q2y * q3x * (q4x + Complex(0, 1) * q4y) - 2 * q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(2, 0, 2, 2, thetaq, wigner) +
                           3 * Sqrt(2) *
                               ((Complex(0, 1) * q2x * q3z * q4x + q2y * q3z * q4x - q2x * q3z * q4y + Complex(0, 1) * q2y * q3z * q4y + q2z * (q3x - Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) + 2 * q2y * q3x * q4z - 2 * q2x * q3y * q4z) * Yket(2, 1, 2, 2, thetaq, wigner) +
                                2 * (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * q4z * Yket(2, 2, 2, 2, thetaq, wigner))))) /
            (6. * Sqrt(5));
        return mtx;
    }
    else if (idxbra == 12 && idxket == 1)
    {
        mtx = -2 * Ftpe3 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 0, 0, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 12 && idxket == 7)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               (Sqrt(2) *
                    (q2z * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, -3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y - Complex(0, 2) * q3y * q4y - Complex(0, 2) * q3z * q4z) +
                     q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Yket(1, -1, 1, 1, thetaq, wigner) -
                Complex(0, 2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                Sqrt(2) *
                    (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                     q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Yket(1, 1, 1, 1, thetaq, wigner))) /
              3.;
        return mtx;
    }
    else if (idxbra == 12 && idxket == 8)
    {
        mtx =
            -0.3333333333333333 *
            (Ftpe3 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
             ((q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) -
              Complex(0, 1) * Sqrt(2) * (q2z * q3x * q4x + q2x * q3z * q4x + q2z * q3y * q4y + q2y * q3z * q4y - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, 0, 1, 1, thetaq, wigner) +
              (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                  Yket(1, 1, 1, 1, thetaq, wigner)));
        return mtx;
    }
    else if (idxbra == 12 && idxket == 9)
    {
        mtx = (Ftpe3 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
               ((Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Yket(2, -2, 1, 1, thetaq, wigner) +
                (q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (-(q3y * q4x) + 2 * q3x * q4y + Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (-(q3x * q4x) - Complex(0, 2) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3z * q4z)) *
                    Yket(2, -1, 1, 1, thetaq, wigner) -
                Sqrt(6) * q2z * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2z * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3z * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + q2y * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) +
                q2x * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 2) * q2y * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * q2z * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - 2 * q2x * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) +
                Complex(0, 1) * q2y * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 1) * q2x * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + 2 * q2z * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * q2z * q3x * q4z * Yket(2, 1, 1, 1, thetaq, wigner) -
                q2z * q3y * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * q2x * q3z * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - q2y * q3z * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * q2z * q3x * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - q2z * q3y * q4x * Yket(2, 2, 1, 1, thetaq, wigner) -
                Complex(0, 1) * q2x * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - q2y * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - q2z * q3x * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 1) * q2z * q3y * q4y * Yket(2, 2, 1, 1, thetaq, wigner) - q2x * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) +
                Complex(0, 1) * q2y * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 2) * q2x * q3x * q4z * Yket(2, 2, 1, 1, thetaq, wigner) + 2 * q2y * q3x * q4z * Yket(2, 2, 1, 1, thetaq, wigner) + 2 * q2x * q3y * q4z * Yket(2, 2, 1, 1, thetaq, wigner) -
                Complex(0, 2) * q2y * q3y * q4z * Yket(2, 2, 1, 1, thetaq, wigner))) /
              Sqrt(5);
        return mtx;
    }
    else if (idxbra == 12 && idxket == 10)
    {
        mtx = -2 * Ftpe3 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(0, 0, 2, 2, thetaq, wigner);
        return mtx;
    }
    else if (idxbra == 12 && idxket == 11)
    {
        mtx = -((Ftpe3 * Ybra(0, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                 (Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                  Complex(0, 2) * (q3z * (q2x * q4x + q2y * q4y) - q2z * (q3x * q4x + q3y * q4y)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                  Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner))) /
                Sqrt(3));
        return mtx;
    }
    else if (idxbra == 13 && idxket == 1)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (-(q3z * (q2x * q4x + q2y * q4y)) + q2z * (q3x * q4x + q3y * q4y)) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 13 && idxket == 6)
    {
        mtx = (Ftpe3 *
               (-(Sqrt(2) *
                  (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                   q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                  Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) -
                Complex(0, 2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) *
                    (q2z * (Complex(0, -1) * q3z * q4x + q3z * q4y + Complex(0, 2) * q3x * q4z - 2 * q3y * q4z) + q2x * (-2 * q3y * q4x + Complex(0, 2) * q3y * q4y + q3x * (Complex(0, 3) * q4x + q4y) + Complex(0, 2) * q3z * q4z) -
                     q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                    Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              3.;
        return mtx;
    }
    else if (idxbra == 13 && idxket == 7)
    {
        mtx = (Ftpe3 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, 2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) -
                                                                                     Sqrt(2) *
                                                                                         (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                                                                                          q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                                                                                         Yket(1, 0, 1, 1, thetaq, wigner)) +
                        Sqrt(2) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((q2z * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, -3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y - Complex(0, 2) * q3y * q4y - Complex(0, 2) * q3z * q4z) +
                              q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                                 Yket(1, -1, 1, 1, thetaq, wigner) -
                             (q2z * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, 3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y + Complex(0, 2) * q3y * q4y + Complex(0, 2) * q3z * q4z) +
                              q2y * (2 * q3x * q4x - Complex(0, 1) * q3y * q4x + Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                                 Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) *
                                                                                        (q2z * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z) + q2x * (Complex(0, -3) * q3x * q4x + 2 * q3y * q4x - q3x * q4y - Complex(0, 2) * q3y * q4y - Complex(0, 2) * q3z * q4z) +
                                                                                         q2y * (2 * q3x * q4x + Complex(0, 1) * q3y * q4x - Complex(0, 2) * q3x * q4y + 3 * q3y * q4y + 2 * q3z * q4z)) *
                                                                                        Yket(1, 0, 1, 1, thetaq, wigner) -
                                                                                    Complex(0, 2) * (2 * q2x * q3z * q4x + 2 * q2y * q3z * q4y - q2x * q3x * q4z - q2y * q3y * q4z + q2z * (2 * q3x * q4x + 2 * q3y * q4y + 3 * q3z * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (3. * Sqrt(3));
        return mtx;
    }
    else if (idxbra == 13 && idxket == 8)
    {
        mtx = (Ftpe3 * (Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (3 * Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) - 2 * (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             (2 * q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 8) * q3y * q4x + 4 * q3x * (q4x - Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (4 * q3y * q4x - 8 * q3x * q4y - Complex(0, 4) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                                 Yket(1, 0, 1, 1, thetaq, wigner) -
                             Complex(0, 1) * Sqrt(2) * (q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, 1) * Sqrt(2) * (q2z * q3x * (q4x + Complex(0, 3) * q4y) + q2x * q3z * (q4x + Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, -3) * q4x + q4y) + q2y * q3z * (Complex(0, -3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             2 * (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (-2 * q3y * (q4x + Complex(0, 1) * q4y) + 4 * q3x * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 4) * q3y * q4x - 2 * q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                                 Yket(1, 0, 1, 1, thetaq, wigner) +
                             3 * Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)) -
                        2 * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((2 * q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                                 Yket(1, -1, 1, 1, thetaq, wigner) +
                             3 * Sqrt(2) * (q2z * q3y * q4x + q2y * q3z * q4x - q2z * q3x * q4y - q2x * q3z * q4y) * Yket(1, 0, 1, 1, thetaq, wigner) -
                             (2 * q2z * (Complex(0, 2) * q3z * q4x + 2 * q3z * q4y - Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                                 Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (6. * Sqrt(3));
        return mtx;
    }
    else if (idxbra == 13 && idxket == 9)
    {
        mtx = (Ftpe3 *
               (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                    (6 * Sqrt(2) * (Complex(0, 1) * (q2x + Complex(0, 1) * q2y) * q3z * q4z + q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z)) * Yket(2, -2, 1, 1, thetaq, wigner) +
                     Complex(0, 3) * Sqrt(2) * (q2z * q3x * (q4x + Complex(0, 3) * q4y) + q2x * q3z * (q4x + Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, -3) * q4x + q4y) + q2y * q3z * (Complex(0, -3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(2, -1, 1, 1, thetaq, wigner) -
                     4 * Sqrt(3) * q2y * q3x * q4x * Yket(2, 0, 1, 1, thetaq, wigner) - 4 * Sqrt(3) * q2x * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2y * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2z * q3z * q4x * Yket(2, 0, 1, 1, thetaq, wigner) +
                     8 * Sqrt(3) * q2x * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2y * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2x * q3y * q4y * Yket(2, 0, 1, 1, thetaq, wigner) - 4 * Sqrt(3) * q2z * q3z * q4y * Yket(2, 0, 1, 1, thetaq, wigner) +
                     Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Yket(2, 0, 1, 1, thetaq, wigner) + 2 * Sqrt(3) * q2z * q3y * q4z * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Yket(2, 0, 1, 1, thetaq, wigner) + 2 * Sqrt(3) * q2y * q3z * q4z * Yket(2, 0, 1, 1, thetaq, wigner) +
                     Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) +
                     3 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner) -
                     Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2y * q3x * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2x * q3y * q4z * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Yket(2, 1, 1, 1, thetaq, wigner)) +
                2 * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                    ((Complex(0, 3) * q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + 3 * (q2x + Complex(0, 1) * q2y) * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Yket(2, -2, 1, 1, thetaq, wigner) -
                     3 * (Complex(0, 2) * q2y * q3y * q4x + q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2x * q3y * (q4x - Complex(0, 1) * q4y) - 2 * q2x * q3x * q4y) * Yket(2, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2z * q3x * q4x * Yket(2, 0, 1, 1, thetaq, wigner) +
                     Complex(0, 1) * Sqrt(6) * q2x * q3z * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2z * q3y * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * q4y * Yket(2, 0, 1, 1, thetaq, wigner) -
                     Complex(0, 2) * Sqrt(6) * q2x * q3x * q4z * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3y * q4z * Yket(2, 0, 1, 1, thetaq, wigner) - 3 * q2y * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - 3 * q2x * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) +
                     Complex(0, 6) * q2y * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + 6 * q2x * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2y * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2x * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) +
                     Complex(0, 3) * q2z * q3x * q4x * Yket(2, 2, 1, 1, thetaq, wigner) + 3 * q2z * q3y * q4x * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 3) * q2x * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) + 3 * q2y * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) +
                     3 * q2z * q3x * q4y * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + 3 * q2x * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 3) * q2y * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) -
                     Complex(0, 6) * q2x * q3x * q4z * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * q2y * q3x * q4z * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * q2x * q3y * q4z * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 6) * q2y * q3y * q4z * Yket(2, 2, 1, 1, thetaq, wigner)) +
                Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                    (3 * Sqrt(2) * (Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Yket(2, -1, 1, 1, thetaq, wigner) +
                     2 * Sqrt(3) * (q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (-2 * q3y * q4x + 4 * q3x * q4y + Complex(0, 2) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (-2 * q3x * q4x - Complex(0, 4) * q3y * q4x + Complex(0, 2) * q3x * q4y + q3z * q4z)) *
                         Yket(2, 0, 1, 1, thetaq, wigner) +
                     Complex(0, 3) * Sqrt(2) *
                         ((q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(2, 1, 1, 1, thetaq, wigner) -
                          2 * ((q2x - Complex(0, 1) * q2y) * q3z * q4z + q2z * (-2 * q3z * q4x + Complex(0, 2) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z)) * Yket(2, 2, 1, 1, thetaq, wigner))))) /
              (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 13 && idxket == 10)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (-(q3z * (q2x * q4x + q2y * q4y)) + q2z * (q3x * q4x + q3y * q4y)) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 13 && idxket == 11)
    {
        mtx = (Ftpe3 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * (q2z * (q3x - Complex(0, 1) * q3y) * (Complex(0, -1) * q4x + q4y) + q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y - q3x * q4z) + q2x * (Complex(0, 1) * q3z * q4x - q3z * q4y + q3y * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                                                                                     Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * q4z + q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) - q2y * (q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                    2 * (q2z * (q3x + Complex(0, 1) * q3y) * (Complex(0, 1) * q4x + q4y) + q2y * (q3z * q4x - Complex(0, 1) * q3z * q4y - q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x - q3z * q4y + q3y * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                        Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * q4z + q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) - q2y * (q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) -
                                                                                    2 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y + q2y * q3x * q4z - q2x * q3y * q4z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                    Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
              3.;
        return mtx;
    }
    else if (idxbra == 14 && idxket == 1)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (-(q3z * (q2x * q4x + q2y * q4y)) + q2z * (q3x * q4x + q3y * q4y)) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              Sqrt(6);
        return mtx;
    }
    else if (idxbra == 14 && idxket == 6)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                    Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (q2z * q3x * q4x + q2x * q3z * q4x + q2z * q3y * q4y + q2y * q3z * q4y - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (Complex(0, 2) * q3z * q4x - 2 * q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) + q3z * q4z)) *
                    Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 7)
    {
        mtx = (Ftpe3 * (Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((Complex(0, 3) * q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + 3 * (q2x + Complex(0, 1) * q2y) * (Complex(0, 1) * q3z * q4x - q3z * q4y - Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) -
                             Sqrt(2) * (2 * q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 2) * q3y * q4x + q3x * (q4x - Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y - Complex(0, 1) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                                 Yket(1, 0, 1, 1, thetaq, wigner) -
                             Complex(0, 1) * (q2z * q3x * (q4x + Complex(0, 3) * q4y) + q2x * q3z * (q4x + Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, -3) * q4x + q4y) + q2y * q3z * (Complex(0, -3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, 1) * (q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (2 * q2z * (Complex(0, 2) * q3z * q4x + 2 * q3z * q4y - Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) - 2 * q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                                 Yket(1, 0, 1, 1, thetaq, wigner) +
                             3 * (Complex(0, -1) * q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (Complex(0, -1) * q3z * q4x - q3z * q4y + Complex(0, 2) * q3x * q4z + 2 * q3y * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Sqrt(2) * (q2z * (Complex(0, -2) * q3z * q4x + 2 * q3z * q4y + Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, 4) * q3y * q4x + 2 * q3x * (q4x - Complex(0, 1) * q4y) - q3z * q4z) + q2x * (2 * q3y * q4x - 4 * q3x * q4y - Complex(0, 2) * q3y * q4y + Complex(0, 1) * q3z * q4z)) *
                                 Yket(1, -1, 1, 1, thetaq, wigner) +
                             6 * (q2z * q3y * q4x + q2y * q3z * q4x - q2z * q3x * q4y - q2x * q3z * q4y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (-2 * q3y * (q4x + Complex(0, 1) * q4y) + 4 * q3x * q4y + Complex(0, 1) * q3z * q4z) + q2y * (Complex(0, 4) * q3y * q4x - 2 * q3x * (q4x + Complex(0, 1) * q4y) + q3z * q4z)) *
                                 Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (3. * Sqrt(6));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 8)
    {
        mtx = (Ftpe3 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Complex(0, -2) * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y + 3 * q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                                                                                     Sqrt(2) *
                                                                                         (q2z * (Complex(0, 1) * q3z * q4x + q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y + q3x * (Complex(0, 3) * q4x + q4y) + Complex(0, 1) * q3z * q4z) +
                                                                                          q2y * (q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + 3 * q3y * q4y + q3z * q4z)) *
                                                                                         Yket(1, 0, 1, 1, thetaq, wigner)) +
                        Complex(0, 1) * (Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) *
                                                                                                         (q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y + Complex(0, 1) * q3z * q4z) +
                                                                                                          q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y + q3z * q4z)) *
                                                                                                         Yket(1, 0, 1, 1, thetaq, wigner) +
                                                                                                     2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y + 3 * q3z * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                                         Sqrt(2) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                                             ((q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y + Complex(0, 1) * q3z * q4z) +
                                               q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y + q3z * q4z)) *
                                                  Yket(1, -1, 1, 1, thetaq, wigner) +
                                              (q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y + q3z * q4z) -
                                               Complex(0, 1) * (q2z * (Complex(0, 1) * q3z * q4x + q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z) + q2y * (q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + 3 * q3y * q4y + q3z * q4z))) *
                                                  Yket(1, 1, 1, 1, thetaq, wigner))))) /
              (6. * Sqrt(3));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 9)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (Complex(0, 6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Yket(2, -2, 1, 1, thetaq, wigner) -
                           Complex(0, 6) * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y - 3 * q3z * q4z)) * Yket(2, -1, 1, 1, thetaq, wigner) + Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Yket(2, 0, 1, 1, thetaq, wigner) +
                           Sqrt(6) * q2y * q3x * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Yket(2, 0, 1, 1, thetaq, wigner) +
                           Sqrt(6) * q2x * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + 3 * Sqrt(6) * q2y * q3y * q4y * Yket(2, 0, 1, 1, thetaq, wigner) -
                           5 * Sqrt(6) * q2z * q3z * q4y * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Yket(2, 0, 1, 1, thetaq, wigner) - 5 * Sqrt(6) * q2z * q3y * q4z * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Yket(2, 0, 1, 1, thetaq, wigner) -
                           5 * Sqrt(6) * q2y * q3z * q4z * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 12) * q2z * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + 12 * q2z * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 12) * q2x * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) +
                           12 * q2y * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + 12 * q2z * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2z * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + 12 * q2x * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner) -
                           Complex(0, 12) * q2y * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 12) * q2x * q3x * q4z * Yket(2, 1, 1, 1, thetaq, wigner) + 12 * q2y * q3x * q4z * Yket(2, 1, 1, 1, thetaq, wigner) + 12 * q2x * q3y * q4z * Yket(2, 1, 1, 1, thetaq, wigner) -
                           Complex(0, 12) * q2y * q3y * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 18) * q2x * q3x * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - 18 * q2y * q3x * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - 18 * q2x * q3y * q4x * Yket(2, 2, 1, 1, thetaq, wigner) +
                           Complex(0, 18) * q2y * q3y * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - 18 * q2x * q3x * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 18) * q2y * q3x * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 18) * q2x * q3y * q4y * Yket(2, 2, 1, 1, thetaq, wigner) +
                           18 * q2y * q3y * q4y * Yket(2, 2, 1, 1, thetaq, wigner)) +
                      Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (6 * Sqrt(2) * (q2z * (Complex(0, -1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) + (Complex(0, -1) * q2x + q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Yket(2, -2, 1, 1, thetaq, wigner) +
                           3 * Sqrt(2) *
                               (3 * q2z * (Complex(0, -1) * q3z * q4x + q3z * q4y - Complex(0, 1) * q3x * q4z + q3y * q4z) + Complex(0, 1) * q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 3 * q3z * q4z) +
                                q2y * (-(q3x * q4x) + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y - 3 * q3y * q4y + 3 * q3z * q4z)) *
                               Yket(2, -1, 1, 1, thetaq, wigner) +
                           Complex(0, 8) * Sqrt(3) * q2z * q3x * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2x * q3z * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2z * q3y * q4y * Yket(2, 0, 1, 1, thetaq, wigner) +
                           Complex(0, 8) * Sqrt(3) * q2y * q3z * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2x * q3x * q4z * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2y * q3y * q4z * Yket(2, 0, 1, 1, thetaq, wigner) -
                           Complex(0, 12) * Sqrt(3) * q2z * q3z * q4z * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2x * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2y * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2x * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * q2y * q3y * q4x * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2z * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2x * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - 9 * Sqrt(2) * q2y * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) + 9 * Sqrt(2) * q2z * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner) +
                           Complex(0, 9) * Sqrt(2) * q2z * q3x * q4z * Yket(2, 1, 1, 1, thetaq, wigner) + 9 * Sqrt(2) * q2z * q3y * q4z * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2x * q3z * q4z * Yket(2, 1, 1, 1, thetaq, wigner) + 9 * Sqrt(2) * q2y * q3z * q4z * Yket(2, 1, 1, 1, thetaq, wigner) -
                           Complex(0, 6) * Sqrt(2) * q2z * q3x * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2z * q3y * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2x * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2y * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) -
                           6 * Sqrt(2) * q2z * q3x * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2z * q3y * q4y * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2x * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2y * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) -
                           Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2y * q3x * q4z * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2x * q3y * q4z * Yket(2, 2, 1, 1, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Yket(2, 2, 1, 1, thetaq, wigner)) +
                      Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (Complex(0, 18) * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Yket(2, -2, 1, 1, thetaq, wigner) +
                           Complex(0, 12) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z)) * Yket(2, -1, 1, 1, thetaq, wigner) -
                           Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2y * q3x * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Yket(2, 0, 1, 1, thetaq, wigner) +
                           Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Yket(2, 0, 1, 1, thetaq, wigner) + Sqrt(6) * q2x * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Yket(2, 0, 1, 1, thetaq, wigner) -
                           Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + 3 * Sqrt(6) * q2y * q3y * q4y * Yket(2, 0, 1, 1, thetaq, wigner) - 5 * Sqrt(6) * q2z * q3z * q4y * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Yket(2, 0, 1, 1, thetaq, wigner) -
                           5 * Sqrt(6) * q2z * q3y * q4z * Yket(2, 0, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Yket(2, 0, 1, 1, thetaq, wigner) - 5 * Sqrt(6) * q2y * q3z * q4z * Yket(2, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3x * q4x * Yket(2, 1, 1, 1, thetaq, wigner) -
                           Complex(0, 6) * q2x * q3z * q4x * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3y * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2y * q3z * q4y * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2x * q3x * q4z * Yket(2, 1, 1, 1, thetaq, wigner) -
                           Complex(0, 6) * q2y * q3y * q4z * Yket(2, 1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2z * q3z * q4z * Yket(2, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3z * q4x * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * q2z * q3z * q4y * Yket(2, 2, 1, 1, thetaq, wigner) -
                           Complex(0, 6) * q2z * q3x * q4z * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * q2z * q3y * q4z * Yket(2, 2, 1, 1, thetaq, wigner) - Complex(0, 6) * q2x * q3z * q4z * Yket(2, 2, 1, 1, thetaq, wigner) - 6 * q2y * q3z * q4z * Yket(2, 2, 1, 1, thetaq, wigner)))) /
            (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 14 && idxket == 10)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2y * q3z * q4z + q2x * (-(q3y * q4x) + Complex(0, 1) * q3y * q4y + Complex(0, 1) * q3z * q4z)) * Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * (-(q3z * (q2x * q4x + q2y * q4y)) + q2z * (q3x * q4x + q3y * q4y)) * Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2y * q3x * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + q2y * q3z * q4z - q2x * (q3y * (q4x + Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              Sqrt(6);
        return mtx;
    }
    else if (idxbra == 14 && idxket == 11)
    {
        mtx =
            (Ftpe3 * (Ybra(1, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * ((q2z * (Complex(0, 1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) - q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x + q3z * q4y + 2 * q3y * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                                                                                   Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + 2 * q2z * (Complex(0, 1) * q3x + q3y) * q4z - 2 * q2y * q3z * q4z - q2x * (q3y * (q4x - Complex(0, 1) * q4y) + Complex(0, 2) * q3z * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                   3 * (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * (q4x - Complex(0, 1) * q4y) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                      Ybra(1, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                          (3 * (Complex(0, 1) * q2z * q3x - q2z * q3y - Complex(0, 1) * q2x * q3z + q2y * q3z) * (q4x + Complex(0, 1) * q4y) * Yket(1, -1, 2, 2, thetaq, wigner) +
                           Sqrt(2) * (-(q2y * q3x * (q4x + Complex(0, 1) * q4y)) + Complex(0, 2) * q2z * (q3x + Complex(0, 1) * q3y) * q4z + 2 * q2y * q3z * q4z + q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y - Complex(0, 2) * q3z * q4z)) * Yket(1, 0, 2, 2, thetaq, wigner) +
                           (q2z * (Complex(0, -1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) - q2y * (q3z * q4x - Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, 1) * q3z * q4x + q3z * q4y + 2 * q3y * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)) +
                      Ybra(1, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (Sqrt(2) * (2 * q2y * q3x * (q4x + Complex(0, 1) * q4y) - 2 * q2x * q3y * (q4x + Complex(0, 1) * q4y) + q2z * (Complex(0, -1) * q3x + q3y) * q4z + Complex(0, 1) * q2x * q3z * q4z - q2y * q3z * q4z) * Yket(1, -1, 2, 2, thetaq, wigner) -
                                                                                  2 * (q2z * q3y * q4x - q2y * q3z * q4x - q2z * q3x * q4y + q2x * q3z * q4y - 2 * q2y * q3x * q4z + 2 * q2x * q3y * q4z) * Yket(1, 0, 2, 2, thetaq, wigner) +
                                                                                  Sqrt(2) * (q2z * (Complex(0, -1) * q3x - q3y) * q4z + q2x * (2 * q3y * q4x - Complex(0, 2) * q3y * q4y + Complex(0, 1) * q3z * q4z) + q2y * (-2 * q3x * q4x + Complex(0, 2) * q3x * q4y + q3z * q4z)) * Yket(1, 1, 2, 2, thetaq, wigner)))) /
            (3. * Sqrt(2));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 1)
    {
        mtx = (Ftpe3 *
               (3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * (Complex(0, 1) * q4x + q4y) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z - q2y * q3z * q4z - q2x * (q3y * (q4x - Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 4 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 4 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 15 && idxket == 6)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (Complex(0, 2) * q3z * q4x + 2 * q3z * q4y - Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) - q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z)) *
                    Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 2) * Sqrt(2) * q2z * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2z * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              Sqrt(10);
        return mtx;
    }
    else if (idxbra == 15 && idxket == 7)
    {
        mtx = (Ftpe3 * (2 * Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 4) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        4 * Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2z * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 3 * q2z * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2x * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 3 * q2y * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        3 * q2z * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        3 * q2x * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2y * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 6 * q2y * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        6 * q2x * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 4) * Sqrt(3) * q2x * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2y * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * ((Complex(0, 1) * q2x + q2y) * q3z * q4z + q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                        2 * Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 4) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        4 * Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2z * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 9 * q2z * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2x * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 9 * q2y * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        9 * q2z * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 3) * q2z * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        9 * q2x * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 3) * q2y * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2z * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 12 * q2z * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3x * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2z * q3y * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2y * q3z * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        3 * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, 1) * (q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (Complex(0, -2) * q2y * q3y * q4x + q2y * q3x * (q4x + Complex(0, 1) * q4y) + q2x * q3y * (q4x + Complex(0, 1) * q4y) - 2 * q2x * q3x * q4y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                             (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (3. * Sqrt(30));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 8)
    {
        mtx = (Ftpe3 * (Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        3 * Sqrt(6) * q2y * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 5 * Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2z * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 12 * q2z * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2x * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 12 * q2y * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        12 * q2z * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2z * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        12 * q2x * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2y * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2x * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 12 * q2y * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        12 * q2x * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2y * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 18) * q2x * q3x * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 18 * q2y * q3x * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        18 * q2x * q3y * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2y * q3y * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        18 * q2x * q3x * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2y * q3x * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 18) * q2x * q3y * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 18 * q2y * q3y * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 8) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 8) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 8) * Sqrt(3) * q2x * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2y * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 12) * Sqrt(3) * q2z * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2x * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2z * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 9 * Sqrt(2) * q2y * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        9 * Sqrt(2) * q2z * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        9 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        9 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 3 * Sqrt(6) * q2y * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        5 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        5 * Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        5 * Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2y * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2z * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2z * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3x * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2z * q3y * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2y * q3z * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        6 * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((q2z * q3z * (Complex(0, 1) * q4x + q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z + (Complex(0, 1) * q2x + q2y) * q3z * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) -
                             Complex(0, 1) * Sqrt(2) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                             3 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (Complex(0, 1) * q4x + q4y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        3 * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, -2) * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y - 3 * q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) *
                                 (Complex(0, -3) * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2y * (q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + 3 * q3y * q4y - 3 * q3z * q4z) +
                                  q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y + q3x * (Complex(0, 3) * q4x + q4y) - Complex(0, 3) * q3z * q4z)) *
                                 Yket(1, 0, 1, 1, thetaq, wigner) +
                             4 * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 15 && idxket == 9)
    {
        mtx = (Complex(0, -0.1) * Ftpe3 *
               (12 * q2z * q3z * q4z * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                3 * q2z * q3z * q4z * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                3 * q2z * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                12 * q2z * q3z * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 15 && idxket == 10)
    {
        mtx = (Ftpe3 *
               (3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * (Complex(0, 1) * q4x + q4y) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z - q2y * q3z * q4z - q2x * (q3y * (q4x - Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 4 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 4 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 15 && idxket == 11)
    {
        mtx = (Ftpe3 * (2 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - 2 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(2) * q2z * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(3) * q2z * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(3) * q2y * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2z * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(3) * q2x * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * q2z * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * q2x * q3z * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * q2z * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * q2y * q3z * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2y * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2y * q3z * q4x * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Sqrt(6) * q2z * q3x * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * q4y * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(3) * (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * Ybra(2, -2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * q4z * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * (q4x - Complex(0, 1) * q4y) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        2 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2z * q3y * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3z * q4z * Ybra(2, 0, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Sqrt(3) * q2z * q3y * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(3) * q2y * q3z * q4x * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(3) * q2z * q3x * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        2 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * Sqrt(3) * q2z * q3y * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 2 * Sqrt(3) * q2y * q3z * q4z * Ybra(2, 2, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(3) * Ybra(2, -1, 2, 0, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((q2z * (Complex(0, 1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) - q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x + q3z * q4y + 2 * q3y * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                             (q4x - Complex(0, 1) * q4y) * (Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 2, 2, thetaq, wigner) + (q2z * (Complex(0, 1) * q3x + q3y) - Complex(0, 1) * q2x * q3z - q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner))))) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 16 && idxket == 1)
    {
        mtx = (Ftpe3 *
               (3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * (Complex(0, 1) * q4x + q4y) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z - q2y * q3z * q4z - q2x * (q3y * (q4x - Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 4 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 4 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 0, 0, thetaq, wigner)) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 16 && idxket == 6)
    {
        mtx = (Ftpe3 *
               (Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * (q2z * (Complex(0, 2) * q3z * q4x + 2 * q3z * q4y - Complex(0, 1) * q3x * q4z - q3y * q4z) + q2y * (Complex(0, -2) * q3y * q4x + q3x * (q4x + Complex(0, 1) * q4y) - q3z * q4z) + q2x * (q3y * q4x - 2 * q3x * q4y + Complex(0, 1) * q3y * q4y - Complex(0, 1) * q3z * q4z)) *
                    Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 2) * Sqrt(2) * q2z * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2z * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 1) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 1) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 1) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 1) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 2) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 2 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 2) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 1, 1, thetaq, wigner)) /
              Sqrt(10);
        return mtx;
    }
    else if (idxbra == 16 && idxket == 7)
    {
        mtx = (Ftpe3 * (2 * Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 4) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        4 * Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 2) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2z * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 3 * q2z * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2x * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 3 * q2y * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        3 * q2z * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2z * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        3 * q2x * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 3) * q2y * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 6 * q2y * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        6 * q2x * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 4) * Sqrt(3) * q2x * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 4) * Sqrt(3) * q2y * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - 6 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (2 * ((Complex(0, 1) * q2x + q2y) * q3z * q4z + q2z * (Complex(0, -2) * q3z * q4x - 2 * q3z * q4y + Complex(0, 1) * q3x * q4z + q3y * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner)) +
                        2 * Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 4) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        4 * Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 2) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 2 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2z * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 9 * q2z * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 3) * q2x * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 9 * q2y * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        9 * q2z * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 3) * q2z * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        9 * q2x * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 3) * q2y * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2z * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 12 * q2z * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3x * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2z * q3y * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2y * q3z * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        3 * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, 1) * (q2z * q3x * (q4x - Complex(0, 3) * q4y) + q2x * q3z * (q4x - Complex(0, 3) * q4y) + q2z * q3y * (Complex(0, 3) * q4x + q4y) + q2y * q3z * (Complex(0, 3) * q4x + q4y) - 2 * q2x * q3x * q4z - 2 * q2y * q3y * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) * (Complex(0, -2) * q2y * q3y * q4x + q2y * q3x * (q4x + Complex(0, 1) * q4y) + q2x * q3y * (q4x + Complex(0, 1) * q4y) - 2 * q2x * q3x * q4y) * Yket(1, 0, 1, 1, thetaq, wigner) +
                             (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * q4x - Complex(0, 1) * q3z * q4y - 2 * q3x * q4z + Complex(0, 2) * q3y * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (3. * Sqrt(30));
        return mtx;
    }
    else if (idxbra == 16 && idxket == 8)
    {
        mtx = (Ftpe3 * (Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        3 * Sqrt(6) * q2y * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 5 * Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 5 * Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2z * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 12 * q2z * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2x * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 12 * q2y * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        12 * q2z * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2z * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        12 * q2x * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2y * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 12) * q2x * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 12 * q2y * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        12 * q2x * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - Complex(0, 12) * q2y * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) -
                        Complex(0, 18) * q2x * q3x * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + 18 * q2y * q3x * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        18 * q2x * q3y * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2y * q3y * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        18 * q2x * q3x * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2y * q3x * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 18) * q2x * q3y * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) - 18 * q2y * q3y * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 1, 1, thetaq, wigner) +
                        Complex(0, 8) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 8) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 8) * Sqrt(3) * q2x * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 8) * Sqrt(3) * q2y * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 12) * Sqrt(3) * q2z * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 9) * Sqrt(2) * q2x * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2y * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2z * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        3 * Sqrt(2) * q2x * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 9 * Sqrt(2) * q2y * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        9 * Sqrt(2) * q2z * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        9 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + Complex(0, 9) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) -
                        9 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 6) * Sqrt(2) * q2x * q3x * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        6 * Sqrt(2) * q2y * q3x * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) + 6 * Sqrt(2) * q2x * q3y * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) +
                        Complex(0, 6) * Sqrt(2) * q2y * q3y * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 1, 1, thetaq, wigner) - Complex(0, 3) * Sqrt(6) * q2x * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2y * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2y * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Sqrt(6) * q2x * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - 3 * Sqrt(6) * q2y * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        5 * Sqrt(6) * q2z * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2z * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        5 * Sqrt(6) * q2z * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 5) * Sqrt(6) * q2x * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        5 * Sqrt(6) * q2y * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2z * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2y * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) - Complex(0, 6) * q2x * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2y * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + Complex(0, 18) * q2z * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2z * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2z * q3x * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2z * q3y * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) -
                        Complex(0, 6) * q2x * q3z * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) + 6 * q2y * q3z * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 1, 1, thetaq, wigner) +
                        6 * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((q2z * q3z * (Complex(0, 1) * q4x + q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z + (Complex(0, 1) * q2x + q2y) * q3z * q4z) * Yket(1, -1, 1, 1, thetaq, wigner) -
                             Complex(0, 1) * Sqrt(2) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Yket(1, 0, 1, 1, thetaq, wigner) +
                             3 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (Complex(0, 1) * q4x + q4y) * Yket(1, 1, 1, 1, thetaq, wigner)) +
                        3 * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            (Complex(0, -2) * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z + q2z * (q3x * q4x + q3y * q4y - 3 * q3z * q4z)) * Yket(1, -1, 1, 1, thetaq, wigner) +
                             Sqrt(2) *
                                 (Complex(0, -3) * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2y * (q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + 3 * q3y * q4y - 3 * q3z * q4z) +
                                  q2x * (q3y * q4x + Complex(0, 1) * q3y * q4y + q3x * (Complex(0, 3) * q4x + q4y) - Complex(0, 3) * q3z * q4z)) *
                                 Yket(1, 0, 1, 1, thetaq, wigner) +
                             4 * (q2z * (Complex(0, 1) * q3x + q3y) * (q4x - Complex(0, 1) * q4y) + (Complex(0, 1) * q2x + q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Yket(1, 1, 1, 1, thetaq, wigner)))) /
              (6. * Sqrt(15));
        return mtx;
    }
    else if (idxbra == 16 && idxket == 9)
    {
        mtx = (Complex(0, -0.1) * Ftpe3 *
               (12 * q2z * q3z * q4z * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -2, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                3 * q2z * q3z * q4z * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) -
                6 * (q2x + Complex(0, 1) * q2y) * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, -1, 1, 1, thetaq, wigner) +
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x + Complex(0, 1) * q3z * q4y + q3x * q4z + Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x + Complex(0, 1) * q3y * q4x + Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, 1) * q3x * q4x + q3y * q4x + q3x * q4y + Complex(0, 3) * q3y * q4y - Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x + Complex(0, 1) * q3y) * (q4x + Complex(0, 1) * q4y) + (q2x + Complex(0, 1) * q2y) * (q3z * (q4x + Complex(0, 1) * q4y) + (q3x + Complex(0, 1) * q3y) * q4z)) * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 0, 1, 1, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                Sqrt(6) * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                Sqrt(6) *
                    (-2 * q2z * (q3z * q4x - Complex(0, 1) * q3z * q4y + q3x * q4z - Complex(0, 1) * q3y * q4z) + q2x * (3 * q3x * q4x - Complex(0, 1) * q3y * q4x - Complex(0, 1) * q3x * q4y + q3y * q4y - 2 * q3z * q4z) +
                     q2y * (Complex(0, -1) * q3x * q4x + q3y * q4x + q3x * q4y - Complex(0, 3) * q3y * q4y + Complex(0, 2) * q3z * q4z)) *
                    Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                3 * q2z * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                3 * (2 * (q2x * q3z * q4x + q2y * q3z * q4y + q2x * q3x * q4z + q2y * q3y * q4z) + q2z * (2 * q3x * q4x + 2 * q3y * q4y - 3 * q3z * q4z)) * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x + Complex(0, 1) * q4y) + q2z * (q3x + Complex(0, 1) * q3y) * q4z + (q2x + Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 1, 1, 1, thetaq, wigner) -
                6 * (q2x - Complex(0, 1) * q2y) * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                2 * Sqrt(6) * (q2z * (q3x - Complex(0, 1) * q3y) * (q4x - Complex(0, 1) * q4y) + (q2x - Complex(0, 1) * q2y) * (q3z * (q4x - Complex(0, 1) * q4y) + (q3x - Complex(0, 1) * q3y) * q4z)) * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                6 * (q2z * q3z * (q4x - Complex(0, 1) * q4y) + q2z * (q3x - Complex(0, 1) * q3y) * q4z + (q2x - Complex(0, 1) * q2y) * q3z * q4z) * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner) -
                12 * q2z * q3z * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(2, 2, 1, 1, thetaq, wigner))) /
              Sqrt(3);
        return mtx;
    }
    else if (idxbra == 16 && idxket == 10)
    {
        mtx = (Ftpe3 *
               (3 * Sqrt(2) * (q2z * q3x - Complex(0, 1) * q2z * q3y - q2x * q3z + Complex(0, 1) * q2y * q3z) * (Complex(0, 1) * q4x + q4y) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * (q2y * q3x * (q4x - Complex(0, 1) * q4y) + q2z * (Complex(0, 1) * q3x + q3y) * q4z - q2y * q3z * q4z - q2x * (q3y * (q4x - Complex(0, 1) * q4y) + Complex(0, 1) * q3z * q4z)) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2z * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2y * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 2 * Sqrt(3) * q2z * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                2 * Sqrt(3) * q2x * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 4 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 4 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                3 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                Complex(0, 3) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2z * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + 3 * Sqrt(2) * q2y * q3z * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - Complex(0, 3) * Sqrt(2) * q2z * q3x * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3y * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2x * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2y * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) +
                3 * Sqrt(2) * q2z * q3x * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) + Complex(0, 3) * Sqrt(2) * q2z * q3y * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) - 3 * Sqrt(2) * q2x * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) -
                Complex(0, 3) * Sqrt(2) * q2y * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner)) *
               Yket(0, 0, 2, 2, thetaq, wigner)) /
              Sqrt(30);
        return mtx;
    }
    else if (idxbra == 16 && idxket == 11)
    {
        mtx = (Ftpe3 * (2 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - 2 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 2) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(2) * q2z * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(3) * q2z * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Sqrt(3) * q2y * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2z * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Sqrt(3) * q2x * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, -1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * q2z * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * q2x * q3z * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * q2z * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 2) * q2y * q3z * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2y * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2x * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2y * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2x * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(6) * q2z * q3x * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Sqrt(6) * q2z * q3y * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(6) * q2x * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Sqrt(6) * q2y * q3z * q4x * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) -
                        Sqrt(6) * q2z * q3x * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(6) * q2z * q3y * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(6) * q2x * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(6) * q2y * q3z * q4y * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 0, 2, 2, thetaq, wigner) +
                        Sqrt(3) * (Complex(0, -1) * q2z * q3x - q2z * q3y + Complex(0, 1) * q2x * q3z + q2y * q3z) * Ybra(2, -2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * (2 * q4z * Yket(1, -1, 2, 2, thetaq, wigner) - Sqrt(2) * (q4x - Complex(0, 1) * q4y) * Yket(1, 0, 2, 2, thetaq, wigner)) +
                        2 * Sqrt(2) * q2y * q3x * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * Sqrt(2) * q2x * q3y * q4x * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(2) * q2y * q3x * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 2) * Sqrt(2) * q2x * q3y * q4y * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(2) * q2z * q3x * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(2) * q2z * q3y * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(2) * q2x * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Sqrt(2) * q2y * q3z * q4z * Ybra(2, 0, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 1) * Sqrt(3) * q2z * q3x * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Sqrt(3) * q2z * q3y * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 1) * Sqrt(3) * q2x * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Sqrt(3) * q2y * q3z * q4x * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(3) * q2z * q3x * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + Complex(0, 1) * Sqrt(3) * q2z * q3y * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Sqrt(3) * q2x * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - Complex(0, 1) * Sqrt(3) * q2y * q3z * q4y * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        2 * Sqrt(3) * q2y * q3x * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * Sqrt(3) * q2x * q3y * q4z * Ybra(2, 1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Complex(0, 2) * Sqrt(3) * q2z * q3x * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) - 2 * Sqrt(3) * q2z * q3y * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) -
                        Complex(0, 2) * Sqrt(3) * q2x * q3z * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) + 2 * Sqrt(3) * q2y * q3z * q4z * Ybra(2, 2, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) * Yket(1, 1, 2, 2, thetaq, wigner) +
                        Sqrt(3) * Ybra(2, -1, 2, 2, thetapp, phipp, thetaqp, phiqp, wigner) *
                            ((q2z * (Complex(0, 1) * q3x + q3y) * (q4x + Complex(0, 1) * q4y) - q2y * (q3z * q4x + Complex(0, 1) * q3z * q4y + 2 * q3x * q4z) + q2x * (Complex(0, -1) * q3z * q4x + q3z * q4y + 2 * q3y * q4z)) * Yket(1, -1, 2, 2, thetaq, wigner) +
                             (q4x - Complex(0, 1) * q4y) * (Sqrt(2) * (q2y * q3x - q2x * q3y) * Yket(1, 0, 2, 2, thetaq, wigner) + (q2z * (Complex(0, 1) * q3x + q3y) - Complex(0, 1) * q2x * q3z - q2y * q3z) * Yket(1, 1, 2, 2, thetaq, wigner))))) /
              Sqrt(30);
        return mtx;
    }

    else
    {
        // std::cerr << "unkown channel index: (" << idxbra << "," << idxket << ") !" << std::endl;
        return 0;
    }
}
} // end namespace aPWD3_c4
