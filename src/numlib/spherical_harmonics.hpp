#pragma once
#ifndef SPHERICAL_HARMONICS_HPP
#define SPHERICAL_HARMONICS_HPP

#include <cmath>
#include <cstring>
#include <complex>
#include <iostream>
#include <string>
#include <cstdlib>
#include "WignerSymbol.hpp"
#include "./../constants.hpp"

// realization of spherical harmonics function.
namespace spherical_harmonics
{
    using namespace constants;

    // n! function.
    double factorial(int n)
    {
        if (n < 0)
        {
            std::cerr << "Error: Negative number for factorial." << std::endl;
            exit(EXIT_FAILURE);
        }
        return (n == 1 || n == 0) ? 1.0 : factorial(n - 1) * n;
    }

    // sign function, returning double values.
    double sign_double(int n)
    {
        if (n % 2 == 0)
        {
            return 1.0;
        }
        else
        {
            return -1.0;
        }
    }

    // associate legendre functions P(l,m,x).
    // std::assoc_legendre only support positive m.
    // TODO: check the Condon-Shortley phase term (-1)^m is omitted from this definition?
    double asso_poly(const int l, const int m, const double x)
    {
        if (m > 0)
        {
            const double plm = std::assoc_legendre(l, m, x);
            return plm;
        }
        else
        {
            const double phase = sign_double(-m) * factorial(l + m) / factorial(l - m);
            const double plm = std::assoc_legendre(l, -m, x);
            return phase * plm;
        }
    }

    std::complex<double> Ylm(const int l, const int m, const double theta, const double phi)
    {
        const double fourpi_sqrt_inv = 0.28209479177387814;
        const double prefac = sqrt(2 * l + 1) * fourpi_sqrt_inv * sqrt(factorial(l - m) / factorial(l + m));
        const double plm = asso_poly(l, m, cos(theta));
        std::complex<double> img(cos(m * phi), sin(m * phi));
        return prefac * plm * img;
    }

    // coupled spherical harmonics,
    // see equ.(14) in "A new way to perform partial wave decompositions of few-nucleon forces".
    std::complex<double> Y_coupled(const int L, const int mL, const int l, const int lam, const double thetap, const double phip, const double thetaq, const double phiq, util::WignerSymbols &wigner)
    {
        std::complex<double> temp(0, 0);
        for (int ml = -l; ml <= l; ml = ml + 1)
        {
            if ((abs(ml) > l) || (abs(mL - ml) > lam))
            {
                continue;
            }
            double cg = wigner.CG(2 * l, 2 * lam, 2 * L, 2 * ml, 2 * (mL - ml), 2 * mL);
            std::complex<double> yp = Ylm(l, ml, thetap, phip);
            std::complex<double> yq = Ylm(lam, mL - ml, thetaq, phiq);
            temp = temp + cg * yp * yq;
        }
        return temp;
    }

    // functions used for aPWD,
    // see equ.(13) in "A new way to perform partial wave decompositions of few-nucleon forces".
    std::complex<double>
    Yppqpstar(int L, int mL, int l, int lam, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp, util::WignerSymbols &wigner)
    {
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];
        return std::conj(Y_coupled(L, mL, l, lam, thetapp, phipp, thetaqp, phiqp, wigner));
    }

    std::complex<double> Ypq(int L, int mL, int l, int lam, int idx_thetaq, util::WignerSymbols &wigner)
    {
        double thetap = 0;
        double thetaq = mesh_theta[idx_thetaq];
        double phip = 0;
        double phiq = 0;
        return Y_coupled(L, mL, l, lam, thetap, phip, thetaq, phiq, wigner);
    }

} // end namespace spherical_harmonics

#endif // SPHERICAL_HARMONICS_HPP