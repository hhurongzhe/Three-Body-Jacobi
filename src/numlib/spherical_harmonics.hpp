#pragma once
#ifndef SPHERICAL_HARMONICS_HPP
#define SPHERICAL_HARMONICS_HPP

#include <complex>
#include "WignerSymbol.hpp"

// realization of spherical harmonics function.
namespace spherical_harmonics
{
    std::complex<double> Ylm(const int l, const int m, const double theta, const double phi);

    std::complex<double> Ybra(const int l, const int m, const int la, const int lb, const double theta_a, const double phi_a, const double theta_b, const double phi_b, util::WignerSymbols &wigner);

    std::complex<double> Yket(const int l, const int m, const int la, const int lb, const double theta_b, util::WignerSymbols &wigner);

} // end namespace spherical_harmonics

#endif // SPHERICAL_HARMONICS_HPP