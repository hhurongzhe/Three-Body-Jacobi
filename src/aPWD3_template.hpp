#pragma once
#ifndef aPWD3_HPP
#define aPWD3_HPP

#include <iostream>
#include <tuple>
#include <vector>
#include <cmath>

#include "constants.hpp"
#include "precalculate.hpp"
#include "numlib/WignerSymbol.hpp"
#include "numlib/spherical_harmonics.hpp"

// to treat "integer \times complex".
std::complex<double> operator*(const int &c1, const std::complex<double> &c2)
{
    return (c1 * 1.0) * c2;
}

namespace aPWD3
{
    using std::conj;
    using std::exp;
    using std::imag;
    using std::pow;
    using std::real;
    using namespace constants;

    // nonlocal regulator function.
    double get_regulator_nonlocal(double pmag, double qmag)
    {
        int n = regulator_power;
        double regulator = exp(-pow((pmag * pmag + 0.75 * qmag * qmag) / (Lambda * Lambda), n));
        return regulator;
    }

    // a simple complex function.
    std::complex<double> Complex(double x, double y)
    {
        std::complex<double> temp(x, y);
        return temp;
    }

    // Partial-wave decomposition of 3N forces.
    // In fact it should be a real number.
    // However, due to calculation of spherical harmonics and other things numerically, it returns a complex number.
    // If handling truncations and precisions appropriatly, it's imaginary part should be very small. So we keep it for check.
    // You can just take it's real part.
    double GG(int idx_angle, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp,
              int lp, int lamp, int Lp, int sp, int twoSp, int tp, int l, int lam, int L, int s, int twoS, int t, int twoT, int twoJ,
              int idxp, int idxq, int idxpp, int idxqp,
              precalculate::Ybra_Container &Ybracontainer, precalculate::Yket_Container &Yketcontainer,
              precalculate::F_TPE1_12_Container &F_TPE1_12_container, precalculate::F_TPE1_13_Container &F_TPE1_13_container, precalculate::F_TPE1_23_Container &F_TPE1_23_container,
              precalculate::F_TPE2_12_Container &F_TPE2_12_container, precalculate::F_TPE2_13_Container &F_TPE2_13_container, precalculate::F_TPE2_23_Container &F_TPE2_23_container,
              precalculate::F_OPE_1_Container &F_OPE_1_container, precalculate::F_OPE_2_Container &F_OPE_2_container, precalculate::F_OPE_3_Container &F_OPE_3_container)
    {
        std::tuple<int, int, int, int, int, int, int, int, int, int, int, int, int, int> spin_tuple(lp, lamp, Lp, sp, twoSp, tp, l, lam, L, s, twoS, t, twoT, twoJ);

        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_q[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        // single momentum transfer of nucleon 1,2,3.
        // in principle we can store these values to speed up, but memory will be doubled at least.
        double q1x = precalculate::get_q1x(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double q1y = precalculate::get_q1y(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double q1z = precalculate::get_q1z(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double q2x = precalculate::get_q2x(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double q2y = precalculate::get_q2y(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double q2z = precalculate::get_q2z(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double q3x = precalculate::get_q3x(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double q3y = precalculate::get_q3y(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double q3z = precalculate::get_q3z(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);

        // two-pion exchange functions.
        // (c1, c3, c4) dependent.
        double Ftpeone12 = F_TPE1_12_container.get_value(idxp, idxq, idxpp, idxqp, idx_angle);
        double Ftpeone13 = F_TPE1_13_container.get_value(idxp, idxq, idxpp, idxqp, idx_angle);
        double Ftpeone23 = F_TPE1_23_container.get_value(idxp, idxq, idxpp, idxqp, idx_angle);
        double Ftpetwo12 = F_TPE2_12_container.get_value(idxp, idxq, idxpp, idxqp, idx_angle);
        double Ftpetwo13 = F_TPE2_13_container.get_value(idxp, idxq, idxpp, idxqp, idx_angle);
        double Ftpetwo23 = F_TPE2_23_container.get_value(idxp, idxq, idxpp, idxqp, idx_angle);

        // one-pion exchange functions.
        // (cD) dependent.
        double Fope1 = F_OPE_1_container.get_value(idxp, idxq, idxpp, idxqp, idx_angle);
        double Fope2 = F_OPE_2_container.get_value(idxp, idxq, idxpp, idxqp, idx_angle);
        double Fope3 = F_OPE_3_container.get_value(idxp, idxq, idxpp, idxqp, idx_angle);

        // contact functions.
        // (cE) dependent.
        double Fcontact = precalculate::get_F_contact();

        // 3bme in specific partial-wave channels.
        // automated partial-wave projection (aPWD) is done by Mathematica,
        // which generates codes below.
        // I only provide a template in the lowest channel.
        std::complex<double> mtx = -((6 * Fcontact + Fope2 * (pow(q2x, 2) + pow(q2y, 2) + pow(q2z, 2)) + 2 * Ftpeone23 * q2x * q3x + Fope3 * pow(q3x, 2) + 2 * Ftpeone23 * q2y * q3y + Fope3 * pow(q3y, 2) + 2 * Ftpeone23 * q2z * q3z + Fope3 * pow(q3z, 2)) * Ybracontainer.get_value(0, 0, 0, 0, idx_angle) * Yketcontainer.get_value(0, 0, 0, 0, idx_angle));
        return real(mtx) * get_regulator_nonlocal(pmag, qmag) * get_regulator_nonlocal(ppmag, qpmag);
        if (spin_tuple == std::make_tuple(0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1))
        {
            std::complex<double> mtx = -((6 * Fcontact + Fope2 * (pow(q2x, 2) + pow(q2y, 2) + pow(q2z, 2)) + 2 * Ftpeone23 * q2x * q3x + Fope3 * pow(q3x, 2) + 2 * Ftpeone23 * q2y * q3y + Fope3 * pow(q3y, 2) + 2 * Ftpeone23 * q2z * q3z + Fope3 * pow(q3z, 2)) * Ybracontainer.get_value(0, 0, 0, 0, idx_angle) * Yketcontainer.get_value(0, 0, 0, 0, idx_angle));
            return real(mtx) * get_regulator_nonlocal(pmag, qmag) * get_regulator_nonlocal(ppmag, qpmag);
        }
        else
        {
            std::complex<double> mtx = 0.0;
            return real(mtx) * get_regulator_nonlocal(pmag, qmag) * get_regulator_nonlocal(ppmag, qpmag);
        }
    }
} // end namespace aPWD3

#endif // aPWD3_HPP