#pragma once
#ifndef PRECALCULATE_HPP
#define PRECALCULATE_HPP

#include <iostream>
#include <tuple>
#include <vector>
#include <cmath>
#include <unordered_map>

#include "constants.hpp"
#include "numlib/WignerSymbol.hpp"
#include "numlib/spherical_harmonics.hpp"

// in this namespace there are some functions and stuctures,
// used for precalculating and storing spin-independent functions.
// this is necessary to speed up the code, with the cost of big memory requirement.
namespace precalculate
{
    using namespace constants;
    using std::cos;
    using std::pow;
    using std::sin;

    // weight for 5-dim angular integration.
    double get_integration_weight(int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double w_thetaq = weight_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double w_thetapp = weight_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double w_thetaqp = weight_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double w_phipp = weight_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];
        double w_phiqp = weight_phi[idx_phiqp];
        return w_thetaq * w_thetapp * w_thetaqp * w_phipp * w_phiqp * sin(thetaq) * sin(thetapp) * sin(thetaqp) * 8 * PI * PI; //* maybe a missing 1/(2Pi)^6 factor.
    }

    struct integration_weight_Container
    {
    private:
        double *container;

    public:
        integration_weight_Container()
        {
            container = new double[right_angular_index - left_angular_index];
        }

        void store_all_value()
        {
            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
            {
                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                int idx_phiqp = idx_angle % Nmesh_angle;
                double weight = get_integration_weight(idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                container[idx_angle - left_angular_index] = weight;
            }
        }

        double get_value(int idx_angle)
        {
            return container[idx_angle - left_angular_index];
        }

        ~integration_weight_Container()
        {
            delete[] container;
        }
    };

    double get_q1x(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q1x = -(qmag * sin(thetaq)) + qpmag * cos(phiqp) * sin(thetaqp);
        return q1x;
    }

    double get_q1y(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q1y = qpmag * sin(phiqp) * sin(thetaqp);
        return q1y;
    }

    double get_q1z(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q1z = -(qmag * cos(thetaq)) + qpmag * cos(thetaqp);
        return q1z;
    }

    double get_q2x(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q2x = ppmag * cos(phipp) * sin(thetapp) + (qmag * sin(thetaq)) / 2. - (qpmag * cos(phiqp) * sin(thetaqp)) / 2.;
        return q2x;
    }

    double get_q2y(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q2y = ppmag * sin(phipp) * sin(thetapp) - (qpmag * sin(phiqp) * sin(thetaqp)) / 2.;
        return q2y;
    }

    double get_q2z(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q2z = -pmag + ppmag * cos(thetapp) + (qmag * cos(thetaq)) / 2. - (qpmag * cos(thetaqp)) / 2.;
        return q2z;
    }

    double get_q3x(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q3x = -(ppmag * cos(phipp) * sin(thetapp)) + (qmag * sin(thetaq)) / 2. - (qpmag * cos(phiqp) * sin(thetaqp)) / 2.;
        return q3x;
    }

    double get_q3y(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q3y = -(ppmag * sin(phipp) * sin(thetapp)) - (qpmag * sin(phiqp) * sin(thetaqp)) / 2.;
        return q3y;
    }

    double get_q3z(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q3z = pmag - ppmag * cos(thetapp) + (qmag * cos(thetaq)) / 2. - (qpmag * cos(thetaqp)) / 2.;
        return q3z;
    }

    double get_q1q1(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q1q1 = pow(qmag, 2) + pow(qpmag, 2) - 2 * qmag * qpmag * cos(thetaq) * cos(thetaqp) -
                      2 * qmag * qpmag * cos(phiqp) * sin(thetaq) * sin(thetaqp);
        return q1q1;
    }

    double get_q2q2(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q2q2 = (pow(
                           2 * pmag - 2 * ppmag * cos(thetapp) - qmag * cos(thetaq) + qpmag * cos(thetaqp),
                           2) +
                       pow(
                           2 * ppmag * cos(phipp) * sin(thetapp) + qmag * sin(thetaq) - qpmag * cos(phiqp) * sin(thetaqp),
                           2) +
                       pow(
                           -2 * ppmag * sin(phipp) * sin(thetapp) + qpmag * sin(phiqp) * sin(thetaqp),
                           2)) /
                      4.0;
        return q2q2;
    }

    double get_q3q3(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q3q3 = (pow(
                           2 * pmag - 2 * ppmag * cos(thetapp) + qmag * cos(thetaq) - qpmag * cos(thetaqp),
                           2) +
                       pow(
                           2 * ppmag * cos(phipp) * sin(thetapp) - qmag * sin(thetaq) + qpmag * cos(phiqp) * sin(thetaqp),
                           2) +
                       pow(
                           2 * ppmag * sin(phipp) * sin(thetapp) + qpmag * sin(phiqp) * sin(thetaqp),
                           2)) /
                      4.0;
        return q3q3;
    }

    double get_q1q2(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q1q2 = -0.5 * pow(qmag, 2) - pow(qpmag, 2) / 2. - qpmag * (pmag - ppmag * cos(thetapp)) * cos(thetaqp) +
                      qmag * cos(thetaq) * (pmag - ppmag * cos(thetapp) + qpmag * cos(thetaqp)) -
                      ppmag * qmag * cos(phipp) * sin(thetapp) * sin(thetaq) +
                      ppmag * qpmag * cos(phipp) * cos(phiqp) * sin(thetapp) * sin(thetaqp) +
                      ppmag * qpmag * sin(phipp) * sin(phiqp) * sin(thetapp) * sin(thetaqp) +
                      qmag * qpmag * cos(phiqp) * sin(thetaq) * sin(thetaqp);
        return q1q2;
    }

    double get_q1q3(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q1q3 = -0.5 * pow(qmag, 2) - pow(qpmag, 2) / 2. + qpmag * (pmag - ppmag * cos(thetapp)) * cos(thetaqp) +
                      qmag * cos(thetaq) * (-pmag + ppmag * cos(thetapp) + qpmag * cos(thetaqp)) +
                      ppmag * qmag * cos(phipp) * sin(thetapp) * sin(thetaq) -
                      ppmag * qpmag * cos(phipp) * cos(phiqp) * sin(thetapp) * sin(thetaqp) -
                      ppmag * qpmag * sin(phipp) * sin(phiqp) * sin(thetapp) * sin(thetaqp) +
                      qmag * qpmag * cos(phiqp) * sin(thetaq) * sin(thetaqp);
        return q1q3;
    }

    double get_q2q3(double pmag, double qmag, double ppmag, double qpmag, double thetaq, double thetapp, double phipp, double thetaqp, double phiqp)
    {
        double q2q3 = (-4 * pow(pmag, 2) - 4 * pow(ppmag, 2) + pow(qmag, 2) + pow(qpmag, 2) + 8 * pmag * ppmag * cos(thetapp) - 2 * qmag * qpmag * cos(thetaq) * cos(thetaqp) - 2 * qmag * qpmag * cos(phiqp) * sin(thetaq) * sin(thetaqp)) / 4.0;
        return q2q3;
    }

    double get_F_TPE1_12(int idxp, int idxq, int idxpp, int idxqp, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_p[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        double qiqi = get_q1q1(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double qjqj = get_q2q2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double qiqj = get_q1q2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double part1 = 0.5 * pow((gA / (2 * fpi)), 2);
        double part2 = 1 / (qiqi + Mpi * Mpi) / (qjqj + Mpi * Mpi);
        double part3 = -4 * c1 * Mpi * Mpi / (fpi * fpi) + 2 * c3 / (fpi * fpi) * qiqj;
        double F_TPE1_ij = part1 * part2 * part3;
        return F_TPE1_ij;
    }

    double get_F_TPE1_13(int idxp, int idxq, int idxpp, int idxqp, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_p[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        double qiqi = get_q1q1(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double qjqj = get_q3q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double qiqj = get_q1q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double part1 = 0.5 * pow((gA / (2 * fpi)), 2);
        double part2 = 1 / (qiqi + Mpi * Mpi) / (qjqj + Mpi * Mpi);
        double part3 = -4 * c1 * Mpi * Mpi / (fpi * fpi) + 2 * c3 / (fpi * fpi) * qiqj;
        double F_TPE1_ij = part1 * part2 * part3;
        return F_TPE1_ij;
    }

    double get_F_TPE1_23(int idxp, int idxq, int idxpp, int idxqp, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_p[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        double qiqi = get_q2q2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double qjqj = get_q3q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double qiqj = get_q2q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double part1 = 0.5 * pow((gA / (2 * fpi)), 2);
        double part2 = 1 / (qiqi + Mpi * Mpi) / (qjqj + Mpi * Mpi);
        double part3 = -4 * c1 * Mpi * Mpi / (fpi * fpi) + 2 * c3 / (fpi * fpi) * qiqj;
        double F_TPE1_ij = part1 * part2 * part3;
        return F_TPE1_ij;
    }

    double get_F_TPE2_12(int idxp, int idxq, int idxpp, int idxqp, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_p[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        double qiqi = get_q1q1(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double qjqj = get_q2q2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double part1 = 0.5 * pow((gA / (2 * fpi)), 2);
        double part2 = 1 / (qiqi + Mpi * Mpi) / (qjqj + Mpi * Mpi);
        double part3 = c4 / (fpi * fpi);
        double F_TPE2_ij = part1 * part2 * part3;
        return F_TPE2_ij;
    }

    double get_F_TPE2_13(int idxp, int idxq, int idxpp, int idxqp, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_p[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        double qiqi = get_q1q1(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double qjqj = get_q3q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double part1 = 0.5 * pow((gA / (2 * fpi)), 2);
        double part2 = 1 / (qiqi + Mpi * Mpi) / (qjqj + Mpi * Mpi);
        double part3 = c4 / (fpi * fpi);
        double F_TPE2_ij = part1 * part2 * part3;
        return F_TPE2_ij;
    }

    double get_F_TPE2_23(int idxp, int idxq, int idxpp, int idxqp, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_p[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        double qiqi = get_q2q2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double qjqj = get_q3q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double part1 = 0.5 * pow((gA / (2 * fpi)), 2);
        double part2 = 1 / (qiqi + Mpi * Mpi) / (qjqj + Mpi * Mpi);
        double part3 = c4 / (fpi * fpi);
        double F_TPE2_ij = part1 * part2 * part3;
        return F_TPE2_ij;
    }

    double get_F_OPE_1(int idxp, int idxq, int idxpp, int idxqp, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_p[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        double qjqj = get_q1q1(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double part1 = -gA / (8 * fpi * fpi);
        double part2 = cD / (fpi * fpi * LambdaX);
        double part3 = 1 / (qjqj + Mpi * Mpi);
        double F_OPE_j = part1 * part2 * part3;
        return F_OPE_j;
    }

    double get_F_OPE_2(int idxp, int idxq, int idxpp, int idxqp, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_p[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        double qjqj = get_q2q2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double part1 = -gA / (8 * fpi * fpi);
        double part2 = cD / (fpi * fpi * LambdaX);
        double part3 = 1 / (qjqj + Mpi * Mpi);
        double F_OPE_j = part1 * part2 * part3;
        return F_OPE_j;
    }

    double get_F_OPE_3(int idxp, int idxq, int idxpp, int idxqp, int idx_thetaq, int idx_thetapp, int idx_phipp, int idx_thetaqp, int idx_phiqp)
    {
        double thetaq = mesh_theta[idx_thetaq];
        double thetapp = mesh_theta[idx_thetapp];
        double thetaqp = mesh_theta[idx_thetaqp];
        double phipp = mesh_phi[idx_phipp];
        double phiqp = mesh_phi[idx_phiqp];

        double pmag = mesh_mom_p[idxp];
        double qmag = mesh_mom_q[idxq];
        double ppmag = mesh_mom_p[idxpp];
        double qpmag = mesh_mom_q[idxqp];

        double qjqj = get_q3q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp);
        double part1 = -gA / (8 * fpi * fpi);
        double part2 = cD / (fpi * fpi * LambdaX);
        double part3 = 1 / (qjqj + Mpi * Mpi);
        double F_OPE_j = part1 * part2 * part3;
        return F_OPE_j;
    }

    double get_F_contact()
    {
        double F_contact = cE / (2 * pow(fpi, 4) * LambdaX);
        return F_contact;
    }

    struct F_TPE1_12_Container
    {
    private:
        double *****container;

    public:
        F_TPE1_12_Container()
        {
            container = new double ****[Nmesh_mom_p];
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                container[i] = new double ***[Nmesh_mom_q];
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    container[i][j] = new double **[Nmesh_mom_p];
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        container[i][j][k] = new double *[Nmesh_mom_q];
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            container[i][j][k][l] = new double[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        double get_value(int idxp, int idxq, int idxpp, int idxqp, int idx_angle)
        {
            return container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index];
        }

        void store_all_value()
        {
            for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
            {
                for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                {
                    for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                    {
                        for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                        {
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                double mtx = get_F_TPE1_12(idxp, idxq, idxpp, idxqp, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                                container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index] = mtx;
                            }
                        }
                    }
                }
            }
        }

        ~F_TPE1_12_Container()
        {
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            delete[] container[i][j][k][l];
                        }
                        delete[] container[i][j][k];
                    }
                    delete[] container[i][j];
                }
                delete[] container[i];
            }
            delete[] container;
        }
    };

    struct F_TPE1_13_Container
    {
    private:
        double *****container;

    public:
        F_TPE1_13_Container()
        {
            container = new double ****[Nmesh_mom_p];
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                container[i] = new double ***[Nmesh_mom_q];
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    container[i][j] = new double **[Nmesh_mom_p];
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        container[i][j][k] = new double *[Nmesh_mom_q];
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            container[i][j][k][l] = new double[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        double get_value(int idxp, int idxq, int idxpp, int idxqp, int idx_angle)
        {
            return container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index];
        }

        void store_all_value()
        {
            for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
            {
                for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                {
                    for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                    {
                        for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                        {
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                double mtx = get_F_TPE1_13(idxp, idxq, idxpp, idxqp, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                                container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index] = mtx;
                            }
                        }
                    }
                }
            }
        }

        ~F_TPE1_13_Container()
        {
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            delete[] container[i][j][k][l];
                        }
                        delete[] container[i][j][k];
                    }
                    delete[] container[i][j];
                }
                delete[] container[i];
            }
            delete[] container;
        }
    };

    struct F_TPE1_23_Container
    {
    private:
        double *****container;

    public:
        F_TPE1_23_Container()
        {
            container = new double ****[Nmesh_mom_p];
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                container[i] = new double ***[Nmesh_mom_q];
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    container[i][j] = new double **[Nmesh_mom_p];
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        container[i][j][k] = new double *[Nmesh_mom_q];
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            container[i][j][k][l] = new double[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        double get_value(int idxp, int idxq, int idxpp, int idxqp, int idx_angle)
        {
            return container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index];
        }

        void store_all_value()
        {
            for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
            {
                for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                {
                    for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                    {
                        for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                        {
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                double mtx = get_F_TPE1_23(idxp, idxq, idxpp, idxqp, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                                container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index] = mtx;
                            }
                        }
                    }
                }
            }
        }

        ~F_TPE1_23_Container()
        {
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            delete[] container[i][j][k][l];
                        }
                        delete[] container[i][j][k];
                    }
                    delete[] container[i][j];
                }
                delete[] container[i];
            }
            delete[] container;
        }
    };

    struct F_TPE2_12_Container
    {
    private:
        double *****container;

    public:
        F_TPE2_12_Container()
        {
            container = new double ****[Nmesh_mom_p];
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                container[i] = new double ***[Nmesh_mom_q];
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    container[i][j] = new double **[Nmesh_mom_p];
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        container[i][j][k] = new double *[Nmesh_mom_q];
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            container[i][j][k][l] = new double[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        double get_value(int idxp, int idxq, int idxpp, int idxqp, int idx_angle)
        {
            return container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index];
        }

        void store_all_value()
        {
            for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
            {
                for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                {
                    for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                    {
                        for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                        {
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                double mtx = get_F_TPE2_12(idxp, idxq, idxpp, idxqp, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                                container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index] = mtx;
                            }
                        }
                    }
                }
            }
        }

        ~F_TPE2_12_Container()
        {
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            delete[] container[i][j][k][l];
                        }
                        delete[] container[i][j][k];
                    }
                    delete[] container[i][j];
                }
                delete[] container[i];
            }
            delete[] container;
        }
    };

    struct F_TPE2_13_Container
    {
    private:
        double *****container;

    public:
        F_TPE2_13_Container()
        {
            container = new double ****[Nmesh_mom_p];
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                container[i] = new double ***[Nmesh_mom_q];
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    container[i][j] = new double **[Nmesh_mom_p];
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        container[i][j][k] = new double *[Nmesh_mom_q];
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            container[i][j][k][l] = new double[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        double get_value(int idxp, int idxq, int idxpp, int idxqp, int idx_angle)
        {
            return container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index];
        }

        void store_all_value()
        {
            for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
            {
                for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                {
                    for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                    {
                        for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                        {
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                double mtx = get_F_TPE2_13(idxp, idxq, idxpp, idxqp, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                                container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index] = mtx;
                            }
                        }
                    }
                }
            }
        }

        ~F_TPE2_13_Container()
        {
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            delete[] container[i][j][k][l];
                        }
                        delete[] container[i][j][k];
                    }
                    delete[] container[i][j];
                }
                delete[] container[i];
            }
            delete[] container;
        }
    };

    struct F_TPE2_23_Container
    {
    private:
        double *****container;

    public:
        F_TPE2_23_Container()
        {
            container = new double ****[Nmesh_mom_p];
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                container[i] = new double ***[Nmesh_mom_q];
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    container[i][j] = new double **[Nmesh_mom_p];
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        container[i][j][k] = new double *[Nmesh_mom_q];
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            container[i][j][k][l] = new double[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        double get_value(int idxp, int idxq, int idxpp, int idxqp, int idx_angle)
        {
            return container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index];
        }

        void store_all_value()
        {
            for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
            {
                for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                {
                    for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                    {
                        for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                        {
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                double mtx = get_F_TPE2_23(idxp, idxq, idxpp, idxqp, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                                container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index] = mtx;
                            }
                        }
                    }
                }
            }
        }

        ~F_TPE2_23_Container()
        {
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            delete[] container[i][j][k][l];
                        }
                        delete[] container[i][j][k];
                    }
                    delete[] container[i][j];
                }
                delete[] container[i];
            }
            delete[] container;
        }
    };

    struct F_OPE_1_Container
    {
    private:
        double *****container;

    public:
        F_OPE_1_Container()
        {
            container = new double ****[Nmesh_mom_p];
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                container[i] = new double ***[Nmesh_mom_q];
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    container[i][j] = new double **[Nmesh_mom_p];
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        container[i][j][k] = new double *[Nmesh_mom_q];
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            container[i][j][k][l] = new double[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        double get_value(int idxp, int idxq, int idxpp, int idxqp, int idx_angle)
        {
            return container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index];
        }

        void store_all_value()
        {
            for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
            {
                for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                {
                    for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                    {
                        for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                        {
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                double mtx = get_F_OPE_1(idxp, idxq, idxpp, idxqp, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                                container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index] = mtx;
                            }
                        }
                    }
                }
            }
        }

        ~F_OPE_1_Container()
        {
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            delete[] container[i][j][k][l];
                        }
                        delete[] container[i][j][k];
                    }
                    delete[] container[i][j];
                }
                delete[] container[i];
            }
            delete[] container;
        }
    };

    struct F_OPE_2_Container
    {
    private:
        double *****container;

    public:
        F_OPE_2_Container()
        {
            container = new double ****[Nmesh_mom_p];
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                container[i] = new double ***[Nmesh_mom_q];
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    container[i][j] = new double **[Nmesh_mom_p];
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        container[i][j][k] = new double *[Nmesh_mom_q];
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            container[i][j][k][l] = new double[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        double get_value(int idxp, int idxq, int idxpp, int idxqp, int idx_angle)
        {
            return container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index];
        }

        void store_all_value()
        {
            for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
            {
                for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                {
                    for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                    {
                        for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                        {
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                double mtx = get_F_OPE_2(idxp, idxq, idxpp, idxqp, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                                container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index] = mtx;
                            }
                        }
                    }
                }
            }
        }

        ~F_OPE_2_Container()
        {
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            delete[] container[i][j][k][l];
                        }
                        delete[] container[i][j][k];
                    }
                    delete[] container[i][j];
                }
                delete[] container[i];
            }
            delete[] container;
        }
    };

    struct F_OPE_3_Container
    {
    private:
        double *****container;

    public:
        F_OPE_3_Container()
        {
            container = new double ****[Nmesh_mom_p];
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                container[i] = new double ***[Nmesh_mom_q];
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    container[i][j] = new double **[Nmesh_mom_p];
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        container[i][j][k] = new double *[Nmesh_mom_q];
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            container[i][j][k][l] = new double[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        double get_value(int idxp, int idxq, int idxpp, int idxqp, int idx_angle)
        {
            return container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index];
        }

        void store_all_value()
        {
            for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
            {
                for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                {
                    for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                    {
                        for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                        {
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                double mtx = get_F_OPE_3(idxp, idxq, idxpp, idxqp, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp);
                                container[idxp][idxq][idxpp][idxqp][idx_angle - left_angular_index] = mtx;
                            }
                        }
                    }
                }
            }
        }

        ~F_OPE_3_Container()
        {
            for (int i = 0; i < Nmesh_mom_p; ++i)
            {
                for (int j = 0; j < Nmesh_mom_q; ++j)
                {
                    for (int k = 0; k < Nmesh_mom_p; ++k)
                    {
                        for (int l = 0; l < Nmesh_mom_q; ++l)
                        {
                            delete[] container[i][j][k][l];
                        }
                        delete[] container[i][j][k];
                    }
                    delete[] container[i][j];
                }
                delete[] container[i];
            }
            delete[] container;
        }
    };

    struct Ybra_Container
    {
    private:
        std::complex<double> *****container;

    public:
        Ybra_Container()
        {
            int L_num = (TWOJ + 3) / 2 + 1;
            container = new std::complex<double> ****[L_num];
            for (int idx_L = 0; idx_L < L_num; ++idx_L)
            {
                int L = idx_L;
                int mL_num = 2 * L + 1;
                container[idx_L] = new std::complex<double> ***[mL_num];
                for (int idx_mL = 0; idx_mL < mL_num; ++idx_mL)
                {
                    int l_num = j_MAX + (TWOJ + 1) / 2 + 2; // consider the biggest value.
                    container[idx_L][idx_mL] = new std::complex<double> **[l_num];
                    for (int idx_l = 0; idx_l < l_num; ++idx_l)
                    {
                        int l = idx_l;
                        int lam_num = 2 * l + 1;
                        container[idx_L][idx_mL][idx_l] = new std::complex<double> *[lam_num];
                        for (int idx_lam = 0; idx_lam < lam_num; ++idx_lam)
                        {
                            container[idx_L][idx_mL][idx_l][idx_lam] = new std::complex<double>[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        std::complex<double> get_value(int L, int mL, int l, int lam, int idx_angle)
        {
            int idx_L = L;
            int idx_mL = mL + L;
            int idx_l = l;
            int idx_lam = lam - abs(L - l);
            return container[idx_L][idx_mL][idx_l][idx_lam][idx_angle - left_angular_index];
        }

        void store_all_value(util::WignerSymbols &wigner)
        {
            int L_num = (TWOJ + 3) / 2 + 1;
            for (int idx_L = 0; idx_L < L_num; ++idx_L)
            {
                int L = idx_L;
                int mL_num = 2 * L + 1;
                for (int idx_mL = 0; idx_mL < mL_num; ++idx_mL)
                {
                    int mL = idx_mL - L;
                    int l_num = j_MAX + (TWOJ + 1) / 2 + 2; // consider the biggest value.
                    for (int idx_l = 0; idx_l < l_num; ++idx_l)
                    {
                        int l = idx_l;
                        int lam_num = 2 * l + 1;
                        for (int idx_lam = 0; idx_lam < lam_num; ++idx_lam)
                        {
                            int lam = idx_lam + abs(L - l);
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                std::complex<double> value = spherical_harmonics::Yppqpstar(L, mL, l, lam, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp, wigner);
                                container[idx_L][idx_mL][idx_l][idx_lam][idx_angle - left_angular_index] = value;
                            }
                        }
                    }
                }
            }
        }

        ~Ybra_Container()
        {
            int L_num = (TWOJ + 3) / 2 + 1;
            for (int idx_L = 0; idx_L < L_num; ++idx_L)
            {
                int L = idx_L;
                int mL_num = 2 * L + 1;
                for (int idx_mL = 0; idx_mL < mL_num; ++idx_mL)
                {
                    int l_num = j_MAX + (TWOJ + 1) / 2 + 2; // consider the biggest value.
                    for (int idx_l = 0; idx_l < l_num; ++idx_l)
                    {
                        int l = idx_l;
                        int lam_num = 2 * l + 1;
                        for (int idx_lam = 0; idx_lam < lam_num; ++idx_lam)
                        {
                            delete[] container[idx_L][idx_mL][idx_l][idx_lam];
                        }
                        delete[] container[idx_L][idx_mL][idx_l];
                    }
                    delete[] container[idx_L][idx_mL];
                }
                delete[] container[idx_L];
            }
            delete[] container;
        }
    };

    struct Yket_Container
    {
    private:
        std::complex<double> *****container;

    public:
        Yket_Container()
        {
            int L_num = (TWOJ + 3) / 2 + 1;
            container = new std::complex<double> ****[L_num];
            for (int idx_L = 0; idx_L < L_num; ++idx_L)
            {
                int L = idx_L;
                int mL_num = 2 * L + 1;
                container[idx_L] = new std::complex<double> ***[mL_num];
                for (int idx_mL = 0; idx_mL < mL_num; ++idx_mL)
                {
                    int l_num = j_MAX + (TWOJ + 1) / 2 + 2; // consider the biggest value.
                    container[idx_L][idx_mL] = new std::complex<double> **[l_num];
                    for (int idx_l = 0; idx_l < l_num; ++idx_l)
                    {
                        int l = idx_l;
                        int lam_num = 2 * l + 1;
                        container[idx_L][idx_mL][idx_l] = new std::complex<double> *[lam_num];
                        for (int idx_lam = 0; idx_lam < lam_num; ++idx_lam)
                        {
                            container[idx_L][idx_mL][idx_l][idx_lam] = new std::complex<double>[right_angular_index - left_angular_index];
                        }
                    }
                }
            }
        }

        std::complex<double> get_value(int L, int mL, int l, int lam, int idx_angle)
        {
            int idx_L = L;
            int idx_mL = mL + L;
            int idx_l = l;
            int idx_lam = lam - abs(L - l);
            return container[idx_L][idx_mL][idx_l][idx_lam][idx_angle - left_angular_index];
        }

        void store_all_value(util::WignerSymbols &wigner)
        {
            int L_num = (TWOJ + 3) / 2 + 1;
            for (int idx_L = 0; idx_L < L_num; ++idx_L)
            {
                int L = idx_L;
                int mL_num = 2 * L + 1;
                for (int idx_mL = 0; idx_mL < mL_num; ++idx_mL)
                {
                    int mL = idx_mL - L;
                    int l_num = j_MAX + (TWOJ + 1) / 2 + 2; // consider the biggest value.
                    for (int idx_l = 0; idx_l < l_num; ++idx_l)
                    {
                        int l = idx_l;
                        int lam_num = 2 * l + 1;
                        for (int idx_lam = 0; idx_lam < lam_num; ++idx_lam)
                        {
                            int lam = idx_lam + abs(L - l);
                            for (int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
                            {
                                int idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
                                int idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
                                int idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
                                int idx_phiqp = idx_angle % Nmesh_angle;
                                std::complex<double> value = spherical_harmonics::Ypq(L, mL, l, lam, idx_thetaq, wigner);
                                container[idx_L][idx_mL][idx_l][idx_lam][idx_angle - left_angular_index] = value;
                            }
                        }
                    }
                }
            }
        }

        ~Yket_Container()
        {
            int L_num = (TWOJ + 3) / 2 + 1;
            for (int idx_L = 0; idx_L < L_num; ++idx_L)
            {
                int L = idx_L;
                int mL_num = 2 * L + 1;
                for (int idx_mL = 0; idx_mL < mL_num; ++idx_mL)
                {
                    int l_num = j_MAX + (TWOJ + 1) / 2 + 2; // consider the biggest value.
                    for (int idx_l = 0; idx_l < l_num; ++idx_l)
                    {
                        int l = idx_l;
                        int lam_num = 2 * l + 1;
                        for (int idx_lam = 0; idx_lam < lam_num; ++idx_lam)
                        {
                            delete[] container[idx_L][idx_mL][idx_l][idx_lam];
                        }
                        delete[] container[idx_L][idx_mL][idx_l];
                    }
                    delete[] container[idx_L][idx_mL];
                }
                delete[] container[idx_L];
            }
            delete[] container;
        }
    };

} // end namespace precalculate

#endif // PRECALCULATE_HPP