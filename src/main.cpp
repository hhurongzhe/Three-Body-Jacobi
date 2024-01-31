#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>
#include <omp.h>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>

#include "constants.hpp"
#include "aPWD3.hpp"
#include "numlib/WignerSymbol.hpp"

using namespace constants;

// 3N interaction matrix elements under Jacobi partial-wave basis in LS-coupling scheme.
// < beta', p', q' | V3N | beta, p, q >
double G(int lp, int lamp, int Lp, int sp, int twoSp, int tp, int l, int lam, int L, int s, int twoS, int t, int twoT, int twoJ,
         int idxp, int idxq, int idxpp, int idxqp,
         precalculate::integration_weight_Container &integration_weight_container,
         precalculate::Ybra_Container &Ybra_container, precalculate::Yket_Container &Yket_container,
         precalculate::F_TPE1_12_Container &F_TPE1_12_container, precalculate::F_TPE1_13_Container &F_TPE1_13_container, precalculate::F_TPE1_23_Container &F_TPE1_23_container,
         precalculate::F_TPE2_12_Container &F_TPE2_12_container, precalculate::F_TPE2_13_Container &F_TPE2_13_container, precalculate::F_TPE2_23_Container &F_TPE2_23_container,
         precalculate::F_OPE_1_Container &F_OPE_1_container, precalculate::F_OPE_2_Container &F_OPE_2_container, precalculate::F_OPE_3_Container &F_OPE_3_container)
{
    double temp = 0;
    int idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp;
    double gg, weight;
    // openmp is used for angular integration only.
#pragma omp parallel for private(idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp, gg, weight) reduction(+ : temp) schedule(dynamic)
    for (unsigned long long int idx_angle = left_angular_index; idx_angle < right_angular_index; idx_angle = idx_angle + 1)
    {
        // use a single loop to do this 5-dim integration, making it simplier for omp.
        idx_thetaq = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle));
        idx_thetapp = (idx_angle / (Nmesh_angle * Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
        idx_phipp = (idx_angle / (Nmesh_angle * Nmesh_angle)) % Nmesh_angle;
        idx_thetaqp = (idx_angle / Nmesh_angle) % Nmesh_angle;
        idx_phiqp = idx_angle % Nmesh_angle;
        gg = aPWD3::GG(idx_angle, idx_thetaq, idx_thetapp, idx_phipp, idx_thetaqp, idx_phiqp,
                       lp, lamp, Lp, sp, twoSp, tp, l, lam, L, s, twoS, t, twoT, twoJ,
                       idxp, idxq, idxpp, idxqp,
                       Ybra_container, Yket_container,
                       F_TPE1_12_container, F_TPE1_13_container, F_TPE1_23_container,
                       F_TPE2_12_container, F_TPE2_13_container, F_TPE2_23_container,
                       F_OPE_1_container, F_OPE_2_container, F_OPE_3_container);
        weight = integration_weight_container.get_value(idx_angle);
        temp = temp + gg * weight;
    }
    return temp;
}

// 3N interaction matrix elements under Jacobi partial-wave basis in JJ-coupling scheme.
// < alpha', p', q' | V3N | alpha, p, q >
double H(int lp, int sp, int jp, int lamp, int twoj3p, int tp, int l, int s, int j, int lam, int twoj3, int t, int twoT, int twoJ, int idxp, int idxq, int idxpp, int idxqp,
         precalculate::integration_weight_Container &integration_weight_container,
         precalculate::Ybra_Container &Ybra_container, precalculate::Yket_Container &Yket_container,
         precalculate::F_TPE1_12_Container &F_TPE1_12_container, precalculate::F_TPE1_13_Container &F_TPE1_13_container, precalculate::F_TPE1_23_Container &F_TPE1_23_container,
         precalculate::F_TPE2_12_Container &F_TPE2_12_container, precalculate::F_TPE2_13_Container &F_TPE2_13_container, precalculate::F_TPE2_23_Container &F_TPE2_23_container,
         precalculate::F_OPE_1_Container &F_OPE_1_container, precalculate::F_OPE_2_Container &F_OPE_2_container, precalculate::F_OPE_3_Container &F_OPE_3_container,
         util::WignerSymbols &wigner)
{
    double temp = 0;
    int Lp_min, Lp_max, L_min, L_max;
    double fac_bra, fac_ket, f9j_bra, f9j_ket, factor, mtx;
    // don't check angular momentum coupling in 9j symbol, because there is not too much cases.
    for (int twoSp = 1; twoSp <= 3; twoSp = twoSp + 2)
    {
        if ((sp == 0) && (twoSp = 3))
        {
            continue;
        }
        Lp_min = abs((twoSp - twoJ) / 2);
        Lp_max = abs((twoSp + twoJ) / 2);
        for (int Lp = Lp_min; Lp <= Lp_max; Lp = Lp + 1)
        {
            for (int twoS = 1; twoS <= 3; twoS = twoS + 2)
            {
                if ((s == 0) && (twoS = 3))
                {
                    continue;
                }
                L_min = abs((twoS - twoJ) / 2);
                L_max = abs((twoS + twoJ) / 2);
                for (int L = L_min; L <= L_max; L = L + 1)
                {
                    // standard transformation from LS coupling to JJ coupling,
                    // by 9j symbol and some "hat" prefactors.
                    fac_bra = sqrt(sqrt((2 * Lp + 1) * (twoSp + 1) * (2 * jp + 1) * (twoj3p + 1)));
                    fac_ket = sqrt(sqrt((2 * L + 1) * (twoS + 1) * (2 * j + 1) * (twoj3 + 1)));
                    f9j_bra = wigner.f9j(2 * lp, 2 * sp, 2 * jp, 2 * lamp, 1, twoj3p, 2 * Lp, twoSp, twoJ);
                    f9j_ket = wigner.f9j(2 * l, 2 * s, 2 * j, 2 * lam, 1, twoj3, 2 * L, twoS, twoJ);
                    factor = fac_bra * fac_ket * f9j_bra * f9j_ket;
                    mtx = G(lp, lamp, Lp, sp, twoSp, tp, l, lam, L, s, twoS, t, twoT, twoJ,
                            idxp, idxq, idxpp, idxqp, integration_weight_container,
                            Ybra_container, Yket_container,
                            F_TPE1_12_container, F_TPE1_13_container, F_TPE1_23_container,
                            F_TPE2_12_container, F_TPE2_13_container, F_TPE2_23_container,
                            F_OPE_1_container, F_OPE_2_container, F_OPE_3_container);
                    temp = temp + factor * mtx;
                }
            }
        }
    }

    return temp;
}

int sign(int n)
{
    if (n % 2 == 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

// struct of 3N partial-wave channel in JJ-coupling scheme.
struct alpha_channel
{
    int index;  // channel index.
    int l;      // relative angular momentum within pair (23).
    int s;      // spin of pair (23).
    int j;      // total momentum of pair (23).
    int lambda; // relative angular momentum between pair (23) and (1).
    int twoj1;  // twice of total momentum of (1).
    int twoJ;   // twice of total J.
    int twoMJ;  // twice of J projection, never used.
    int t;      // isospin of pair (23).
    int twoT;   // twice of total isospin.
    int twoMT;  // twice of T projection, never used.
    int P;      // parity.

    // Constructor for initialization
    alpha_channel(int i, int l_value, int s_value, int j_value, int lambda_value, int twoj1_value, int twoJ_value, int twoMJ_value,
                  int t_value, int twoT_value, int twoMT_value)
        : index(i), l(l_value), s(s_value), j(j_value), lambda(lambda_value), twoj1(twoj1_value), twoJ(twoJ_value), twoMJ(twoMJ_value),
          t(t_value), twoT(twoT_value), twoMT(twoMT_value)
    {
        P = sign(l_value + lambda_value);
    }
};

// build partial-wave 3N channels for fixed J, P and T.
std::vector<alpha_channel> setup_alpha_channels(int twoJ, int P, int twoT)
{
    bool check = (twoJ > 0) && (P == 1 || P == -1) && (twoT == 1 || twoT == 3);
    if (!check)
    {
        throw std::invalid_argument("Invalid quantum number for J, P, T !");
    }
    std::vector<alpha_channel> alpha_channels;
    int idx = 0;
    for (int j = 0; j <= j_MAX; j = j + 1)
    {
        for (int s = 0; s <= 1; s = s + 1)
        {
            int lmin = abs(j - s);
            int lmax = abs(j + s);
            for (int l = lmin; l <= lmax; l = l + 1)
            {
                int twoj3min = abs(twoJ - 2 * j);
                int twoj3max = abs(twoJ + 2 * j);
                for (int twoj3 = twoj3min; twoj3 <= twoj3max; twoj3 = twoj3 + 2)
                {
                    int lammin = abs(twoj3 - 1) / 2;
                    int lammax = abs(twoj3 + 1) / 2;
                    for (int lambda = lammin; lambda <= lammax; lambda = lambda + 1)
                    {
                        if (P * sign(l + lambda) < 0)
                        {
                            continue;
                        }
                        int t;
                        if ((l + s) % 2 == 0)
                        {
                            t = 1;
                        }
                        else
                        {
                            t = 0;
                        }
                        if ((twoT == 3) && (t == 0))
                        {
                            continue;
                        }
                        int twoMJ = twoJ;
                        int twoMT = twoT;
                        alpha_channel chan(idx, l, s, j, lambda, twoj3, twoJ, twoMJ, t, twoT, twoMT);
                        alpha_channels.push_back(chan);
                        idx = idx + 1;
                    }
                }
            }
        }
    }
    return alpha_channels;
}

// calculate 3bme in specific alpha channel, for all momentum mesh points.
void cal_3bme_alpha_channel(alpha_channel chanel_bra, alpha_channel chanel_ket,
                            precalculate::integration_weight_Container &integration_weight_container,
                            precalculate::Ybra_Container &Ybra_container, precalculate::Yket_Container &Yket_container,
                            precalculate::F_TPE1_12_Container &F_TPE1_12_container, precalculate::F_TPE1_13_Container &F_TPE1_13_container, precalculate::F_TPE1_23_Container &F_TPE1_23_container,
                            precalculate::F_TPE2_12_Container &F_TPE2_12_container, precalculate::F_TPE2_13_Container &F_TPE2_13_container, precalculate::F_TPE2_23_Container &F_TPE2_23_container,
                            precalculate::F_OPE_1_Container &F_OPE_1_container, precalculate::F_OPE_2_Container &F_OPE_2_container, precalculate::F_OPE_3_Container &F_OPE_3_container,
                            util::WignerSymbols &wigner)
{
    std::ofstream fp(FILE_NAME, std::ios::app); // append mode.
    if (!fp.is_open())
    {
        std::cerr << "failed to open file: " << FILE_NAME << '\n';
        std::exit(-1);
    }
    fp << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;
    bool check_J = (chanel_bra.twoJ == chanel_ket.twoJ);
    bool check_P = (chanel_bra.P == chanel_ket.P);
    bool check_T = (chanel_bra.twoT == chanel_ket.twoT);
    if (check_J && check_P && check_T)
    {
        int lp = chanel_bra.l;
        int sp = chanel_bra.s;
        int jp = chanel_bra.j;
        int lamp = chanel_bra.lambda;
        int twoj1p = chanel_bra.twoj1;
        int tp = chanel_bra.t;
        int twoTp = chanel_bra.twoT;

        int l = chanel_ket.l;
        int s = chanel_ket.s;
        int j = chanel_ket.j;
        int lam = chanel_ket.lambda;
        int twoj1 = chanel_ket.twoj1;
        int t = chanel_ket.t;
        int twoT = chanel_ket.twoT;

        int twoJ = chanel_bra.twoJ;

        fp << "# (l',s',j',lambda',2j1',t',l,s,j,lambda,2j1,t,2J,P,2T) = (" << lp << "," << sp << "," << jp << "," << lamp << "," << twoj1p << "," << tp << "," << l << "," << s << "," << j << "," << lam << "," << twoj1 << "," << t << "," << twoJ << "," << PARITY << "," << twoT << ")" << std::endl;
        fp << "##" << std::endl;
        for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
        {
            for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
            {
                for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                {
                    for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                    {
                        double mtx = H(lp, sp, jp, lamp, twoj1p, tp, l, s, j, lam, twoj1, t, twoT, twoJ, idxp, idxq, idxpp, idxqp,
                                       integration_weight_container,
                                       Ybra_container, Yket_container,
                                       F_TPE1_12_container, F_TPE1_13_container, F_TPE1_23_container,
                                       F_TPE2_12_container, F_TPE2_13_container, F_TPE2_23_container,
                                       F_OPE_1_container, F_OPE_2_container, F_OPE_3_container,
                                       wigner);
                        const double unit_factor = pow(hbarc, 5);
                        double mtx_with_unit = mtx * unit_factor;
                        if (std::abs(mtx_with_unit) > 1e-16)
                        {
                            fp << "(p',q',p,q) = (" << idxpp << "," << idxqp << "," << idxp << "," << idxq << ")    ";
                            fp << "mtx: " << std::scientific << std::setprecision(17) << mtx_with_unit << std::endl;
                        }
                        // if (std::abs(mtx * pow(hbarc, 5)) > 1e-16)
                        // {
                        //     std::cout << "(p',q',p,q) = (" << idxpp << "," << idxqp << "," << idxp << "," << idxq << ")    ";
                        //     std::cout << "mtx: " << mtx * pow(hbarc, 5) << " fm^(5)" << std::endl;
                        // }
                        // std::cout << "(p',q',p,q) = (" << idxpp << "," << idxqp << "," << idxp << "," << idxq << ")    ";
                        // std::cout << "mtx: " << mtx * pow(hbarc, 5) << " fm^(5)" << std::endl;
                    }
                }
            }
        }
    }
    else
    {
        std::cout << "quantum number violation between alpha channels.    check_J : " << check_J << "    check_P : " << check_P << "    check_T : " << check_T << std::endl;
    }
    fp << "##" << std::endl;
    fp << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;
    fp.close();
}

// writting header into the file, including all parameters.
void write_headers()
{
    std::ofstream fp(FILE_NAME);
    if (!fp.is_open())
    {
        std::cerr << "failed to open file: " << FILE_NAME << '\n';
        std::exit(-1);
    }

    fp << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;
    fp << "! Three-Body-Jacobi: C++ code to calculate matrix elements of chiral 3N interactions in Jacobi momentum space" << std::endl;

    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    char dateStr[100];
    std::strftime(dateStr, sizeof(dateStr), "%Y-%m-%d", std::localtime(&now_c));
    fp << "\n! Current Date: " << dateStr << std::endl;

    // Write constants to the file
    fp << "\n! 3N LECs:" << std::endl;
    fp << "# c1 (MeV^(-1)): " << c1 << std::endl;
    fp << "# c3 (MeV^(-1)): " << c3 << std::endl;
    fp << "# c4 (MeV^(-1)): " << c4 << std::endl;
    fp << "# cD (dimensionless): " << cD << std::endl;
    fp << "# cE (dimensionless): " << cE << std::endl;
    fp << "# regulator_power: " << regulator_power << std::endl;
    fp << "# Lambda (MeV): " << Lambda << std::endl;
    fp << "# LambdaX (MeV): " << LambdaX << std::endl;

    fp << "\n! Partial-Wave Channel and Cuts:" << std::endl;
    fp << "# TWOJ: " << TWOJ << std::endl;
    fp << "# PARITY: " << PARITY << std::endl;
    fp << "# TWOT: " << TWOT << std::endl;
    fp << "# j_MAX: " << j_MAX << std::endl;

    fp << "\n! 5-Dim Angular Index Interval:" << std::endl;
    fp << "# full left_angular_index: " << 0 << std::endl;
    fp << "# full right_angular_index: " << Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle - 1 << std::endl;
    fp << "# scaling left_angular_index: " << left_angular_index << std::endl;
    fp << "# scaling right_angular_index: " << right_angular_index << std::endl;

    fp << "\n! Constants:" << std::endl;
    fp << "# gA: " << gA << std::endl;
    fp << "# fpi: " << fpi << std::endl;
    fp << "# Mpi: " << Mpi << std::endl;
    fp << "# hbarc: " << hbarc << std::endl;
    fp << "# PI: " << PI << std::endl;
    fp << "# NUM_THREADS: " << NUM_THREADS << std::endl;

    fp << "\n! mesh_theta, weight_theta, mesh_phi, weight_phi:" << std::endl;
    fp << "# Nmesh_angle: " << Nmesh_angle << std::endl;
    fp << "##" << std::endl;
    for (int i = 0; i < Nmesh_angle; ++i)
    {
        fp << std::fixed << "\t" << std::scientific << std::setprecision(17) << mesh_theta[i] << "\t" << std::setw(17) << weight_theta[i] << "\t" << std::setw(17) << mesh_phi[i] << "\t" << std::setw(17) << weight_phi[i] << std::endl;
    }
    fp << "##" << std::endl;

    fp << "\n! mesh_mom_p, weight_mom_p:" << std::endl;
    fp << "# Nmesh_mom_p: " << Nmesh_mom_p << std::endl;
    fp << "##" << std::endl;
    for (int i = 0; i < Nmesh_mom_p; ++i)
    {
        fp << std::fixed << "\t" << std::scientific << std::setprecision(17) << mesh_mom_p[i] << "\t" << std::setw(17) << weight_mom_p[i] << std::endl;
    }
    fp << "##" << std::endl;

    fp << "\n! mesh_mom_q, weight_mom_q:" << std::endl;
    fp << "# Nmesh_mom_q: " << Nmesh_mom_q << std::endl;
    fp << "##" << std::endl;
    for (int i = 0; i < Nmesh_mom_q; ++i)
    {
        fp << std::fixed << "\t" << std::scientific << std::setprecision(17) << mesh_mom_q[i] << "\t" << std::setw(17) << weight_mom_q[i] << std::endl;
    }
    fp << "##" << std::endl;

    fp << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n"
       << std::endl;
    fp.close();
    std::cout << "\nWriting header complete !\n"
              << std::endl;
}

int main()
{
    std::cout << "running 3b-jacobi.exe..." << std::endl;
    write_headers();

    // set parallel threads in openpm.
    omp_set_num_threads(NUM_THREADS);

    auto time1 = std::chrono::high_resolution_clock::now();
    auto alpha_channels = setup_alpha_channels(TWOJ, PARITY, TWOT);
    auto time2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    int Nalpha = alpha_channels.size();
    std::cout << "Building alpha channels !    Timing : " << duration << " seconds" << std::endl;

    time1 = std::chrono::high_resolution_clock::now();
    util::WignerSymbols wigner;
    wigner.reserve(200, "Jmax", 9);
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "\nReserve 9j symbols !    Timing : " << duration << " seconds" << std::endl;

    unsigned long long int angle_dimensions = pow(Nmesh_angle, 5);
    unsigned long long int momentum_dimensions = pow(Nmesh_mom_p, 2) * pow(Nmesh_mom_q, 2);
    unsigned long long int tbme_dimensions = angle_dimensions * momentum_dimensions;
    double precalculate_storage = 9 * tbme_dimensions * 8 / pow(1024, 3); // precalculate functions' minimal required storage (in GB).
    std::cout << "\nBeginning precalculate functions......\n";
    std::cout << "Without scaling, angular index interval   : "
              << "[ 0 , " << angle_dimensions - 1 << " ]\n";
    std::cout << "After   scaling, angular index interval   : "
              << "[ " << left_angular_index
              << " , " << right_angular_index - 1 << " ]\n";
    std::cout << "Without scaling, 5-dim angular dimension  : " << angle_dimensions << "\n";
    std::cout << "After   scaling, 5-dim angular dimension  : " << right_angular_index - left_angular_index + 1 << "\n";
    std::cout << "Without scaling, minimal storage required : " << precalculate_storage << " GB\n";
    std::cout << "After   scaling, minimal storage required : " << precalculate_storage * (right_angular_index - left_angular_index) / angle_dimensions << " GB\n"
              << std::endl;

    precalculate::integration_weight_Container integration_weight_container;

    precalculate::Ybra_Container Ybra_container;
    precalculate::Yket_Container Yket_container;

    precalculate::F_TPE1_12_Container F_TPE1_12_container;
    precalculate::F_TPE1_13_Container F_TPE1_13_container;
    precalculate::F_TPE1_23_Container F_TPE1_23_container;

    precalculate::F_TPE2_12_Container F_TPE2_12_container;
    precalculate::F_TPE2_13_Container F_TPE2_13_container;
    precalculate::F_TPE2_23_Container F_TPE2_23_container;

    precalculate::F_OPE_1_Container F_OPE_1_container;
    precalculate::F_OPE_2_Container F_OPE_2_container;
    precalculate::F_OPE_3_Container F_OPE_3_container;

    time1 = std::chrono::high_resolution_clock::now();
    integration_weight_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "5-dim angular integration weights calculated and stored!    Timing : " << duration << " seconds" << std::endl;

    time1 = std::chrono::high_resolution_clock::now();
    F_TPE1_12_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "F_TPE1_12 function calculated and stored!    Timing : " << duration << " seconds" << std::endl;
    time1 = std::chrono::high_resolution_clock::now();
    F_TPE1_13_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "F_TPE1_13 function calculated and stored!    Timing : " << duration << " seconds" << std::endl;
    time1 = std::chrono::high_resolution_clock::now();
    F_TPE1_23_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "F_TPE1_23 function calculated and stored!    Timing : " << duration << " seconds" << std::endl;

    time1 = std::chrono::high_resolution_clock::now();
    F_TPE2_12_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "F_TPE2_12 function calculated and stored!    Timing : " << duration << " seconds" << std::endl;
    time1 = std::chrono::high_resolution_clock::now();
    F_TPE2_13_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "F_TPE2_13 function calculated and stored!    Timing : " << duration << " seconds" << std::endl;
    time1 = std::chrono::high_resolution_clock::now();
    F_TPE2_23_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "F_TPE2_23 function calculated and stored!    Timing : " << duration << " seconds" << std::endl;

    time1 = std::chrono::high_resolution_clock::now();
    F_OPE_1_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "F_OPE_1 function calculated and stored!    Timing : " << duration << " seconds" << std::endl;
    time1 = std::chrono::high_resolution_clock::now();
    F_OPE_2_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "F_OPE_2 function calculated and stored!    Timing : " << duration << " seconds" << std::endl;
    time1 = std::chrono::high_resolution_clock::now();
    F_OPE_3_container.store_all_value();
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "F_OPE_3 function calculated and stored!    Timing : " << duration << " seconds" << std::endl;

    time1 = std::chrono::high_resolution_clock::now();
    Ybra_container.store_all_value(wigner);
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "Ybra function calculated and stored!    Timing : " << duration << " seconds" << std::endl;
    time1 = std::chrono::high_resolution_clock::now();
    Yket_container.store_all_value(wigner);
    time2 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "Yket function calculated and stored!    Timing : " << duration << " seconds" << std::endl;

    std::cout << "\nPrecalculation all done!\n\nBeginning calculating 3bmes......\n"
              << std::endl;

    std::cout << "Calculating 3bmes in channel space     : J = " << TWOJ << "/2"
              << " , P = " << PARITY << " , T = " << TWOT << "/2" << std::endl;
    std::cout << "Basis dimension Nalpha in this channel : Nalpha = " << Nalpha << std::endl;
    std::cout << "Dimension of V3N in this channel       : dim(V3N) = " << Nalpha * Nalpha * momentum_dimensions << std::endl;
    int dummy = 1;
    for (int idxbra = 0; idxbra < Nalpha; idxbra = idxbra + 1)
    {
        for (int idxket = 0; idxket < Nalpha; idxket = idxket + 1)
        {
            auto chanel_bra = alpha_channels[idxbra];
            auto chanel_ket = alpha_channels[idxket];
            time1 = std::chrono::high_resolution_clock::now();
            cal_3bme_alpha_channel(chanel_bra, chanel_ket,
                                   integration_weight_container,
                                   Ybra_container, Yket_container,
                                   F_TPE1_12_container,
                                   F_TPE1_13_container, F_TPE1_23_container,
                                   F_TPE2_12_container, F_TPE2_13_container, F_TPE2_23_container,
                                   F_OPE_1_container, F_OPE_2_container, F_OPE_3_container,
                                   wigner);
            time2 = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
            std::cout << dummy << "/" << Nalpha * Nalpha << "    Timing : " << duration << " seconds" << std::endl;
            dummy = dummy + 1;
        }
    }
    if (dummy == Nalpha * Nalpha + 1)
    {
        std::cout << "\nTerminate successfully !" << std::endl;
    }
    else
    {
        std::cout << "\nTerminate ahead !" << std::endl;
    }
}
