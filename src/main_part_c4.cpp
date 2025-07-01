#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <vector>
#include <chrono>
#include <cmath>

#include "constants.hpp"
#include "aPWD3_part_c4.hpp"
#include "numlib/WignerSymbol.hpp"

using namespace constants;

constexpr int N_channels = 16;  // number of total beta-channels to calculate.
constexpr bool verbose = false; // if print detailed information.

constexpr int angle_dimension_total = Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle;
constexpr int angle_dimension_4 = Nmesh_angle * Nmesh_angle * Nmesh_angle * Nmesh_angle;
constexpr int angle_dimension_3 = Nmesh_angle * Nmesh_angle * Nmesh_angle;
constexpr int angle_dimension_2 = Nmesh_angle * Nmesh_angle;

void store_angle_mesh_weights()
{
    std::cout << "---- calculate and store angular mesh weights..." << std::endl;
    std::ostringstream oss;
    oss << "./data/kernel-angular-mesh-weights.bin"; //---- angular mesh weights file position.
    auto file_name = oss.str();
    std::ofstream fp(file_name, std::ios::binary);
    if (fp)
    {
        int idx_theta_q, idx_theta_pp, idx_phi_pp, idx_theta_qp, idx_phi_qp;
        double theta_q, theta_pp, theta_qp, sin_theta_q, sin_theta_pp, sin_theta_qp, w;
        for (unsigned long long int idx_angle = 0; idx_angle < angle_dimension_total; idx_angle = idx_angle + 1)
        {
            idx_theta_q = (idx_angle / angle_dimension_4);
            idx_theta_pp = (idx_angle / angle_dimension_3) % Nmesh_angle;
            idx_phi_pp = (idx_angle / angle_dimension_2) % Nmesh_angle;
            idx_theta_qp = (idx_angle / Nmesh_angle) % Nmesh_angle;
            idx_phi_qp = idx_angle % Nmesh_angle;
            theta_q = mesh_theta[idx_theta_q];
            theta_pp = mesh_theta[idx_theta_pp];
            theta_qp = mesh_theta[idx_theta_qp];
            sin_theta_q = sin(theta_q);
            sin_theta_pp = sin(theta_pp);
            sin_theta_qp = sin(theta_qp);
            //---- w: mesh weight.
            w = weight_theta[idx_theta_q] * weight_theta[idx_theta_pp] * weight_theta[idx_theta_qp] * weight_phi[idx_phi_pp] * weight_phi[idx_phi_qp] * sin_theta_q * sin_theta_pp * sin_theta_qp * 8 * PI * PI;
            fp.write(reinterpret_cast<const char *>(&w), sizeof(double));
        }
        std::cout << "angular mesh weights stored in: " << file_name << "\n" << std::endl;
    }
    else
    {
        std::cerr << "Failed to create binary file: " << file_name << std::endl;
        std::exit(-1);
    }
}

//---- 3N interaction matrix elements under Jacobi partial-wave basis in LS-coupling scheme:
//---- < beta', p', q' | V3N | beta, p, q >
double G(int idx_beta_bra, int idx_beta_ket, int idxpp, int idxqp, int idxp, int idxq, std::vector<double> &angular_mesh_weights, util::WignerSymbols &wigner)
{
    double ppmag = mesh_mom_p[idxpp];
    double qpmag = mesh_mom_q[idxqp];
    double pmag = mesh_mom_p[idxp];
    double qmag = mesh_mom_q[idxq];
    int idx_theta_q, idx_theta_pp, idx_phi_pp, idx_theta_qp, idx_phi_qp;
    double theta_q, theta_pp, theta_qp, phi_pp, phi_qp, w, gt;
    double temp = 0;

#pragma omp parallel for private(idx_theta_q, idx_theta_pp, idx_phi_pp, idx_theta_qp, idx_phi_qp, theta_q, theta_pp, theta_qp, phi_pp, phi_qp, w, gt) reduction(+ : temp) schedule(static)
    for (unsigned long long int idx_angle = 0; idx_angle < angle_dimension_total; idx_angle = idx_angle + 1)
    {
        //---- use a single loop to do this 5-dim integration, making it more efficient for openmp.
        idx_theta_q = (idx_angle / angle_dimension_4);
        idx_theta_pp = (idx_angle / angle_dimension_3) % Nmesh_angle;
        idx_phi_pp = (idx_angle / angle_dimension_2) % Nmesh_angle;
        idx_theta_qp = (idx_angle / Nmesh_angle) % Nmesh_angle;
        idx_phi_qp = idx_angle % Nmesh_angle;
        theta_q = mesh_theta[idx_theta_q];
        theta_pp = mesh_theta[idx_theta_pp];
        theta_qp = mesh_theta[idx_theta_qp];
        phi_pp = mesh_phi[idx_phi_pp];
        phi_qp = mesh_phi[idx_phi_qp];
        //---- w: mesh weight.
        w = angular_mesh_weights[idx_angle];
        //---- gt: integrated value.
        gt = real(aPWD3_c4::Gt(idx_beta_bra, idx_beta_ket, ppmag, qpmag, pmag, qmag, theta_q, theta_pp, phi_pp, theta_qp, phi_qp, wigner));
        temp = temp + w * gt;
    }
    return temp;
}

void store_mesh()
{
    std::cout << "---- store momentum mesh..." << std::endl;

    std::ostringstream oss_p;
    oss_p << "./data/kernel-c4-pmesh" << Nmesh_mom_p << ".bin";
    auto file_name_p = oss_p.str();
    std::ofstream fp_p(file_name_p, std::ios::binary);
    if (fp_p)
    {
        for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
        {
            double p_value = constants::mesh_mom_p[idxp];
            fp_p.write(reinterpret_cast<const char *>(&p_value), sizeof(double));
        }
        std::cout << "p momentum mesh stored in: " << file_name_p << std::endl;
    }
    else
    {
        std::cerr << "Failed to create binary file: " << file_name_p << std::endl;
        std::exit(-1);
    }

    std::ostringstream oss_q;
    oss_q << "./data/kernel-c4-qmesh" << Nmesh_mom_q << ".bin";
    auto file_name_q = oss_q.str();
    std::ofstream fp_q(file_name_q, std::ios::binary);
    if (fp_q)
    {
        for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
        {
            double q_value = constants::mesh_mom_q[idxq];
            fp_q.write(reinterpret_cast<const char *>(&q_value), sizeof(double));
        }
        std::cout << "q momentum mesh stored in: " << file_name_q << "\n" << std::endl;
    }
    else
    {
        std::cerr << "Failed to create binary file: " << file_name_q << std::endl;
        std::exit(-1);
    }
}

//---- calculate and store on the momentum mesh.
int main()
{
    std::cout << "---- running apwd3...\n\n";

    //---- print current date.
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    char dateStr[100];
    std::strftime(dateStr, sizeof(dateStr), "%Y-%m-%d", std::localtime(&now_c));
    std::cout << "---- Current Date: " << dateStr << "\n" << std::endl;

    std::cout << "---- you are calculating for part c4 !\n" << std::endl;

    //---- estimate memory.
    double memory_min_each = constants::Nmesh_mom_p * constants::Nmesh_mom_p * constants::Nmesh_mom_q * constants::Nmesh_mom_q * 8. / 1024. / 1024.;
    double memory_min_all = N_channels * N_channels * memory_min_each;
    std::cout << "---- each file  size: " << memory_min_each << " MB" << std::endl;
    std::cout << "---- all  files size: " << memory_min_all << " MB\n" << std::endl;

    //---- set parallel threads in openpm.
    const int thread_number = omp_get_num_procs();
    std::cout << "---- number of threads for openmp: " << thread_number << "\n" << std::endl;
    omp_set_num_threads(thread_number);

    //---- reserve wigner symbols.
    util::WignerSymbols wigner;
    wigner.reserve(200, "Jmax", 6);
    std::cout << "---- reserve wigner symbols !\n" << std::endl;

    //---- store momentum mesh.
    store_mesh();

    //---- calculate and store angular mesh weights.
    store_angle_mesh_weights();
    //---- open the file containing precomputed weights.
    std::ifstream fp_angular_mesh("./data/kernel-angular-mesh-weights.bin", std::ios::binary);
    if (!fp_angular_mesh)
    {
        std::cerr << "Failed to open binary file containing angular mesh weights: ./data/kernel-angular-mesh-weights.bin" << std::endl;
        std::exit(-1);
    }
    //---- Read weights into a vector.
    std::vector<double> angular_mesh_weights(angle_dimension_total);
    fp_angular_mesh.read(reinterpret_cast<char *>(angular_mesh_weights.data()), angle_dimension_total * sizeof(double));
    fp_angular_mesh.close();

    auto time1 = std::chrono::high_resolution_clock::now();
    //---- index of channels: 1,2,...,N_channels
    for (int idx_channel_bra = 1; idx_channel_bra <= N_channels; idx_channel_bra = idx_channel_bra + 1)
    {
        for (int idx_channel_ket = 1; idx_channel_ket <= N_channels; idx_channel_ket = idx_channel_ket + 1)
        {
            std::cout << "---- calculating for channel : (beta',beta) = (" << idx_channel_bra << "," << idx_channel_ket << ")" << std::endl;
            auto time_part1 = std::chrono::high_resolution_clock::now();
            std::ostringstream oss;
            oss << "./data/kernel-c4-bra" << idx_channel_bra << "-ket" << idx_channel_ket << "-nmesh" << constants::Nmesh_angle << ".bin";
            auto file_name = oss.str();
            std::ofstream fp(file_name, std::ios::binary);
            if (fp)
            {
                bool if_estimate_timing = false; // a flag to denote if already estimate the timing in this channel.
                for (int idxp = 0; idxp < Nmesh_mom_p; idxp = idxp + 1)
                {
                    for (int idxq = 0; idxq < Nmesh_mom_q; idxq = idxq + 1)
                    {
                        for (int idxpp = 0; idxpp < Nmesh_mom_p; idxpp = idxpp + 1)
                        {
                            auto time_estimate_1 = std::chrono::high_resolution_clock::now();
                            for (int idxqp = 0; idxqp < Nmesh_mom_q; idxqp = idxqp + 1)
                            {
                                double mtx = G(idx_channel_bra, idx_channel_ket, idxpp, idxqp, idxp, idxq, angular_mesh_weights, wigner);
                                if (verbose)
                                {
                                    std::cout << "(p',q',p,q) = (" << idxpp << "," << idxqp << "," << idxp << "," << idxq << ")   mtx:  " << mtx << std::endl;
                                }
                                fp.write(reinterpret_cast<const char *>(&mtx), sizeof(double));
                            }
                            auto time_estimate_2 = std::chrono::high_resolution_clock::now();
                            auto duration_estimate = std::chrono::duration_cast<std::chrono::milliseconds>(time_estimate_2 - time_estimate_1).count();
                            if (!if_estimate_timing)
                            {
                                std::cout << "---- This Channel Timing (estimate) : " << 0.001 * duration_estimate * Nmesh_mom_p * Nmesh_mom_q * Nmesh_mom_p << " seconds" << std::endl;
                                if_estimate_timing = true;
                            }
                        }
                    }
                }
            }
            else
            {
                std::cerr << "Failed to create binary file: " << file_name << std::endl;
                std::exit(-1);
            }
            auto time_part2 = std::chrono::high_resolution_clock::now();
            auto duration_part = std::chrono::duration_cast<std::chrono::seconds>(time_part2 - time_part1).count();
            std::cout << "---- This Channel Timing (actual)   : " << duration_part << " seconds\n" << std::endl;
        }
    }
    auto time2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
    std::cout << "---- All Channel Timing : " << duration << " seconds" << std::endl;

    std::cout << "---- Terminate sucessfully !" << std::endl;
    return 1;
}
