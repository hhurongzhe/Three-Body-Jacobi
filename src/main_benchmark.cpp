#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <vector>
#include <chrono>
#include <cmath>

#include "constants.hpp"
#include "aPWD3_part_c1.hpp"
#include "aPWD3_part_c3.hpp"
#include "aPWD3_part_c4.hpp"
#include "numlib/WignerSymbol.hpp"

using namespace constants;

// example LECs.
constexpr double c1 = -0.81 * 1e-3;
constexpr double c3 = -3.4 * 1e-3;
constexpr double c4 = 3.4 * 1e-3;

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

double G(int idx_beta_bra, int idx_beta_ket, double ppmag, double qpmag, double pmag, double qmag, std::vector<double> &angular_mesh_weights, util::WignerSymbols &wigner)
{
    int idx_theta_q, idx_theta_pp, idx_phi_pp, idx_theta_qp, idx_phi_qp;
    double theta_q, theta_pp, theta_qp, phi_pp, phi_qp, w, gt_c1, gt_c3, gt_c4;
    double temp = 0;

#pragma omp parallel for private(idx_theta_q, idx_theta_pp, idx_phi_pp, idx_theta_qp, idx_phi_qp, theta_q, theta_pp, theta_qp, phi_pp, phi_qp, w, gt_c1, gt_c3, gt_c4) reduction(+ : temp) schedule(static)
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
        gt_c1 = real(aPWD3_c1::Gt(idx_beta_bra, idx_beta_ket, ppmag, qpmag, pmag, qmag, theta_q, theta_pp, phi_pp, theta_qp, phi_qp, wigner));
        gt_c3 = real(aPWD3_c3::Gt(idx_beta_bra, idx_beta_ket, ppmag, qpmag, pmag, qmag, theta_q, theta_pp, phi_pp, theta_qp, phi_qp, wigner));
        gt_c4 = real(aPWD3_c4::Gt(idx_beta_bra, idx_beta_ket, ppmag, qpmag, pmag, qmag, theta_q, theta_pp, phi_pp, theta_qp, phi_qp, wigner));
        temp = temp + w * (c1 * gt_c1 + c3 * gt_c3 + c4 * gt_c4);
    }
    return temp;
}

//* for benchmark.
int main()
{
    std::cout << "----you are benchmarking for (c1,c3,c4) !\n\n";
    const int thread_number = omp_get_num_procs();
    omp_set_num_threads(thread_number);
    util::WignerSymbols wigner;
    wigner.reserve(200, "Jmax", 6);
    store_angle_mesh_weights();
    std::ifstream fp_angular_mesh("./data/kernel-angular-mesh-weights.bin", std::ios::binary);
    if (!fp_angular_mesh)
    {
        std::cerr << "Failed to open binary file containing angular mesh weights: ./data/kernel-angular-mesh-weights.bin" << std::endl;
        std::exit(-1);
    }
    std::vector<double> angular_mesh_weights(angle_dimension_total);
    fp_angular_mesh.read(reinterpret_cast<char *>(angular_mesh_weights.data()), angle_dimension_total * sizeof(double));
    fp_angular_mesh.close();

    double ppmag = 3.0 * constants::hbarc;
    double qpmag = 4.0 * constants::hbarc;
    double pmag = 1.0 * constants::hbarc;
    double qmag = 2.0 * constants::hbarc;

    std::vector<std::vector<int>> idx_braket = {{1, 1}, {2, 1}, {6, 11}, {5, 10}};
    std::vector<double> test_values = {443.618, 1200.219, -5.49290, -5.48626};

    for (int idx_channel = 0; idx_channel < idx_braket.size(); idx_channel = idx_channel + 1)
    {
        auto this_channel_index = idx_braket[idx_channel];
        int idx_bra = this_channel_index[0];
        int idx_ket = this_channel_index[1];
        double mtx = G(idx_bra, idx_ket, ppmag, qpmag, pmag, qmag, angular_mesh_weights, wigner);
        const double unit_factor = pow(hbarc, 5);
        double mtx_with_unit = mtx * unit_factor;
        std::cout << "G(" << idx_bra << "," << idx_ket << ")  calculate: " << mtx_with_unit << "    standard: " << test_values[idx_channel] << "\n";
    }
}
