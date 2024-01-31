#pragma once
#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <iostream>
#include <vector>

namespace constants
{
    // 3N LECs:
    // ci in MeV^(-1),
    // cD, cE dimensionless,
    // Lambda, LambdaX in MeV.
    constexpr double c1 = -0.81 * 1e-3;
    constexpr double c3 = -3.2 * 1e-3;
    constexpr double c4 = 5.4 * 1e-3;
    constexpr double cD = 1.264;
    constexpr double cE = -0.120;
    constexpr int regulator_power = 2;
    constexpr double Lambda = 500;
    constexpr double LambdaX = 700.0;

    // partial-wave channel and cuts, only used when building up 3N partial-wave channels.
    constexpr int TWOJ = 1;
    constexpr int PARITY = 1;
    constexpr int TWOT = 1;
    constexpr int j_MAX = 5; // max total momentum of the unpaired nucleon.

    // interval of 5-dim angular index, used for effective reduction of memory.
    // you should combine all results to get the final right result.
    constexpr unsigned long long int left_angular_index = 0;      // TODO: check if smaller than right and if positive.
    constexpr unsigned long long int right_angular_index = 10000; // TODO: check if smaller than biggest possible value.

    // result file name.
    const std::string FILE_PREFIX = "./data/test_2J_" + std::to_string(TWOJ) + "_P_" + std::to_string(PARITY) + "_2T_" + std::to_string(TWOT) + "_index_" + std::to_string(left_angular_index) + "_" + std::to_string(right_angular_index);
    const std::string FILE_SUFFIX = ".dat";
    const std::string FILE_NAME = FILE_PREFIX + FILE_SUFFIX;

    // Constants.
    constexpr double gA = 1.29;
    constexpr double fpi = 92.4;
    constexpr double Mpi = 138.0390;
    constexpr double hbarc = 197.32698;
    constexpr double PI = 3.141592653589793;
    constexpr int64_t NUM_THREADS = 128;

    // momentum mesh points and weights for Jacobi momentum p, p' and q, q in (MeV).
    constexpr int Nmesh_mom_p = 4;
    constexpr double mesh_mom_p[Nmesh_mom_p] = {100, 200, 300, 400};
    constexpr double weight_mom_p[Nmesh_mom_p] = {100, 200, 300, 400};
    constexpr int Nmesh_mom_q = 4;
    constexpr double mesh_mom_q[Nmesh_mom_q] = {100, 200, 300, 400};
    constexpr double weight_mom_q[Nmesh_mom_q] = {100, 200, 300, 400};
    // constexpr int Nmesh_mom_p = 16;
    // constexpr double mesh_mom_p[Nmesh_mom_p] = {2.6497662520875167, 13.856244231691846, 33.59219940304206, 61.148897911249264, 95.53093889933905, 135.49580558569318, 179.59911230518526, 226.24687254059063, 273.75312745940937, 320.40088769481474, 364.5041944143068, 404.46906110066095, 438.85110208875074, 466.40780059695794, 486.14375576830815, 497.3502337479125};
    // constexpr double weight_mom_p[Nmesh_mom_p] = {6.788114852938509, 15.563380984661926, 23.789627920623147, 31.157242813883506, 37.398997204144194, 42.289129848750655, 45.650853761230906, 47.362652613767146, 47.362652613767146, 45.650853761230906, 42.289129848750655, 37.398997204144194, 31.157242813883506, 23.789627920623147, 15.563380984661926, 6.788114852938509};
    // constexpr int Nmesh_mom_q = 16;
    // constexpr double mesh_mom_q[Nmesh_mom_q] = {2.6497662520875167, 13.856244231691846, 33.59219940304206, 61.148897911249264, 95.53093889933905, 135.49580558569318, 179.59911230518526, 226.24687254059063, 273.75312745940937, 320.40088769481474, 364.5041944143068, 404.46906110066095, 438.85110208875074, 466.40780059695794, 486.14375576830815, 497.3502337479125};
    // constexpr double weight_mom_q[Nmesh_mom_q] = {6.788114852938509, 15.563380984661926, 23.789627920623147, 31.157242813883506, 37.398997204144194, 42.289129848750655, 45.650853761230906, 47.362652613767146, 47.362652613767146, 45.650853761230906, 42.289129848750655, 37.398997204144194, 31.157242813883506, 23.789627920623147, 15.563380984661926, 6.788114852938509};

    // mesh points and weights for angle integration.
    constexpr int Nmesh_angle = 16;
    constexpr double mesh_theta[Nmesh_angle] = {0.016648972382576677, 0.08706135016925809, 0.21106601372504086, 0.3842098569061858, 0.600238591673398, 0.851345254840489, 1.1284545036184366, 1.4215510253423718, 1.7200416282474214, 2.0131381499713568, 2.290247398749304, 2.5413540619163952, 2.7573827966836073, 2.9305266398647523, 3.054531303420535, 3.1249436812072164};
    constexpr double weight_theta[Nmesh_angle] = {0.04265098350743076, 0.09778760673286598, 0.1494746406141286, 0.1957667302604196, 0.2349848297363292, 0.2657104393190798, 0.28683277361277, 0.29758832301187255, 0.29758832301187255, 0.28683277361277, 0.2657104393190798, 0.2349848297363292, 0.1957667302604196, 0.1494746406141286, 0.09778760673286598, 0.04265098350743076};
    constexpr double mesh_phi[Nmesh_angle] = {0.03329794476515335, 0.17412270033851618, 0.4221320274500817, 0.7684197138123716, 1.200477183346796, 1.702690509680978, 2.256909007236873, 2.8431020506847435, 3.4400832564948427, 4.0262762999427135, 4.580494797498608, 5.0827081238327905, 5.514765593367215, 5.8610532797295045, 6.10906260684107, 6.249887362414433};
    constexpr double weight_phi[Nmesh_angle] = {0.08530196701486152, 0.19557521346573195, 0.2989492812282572, 0.3915334605208392, 0.4699696594726584, 0.5314208786381596, 0.57366554722554, 0.5951766460237451, 0.5951766460237451, 0.57366554722554, 0.5314208786381596, 0.4699696594726584, 0.3915334605208392, 0.2989492812282572, 0.19557521346573195, 0.08530196701486152};

} // end namespace constants

#endif // CONSTANTS_HPP