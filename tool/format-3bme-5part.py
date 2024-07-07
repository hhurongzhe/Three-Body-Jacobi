import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

hbarc = 197.32698

# LECs
c1 = -0.81 * 1e-3
c3 = -3.2 * 1e-3
c4 = 5.4 * 1e-3
cD = 1.264
cE = -0.120
lambda_chi = 700
regulator_power = 4
lambda_3nf = 2 * hbarc

# mesh number
number_pmesh = 25
number_qmesh = 25

hbarc5 = hbarc**5
twopi6 = (2.0 * np.pi) ** 6

shape = (number_pmesh, number_qmesh, number_pmesh, number_qmesh)
file_name_c1 = "./data/kernel-c1-bra1-ket1-nmesh12.bin"
file_name_c3 = "./data/kernel-c3-bra1-ket1-nmesh12.bin"
file_name_c4 = "./data/kernel-c4-bra1-ket2-nmesh12.bin"
file_name_cD = "./data/kernel-cD-bra1-ket1-nmesh12.bin"
file_name_cE = "./data/kernel-cE-bra1-ket1-nmesh12.bin"
file_name_pmesh = "./data/kernel-c1-pmesh25.bin"
file_name_qmesh = "./data/kernel-c1-qmesh25.bin"


def read_binary_file(file_name, shape):
    with open(file_name, "rb") as f:
        data = np.fromfile(f, dtype=np.float64)
    return data.reshape(shape)


mtx_c1 = c1 * read_binary_file(file_name_c1, shape) * hbarc5 / twopi6
mtx_c3 = c3 * read_binary_file(file_name_c3, shape) * hbarc5 / twopi6
mtx_c4 = c4 * read_binary_file(file_name_c4, shape) * hbarc5 / twopi6
mtx_cD = cD / lambda_chi * read_binary_file(file_name_cD, shape) * hbarc5 / twopi6
mtx_cE = cE / lambda_chi * read_binary_file(file_name_cE, shape) * hbarc5 / twopi6

mom_mesh_p = read_binary_file(file_name_pmesh, number_pmesh)
mom_mesh_q = read_binary_file(file_name_qmesh, number_qmesh)


def freg(p, q):
    dom = ((p**2 + 0.75 * q**2) / (lambda_3nf**2)) ** regulator_power
    return np.exp(-dom)


# interaction matrix with regulation.
mtx_vc1_reg = np.zeros_like(mtx_c1)
mtx_vc3_reg = np.zeros_like(mtx_c3)
mtx_vc4_reg = np.zeros_like(mtx_c4)
mtx_vcD_reg = np.zeros_like(mtx_cD)
mtx_vcE_reg = np.zeros_like(mtx_cE)

for idxp in range(number_pmesh):
    for idxq in range(number_qmesh):
        p = mom_mesh_p[idxp]
        q = mom_mesh_q[idxq]
        regulator_pq = freg(p, q)
        for idxpp in range(number_pmesh):
            for idxqp in range(number_qmesh):
                pp = mom_mesh_p[idxpp]
                qp = mom_mesh_q[idxqp]
                regulator_ppqp = freg(pp, qp)
                mtx_vc1_reg[idxp, idxq, idxpp, idxqp] = (
                    mtx_c1[idxp, idxq, idxpp, idxqp] * regulator_pq * regulator_ppqp
                )
                mtx_vc3_reg[idxp, idxq, idxpp, idxqp] = (
                    mtx_c3[idxp, idxq, idxpp, idxqp] * regulator_pq * regulator_ppqp
                )
                mtx_vc4_reg[idxp, idxq, idxpp, idxqp] = (
                    mtx_c4[idxp, idxq, idxpp, idxqp] * regulator_pq * regulator_ppqp
                )
                mtx_vcD_reg[idxp, idxq, idxpp, idxqp] = (
                    mtx_cD[idxp, idxq, idxpp, idxqp] * regulator_pq * regulator_ppqp
                )
                mtx_vcE_reg[idxp, idxq, idxpp, idxqp] = (
                    mtx_cE[idxp, idxq, idxpp, idxqp] * regulator_pq * regulator_ppqp
                )

np.save("vc1.npy", mtx_vc1_reg)
np.save("vc3.npy", mtx_vc3_reg)
np.save("vc4.npy", mtx_vc4_reg)
np.save("vcD.npy", mtx_vcD_reg)
np.save("vcE.npy", mtx_vcE_reg)
