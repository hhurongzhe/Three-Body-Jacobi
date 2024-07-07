import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

hbarc = 197.32698

# EM1.8/2.0 LECs
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
file_name_c4 = "./data/kernel-c4-bra1-ket1-nmesh12.bin"
file_name_cD = "./data/kernel-cD-bra1-ket1-nmesh12.bin"
file_name_cE = "./data/kernel-cE-bra1-ket1-nmesh12.bin"
file_name_pmesh = "./data/kernel-c1-pmesh25.bin"
file_name_qmesh = "./data/kernel-c1-qmesh25.bin"


def read_binary_file(file_name, shape):
    with open(file_name, "rb") as f:
        data = np.fromfile(f, dtype=np.float64)
    return data.reshape(shape)


mtx_c1 = read_binary_file(file_name_c1, shape)
mtx_c3 = read_binary_file(file_name_c3, shape)
mtx_c4 = read_binary_file(file_name_c4, shape)
mtx_cD = read_binary_file(file_name_cD, shape)
mtx_cE = read_binary_file(file_name_cE, shape)
mtx_v3n = (
    c1 * mtx_c1
    + c3 * mtx_c3
    + c4 * mtx_c4
    + cD / lambda_chi * mtx_cD
    + cE / lambda_chi * mtx_cE
) / twopi6

mom_mesh_p = read_binary_file(file_name_pmesh, number_pmesh)
mom_mesh_q = read_binary_file(file_name_qmesh, number_qmesh)


def freg(p, q):
    dom = ((p**2 + 0.75 * q**2) / (lambda_3nf**2)) ** regulator_power
    return np.exp(-dom)


# interaction matrix with regulation.
mtx_v3n_reg = np.zeros_like(mtx_v3n)

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
                mtx_v3n_reg[idxp, idxq, idxpp, idxqp] = (
                    mtx_v3n[idxp, idxq, idxpp, idxqp] * regulator_pq * regulator_ppqp
                )

xi_values = np.zeros(number_pmesh)
for idxp in range(number_pmesh):
    xi_values[idxp] = np.sqrt(7.0 / 4.0 * mom_mesh_p[idxp] ** 2) / hbarc

xip, xi = np.meshgrid(xi_values, xi_values)

mtx_plot = np.zeros((number_pmesh, number_pmesh))
for idxpp in range(number_pmesh):
    for idxp in range(number_pmesh):
        mtx_plot[idxpp, idxp] = mtx_v3n_reg[idxpp, idxpp, idxp, idxp] * hbarc5

config = {
    "text.usetex": True,
    "mathtext.fontset": "stix",
    "font.family": "Times New Roman",
}
rcParams.update(config)
fig = plt.figure(dpi=200)

z_min, z_max = -np.abs(mtx_plot).max(), np.abs(mtx_plot).max()
plt.imshow(
    mtx_plot,
    extent=(xi.min(), xi.max(), xip.min(), xip.max()),
    cmap="RdBu_r",
    interpolation="bicubic",
    origin="lower",
    vmin=z_min,
    vmax=z_max,
)
plt.colorbar()
plt.xlabel(r"$\xi\;(\mathrm{fm^{-1}})$")
plt.ylabel(r"$\xi'\;(\mathrm{fm^{-1}})$")

plt.title("EM1.8/2.0")
plt.savefig("plot-3bme-v3n.png", bbox_inches="tight", transparent=False, dpi=200)
plt.show()
