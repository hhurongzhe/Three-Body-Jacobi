import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

hbarc = 197.32698


# mesh number
number_pmesh = 25
number_qmesh = 25


mtx_c1 = np.load("vc1.npy")
mtx_c3 = np.load("vc3.npy")
mtx_c4 = np.load("vc4.npy")
mtx_cD = np.load("vcD.npy")
mtx_cE = np.load("vcE.npy")


def get_matrix(mtx):
    mtx_plot = np.zeros((number_pmesh, number_pmesh))
    for idxpp in range(number_pmesh):
        for idxp in range(number_pmesh):
            mtx_plot[idxpp, idxp] = mtx[idxpp, idxpp, idxp, idxp]
    return mtx_plot


mtx_c1_plot = get_matrix(mtx_c1)
mtx_c3_plot = get_matrix(mtx_c3)
mtx_c4_plot = get_matrix(mtx_c4)
mtx_cD_plot = get_matrix(mtx_cD)
mtx_cE_plot = get_matrix(mtx_cE)
file_name_pmesh = "./data/kernel-c1-pmesh25.bin"
file_name_qmesh = "./data/kernel-c1-qmesh25.bin"


def read_binary_file(file_name, shape):
    with open(file_name, "rb") as f:
        data = np.fromfile(f, dtype=np.float64)
    return data.reshape(shape)


mom_mesh_p = read_binary_file(file_name_pmesh, number_pmesh)
mom_mesh_q = read_binary_file(file_name_qmesh, number_qmesh)

xi_values = np.zeros(number_pmesh)
for idxp in range(number_pmesh):
    xi_values[idxp] = np.sqrt(7.0 / 4.0 * mom_mesh_p[idxp] ** 2) / hbarc

xip, xi = np.meshgrid(xi_values, xi_values)

config = {
    "text.usetex": True,
    "mathtext.fontset": "stix",
    "font.family": "Times New Roman",
}
rcParams.update(config)

# 创建图形和子图布局
fig, axes = plt.subplots(1, 5, figsize=(15, 3), dpi=200, sharey=True)

# 矩阵的颜色范围设置
z_min = min(
    [
        -np.abs(mtx_c1_plot).max(),
        -np.abs(mtx_c3_plot).max(),
        -np.abs(mtx_c4_plot).max(),
        -np.abs(mtx_cD_plot).max(),
        -np.abs(mtx_cE_plot).max(),
    ]
)
z_max = max(
    [
        np.abs(mtx_c1_plot).max(),
        np.abs(mtx_c3_plot).max(),
        np.abs(mtx_c4_plot).max(),
        np.abs(mtx_cD_plot).max(),
        np.abs(mtx_cE_plot).max(),
    ]
)

idx_s = 0
c = axes[idx_s].imshow(
    mtx_c1_plot,
    cmap="RdBu_r",
    interpolation="bicubic",
    extent=(xi.min(), xi.max(), xip.min(), xip.max()),
    origin="lower",
    vmin=z_min,
    vmax=z_max,
)
axes[idx_s].set_xlabel(r"$\xi\;(\mathrm{fm^{-1}})$")
axes[idx_s].set_title(r"$c_1$")

idx_s = 1
c = axes[idx_s].imshow(
    mtx_c3_plot,
    cmap="RdBu_r",
    interpolation="bicubic",
    extent=(xi.min(), xi.max(), xip.min(), xip.max()),
    origin="lower",
    vmin=z_min,
    vmax=z_max,
)
axes[idx_s].set_xlabel(r"$\xi\;(\mathrm{fm^{-1}})$")
axes[idx_s].set_title(r"$c_3$")

idx_s = 2
c = axes[idx_s].imshow(
    mtx_c4_plot,
    cmap="RdBu_r",
    interpolation="bicubic",
    extent=(xi.min(), xi.max(), xip.min(), xip.max()),
    origin="lower",
    vmin=z_min,
    vmax=z_max,
)
axes[idx_s].set_xlabel(r"$\xi\;(\mathrm{fm^{-1}})$")
axes[idx_s].set_title(r"$c_4$")

idx_s = 3
c = axes[idx_s].imshow(
    mtx_cD_plot,
    cmap="RdBu_r",
    interpolation="bicubic",
    extent=(xi.min(), xi.max(), xip.min(), xip.max()),
    origin="lower",
    vmin=z_min,
    vmax=z_max,
)
axes[idx_s].set_xlabel(r"$\xi\;(\mathrm{fm^{-1}})$")
axes[idx_s].set_title(r"$c_D$")

idx_s = 4
c = axes[idx_s].imshow(
    mtx_cE_plot,
    cmap="RdBu_r",
    interpolation="bicubic",
    extent=(xi.min(), xi.max(), xip.min(), xip.max()),
    origin="lower",
    vmin=z_min,
    vmax=z_max,
)
axes[idx_s].set_xlabel(r"$\xi\;(\mathrm{fm^{-1}})$")
axes[idx_s].set_title(r"$c_E$")

axes[0].set_ylabel(r"$\xi'\;(\mathrm{fm^{-1}})$")
plt.subplots_adjust(wspace=0, right=1)
cbar = fig.colorbar(c, ax=axes, orientation="vertical", pad=0.01)
cbar.set_label("$V\,\mathrm{(fm^{5})}$")

plt.savefig("plot-3bme-5part.png", bbox_inches="tight", transparent=False, dpi=200)
plt.show()
