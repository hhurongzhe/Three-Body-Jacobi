import numpy as np


def get_mesh_ori(precision, deg):
    np.set_printoptions(precision)
    temp = np.polynomial.legendre.leggauss(deg)
    xlist = []
    wlist = []
    for i in range(0, len(temp[0])):
        xlist.append(temp[0][i])
        wlist.append(temp[1][i])
    return xlist, wlist


def map(x, w, a, b):
    xx = []
    ww = []
    for _, xi in enumerate(x):
        wi = w[_]
        xxi = (b - a) / 2.0 * xi + (b + a) / 2.0
        wwi = (b - a) / 2.0 * wi
        xx.append(xxi)
        ww.append(wwi)
    return xx, ww


xlist1, wlist1 = get_mesh_ori(precision=16, deg=8)
xlist2, wlist2 = get_mesh_ori(precision=16, deg=16)

mesh_theta, weight_theta = map(xlist1, wlist1, 0, np.pi)
mesh_phi, weight_phi = map(xlist1, wlist1, 0, 2 * np.pi)
mesh_mom, weight_mom = map(xlist2, wlist2, 0, 600)

print("original mesh in [-1,1] :\n")
print(xlist1, sep=",", end="\n\n", flush=False)
print(wlist1, sep=",", end="\n\n", flush=False)

print("theta mesh :\n")
print(mesh_theta, sep=",", end="\n\n", flush=False)
print(weight_theta, sep=",", end="\n\n", flush=False)

print("phi mesh :\n")
print(mesh_phi, sep=",", end="\n\n", flush=False)
print(weight_phi, sep=",", end="\n\n", flush=False)

print("mom mesh :\n")
print(mesh_mom, sep=",", end="\n\n", flush=False)
print(weight_mom, sep=",", end="\n\n", flush=False)
