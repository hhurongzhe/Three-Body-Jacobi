import numpy as np
import pandas as pd
import WignerSymbol as ws
import time
import functools
from functools import lru_cache


hbarc = 197.32698
hbarc5 = hbarc**5
twopi6 = (2.0 * np.pi) ** 6


def timing(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f"function '{func.__name__}' timing: {elapsed_time:.6f} seconds")
        return result

    return wrapper


@lru_cache(maxsize=None)
def read_binary_file(file_name, shape):
    with open(file_name, "rb") as f:
        data = np.fromfile(f, dtype=np.float64)
    return data.reshape(shape)


def write_binary_file(data_array: np.ndarray, file_name: str):
    try:
        data_array.tofile(file_name)
        print(f"data written in file: {file_name}")
    except IOError as e:
        print(f"error writing file: {e}")


def read_beta_pw_info(file_name: str):
    try:
        df = pd.read_csv(file_name, sep=r"\s+", index_col="index", engine="python")
        quantum_dict = {tuple(row): idx for idx, row in df.iterrows()}
        return quantum_dict
    except FileNotFoundError:
        print(f"error: cannot find file: '{file_name}'")
        return {}


# generate 3N partial-wave channels in JJ-scheme
# in the same order of Hebeler's review, https://doi.org/10.1016/j.physrep.2020.08.009, page 102
def gen_alpha_pw_channels(twoJ, P, twoT, jmax):
    pw_channels = []
    for s in [0, 1]:
        for j in range(0, jmax + 1, 1):
            lmin, lmax = abs(j - s), abs(j + s)
            for l in range(lmin, lmax + 1, 1):
                twoj3min, twoj3max = abs(twoJ - 2 * j), abs(twoJ + 2 * j)
                for twoj3 in range(twoj3min, twoj3max + 1, 2):
                    lammin, lammax = int(abs(twoj3 - 1) / 2), int(abs(twoj3 + 1) / 2)
                    for lam in range(lammin, lammax + 1, 1):
                        if P * (-1) ** (l + lam) > 0:
                            if (l + s) % 2 == 0:
                                t = 1
                            else:
                                t = 0
                            if not (t == 0 and twoT == 3):
                                pw_channels.append([l, s, j, t, lam, twoj3])
    pw_channels_sorted = sorted(pw_channels, key=lambda x: x[0])
    return pw_channels_sorted


def save_alpha_channel_info(pw_channels, filename):
    with open(filename, "w") as f:
        f.write("index l s j t lam 2j3\n")
        for idx_chan, chan in enumerate(pw_channels):
            [l, s, j, t, lam, twoj3] = chan
            f.write(f"{idx_chan+1} {l} {s} {j} {t} {lam} {twoj3}\n")
    print(f"Channel information saved to {filename}")


@lru_cache(maxsize=None)
def hat(j):
    two_j = int(2 * j)
    return np.sqrt(two_j + 1)


@lru_cache(maxsize=None)
def f9j(j1, j2, j3, j4, j5, j6, j7, j8, j9):
    dj1 = int(2.0 * j1)
    dj2 = int(2.0 * j2)
    dj3 = int(2.0 * j3)
    dj4 = int(2.0 * j4)
    dj5 = int(2.0 * j5)
    dj6 = int(2.0 * j6)
    dj7 = int(2.0 * j7)
    dj8 = int(2.0 * j8)
    dj9 = int(2.0 * j9)
    result = ws.f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)
    return result


# mesh number
number_pmesh = 16
number_qmesh = 16
shape = (number_pmesh, number_qmesh, number_pmesh, number_qmesh)

# JJ-scheme parameters
[twoJ, P, twoT] = [1, 1, 1]
jmax = 2

# 3NF part and mesh info
part, nmesh = "c1", 8


@timing
def main():
    beta_channel_file = f"data/channel_beta_info_twoJ{twoJ}_P{P}_twoT{twoT}.txt"
    beta_pw_info = read_beta_pw_info(beta_channel_file)
    ws.init(24, "Jmax", 9)
    alpha_channels = gen_alpha_pw_channels(twoJ, P, twoT, jmax)
    alpha_channel_file = f"data/channel_alpha_info_twoJ{twoJ}_P{P}_twoT{twoT}.txt"
    save_alpha_channel_info(alpha_channels, alpha_channel_file)
    for idx_bra, chan_bra in enumerate(alpha_channels):
        for idx_ket, chan_ket in enumerate(alpha_channels):
            (l_bra, s_bra, j_bra, t_bra, lam_bra, twoj3_bra) = chan_bra
            (l_ket, s_ket, j_ket, t_ket, lam_ket, twoj3_ket) = chan_ket
            chan_mtx = np.zeros(shape, dtype=np.float64)
            for twoS_bra in [1, 3]:
                for twoS_ket in [1, 3]:
                    Lmin_bra, Lmax_bra = int(abs(twoJ - twoS_bra) / 2), int(abs(twoJ + twoS_bra) / 2)
                    Lmin_ket, Lmax_ket = int(abs(twoJ - twoS_ket) / 2), int(abs(twoJ + twoS_ket) / 2)
                    for L_bra in range(Lmin_bra, Lmax_bra + 1, 1):
                        for L_ket in range(Lmin_ket, Lmax_ket + 1, 1):
                            beta_key_bra = l_bra, lam_bra, L_bra, s_bra, twoS_bra, t_bra, twoT, twoJ
                            beta_key_ket = l_ket, lam_ket, L_ket, s_ket, twoS_ket, t_ket, twoT, twoJ
                            if beta_key_bra in beta_pw_info and beta_key_ket in beta_pw_info:
                                chan_beta_index_bra = beta_pw_info[beta_key_bra]
                                chan_beta_index_ket = beta_pw_info[beta_key_ket]
                                beta_mtx_name = f"data/kernel-{part}-bra{chan_beta_index_bra}-ket{chan_beta_index_ket}-nmesh{nmesh}.bin"
                                beta_mtx = read_binary_file(beta_mtx_name, shape)
                                factor = np.sqrt(hat(L_bra) * hat(L_ket) * hat(twoS_bra / 2) * hat(twoS_ket / 2) * hat(j_bra) * hat(j_ket) * hat(twoj3_bra / 2) * hat(twoj3_ket / 2))
                                factor *= f9j(l_bra, s_bra, j_bra, lam_bra, 1 / 2, twoj3_bra / 2, L_bra, twoS_bra / 2, twoJ / 2) * f9j(l_ket, s_ket, j_ket, lam_ket, 1 / 2, twoj3_ket / 2, L_ket, twoS_ket / 2, twoJ / 2)
                                chan_mtx += factor * beta_mtx
            alpha_mtx_file = f"data/kernel-alpha-{part}-bra{idx_bra + 1}-ket{idx_ket + 1}-nmesh{nmesh}.bin"
            write_binary_file(chan_mtx, alpha_mtx_file)


if __name__ == "__main__":
    main()
