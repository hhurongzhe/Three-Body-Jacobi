import numpy as np
from scipy.integrate import nquad
from functools import lru_cache
import time

gA = 1.29
fpi = 92.4
mpi = 138.0
c1 = -0.81 * 1e-3
c3 = -3.4 * 1e-3
c4 = 3.4 * 1e-3
hbarc = 197.32698


@lru_cache(maxsize=None)
def get_q2q2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp):
    return (
        np.power(
            2 * pmag
            - 2 * ppmag * np.cos(thetapp)
            - qmag * np.cos(thetaq)
            + qpmag * np.cos(thetaqp),
            2,
        )
        + np.power(
            2 * ppmag * np.cos(phipp) * np.sin(thetapp)
            + qmag * np.sin(thetaq)
            - qpmag * np.cos(phiqp) * np.sin(thetaqp),
            2,
        )
        + np.power(
            -2 * ppmag * np.sin(phipp) * np.sin(thetapp)
            + qpmag * np.sin(phiqp) * np.sin(thetaqp),
            2,
        )
    ) / 4.0


@lru_cache(maxsize=None)
def get_q3q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp):
    return (
        np.power(
            2 * pmag
            - 2 * ppmag * np.cos(thetapp)
            + qmag * np.cos(thetaq)
            - qpmag * np.cos(thetaqp),
            2,
        )
        + np.power(
            2 * ppmag * np.cos(phipp) * np.sin(thetapp)
            - qmag * np.sin(thetaq)
            + qpmag * np.cos(phiqp) * np.sin(thetaqp),
            2,
        )
        + np.power(
            2 * ppmag * np.sin(phipp) * np.sin(thetapp)
            + qpmag * np.sin(phiqp) * np.sin(thetaqp),
            2,
        )
    ) / 4.0


@lru_cache(maxsize=None)
def get_q2q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp):
    return (
        -4 * np.power(pmag, 2)
        - 4 * np.power(ppmag, 2)
        + np.power(qmag, 2)
        + np.power(qpmag, 2)
        + 8 * pmag * ppmag * np.cos(thetapp)
        - 2 * qmag * qpmag * np.cos(thetaq) * np.cos(thetaqp)
        - 2 * qmag * qpmag * np.cos(phiqp) * np.sin(thetaq) * np.sin(thetaqp)
    ) / 4.0


@lru_cache(maxsize=None)
def get_F1(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp):
    global gA, fpi, mpi, c1, c3, c4
    q2q2 = get_q2q2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    q3q3 = get_q3q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    q2q3 = get_q2q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    part1 = (gA / (2 * fpi)) ** 2
    part2 = 1 / (q2q2 + mpi * mpi) / (q3q3 + mpi * mpi)
    part3 = -4 * c1 * mpi * mpi / (fpi * fpi) + 2 * c3 / (fpi * fpi) * q2q3
    F1 = part1 * part2 * part3
    return F1


@lru_cache(maxsize=None)
def get_F2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp):
    global gA, fpi, mpi, c1, c3, c4
    q2q2 = get_q2q2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    q3q3 = get_q3q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    part1 = (gA / (2 * fpi)) ** 2
    part2 = 1 / (q2q2 + mpi * mpi) / (q3q3 + mpi * mpi)
    part3 = c4 / (fpi * fpi)
    F2 = part1 * part2 * part3
    return F2


@lru_cache(maxsize=None)
def GG(
    thetaq,
    thetapp,
    phipp,
    thetaqp,
    phiqp,
    lp,
    lamp,
    Lp,
    sp,
    twoSp,
    l,
    lam,
    L,
    s,
    twoS,
    twoJ,
    pmag,
    qmag,
    ppmag,
    qpmag,
):
    F1 = get_F1(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    # F2 = get_F2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    spin_tuple = (lp, lamp, Lp, sp, twoSp, l, lam, L, s, twoS, twoJ)
    if spin_tuple == (0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1):
        mtx = (
            F1
            * (
                4 * np.power(pmag, 2)
                + 4 * np.power(ppmag, 2)
                - np.power(qmag, 2)
                - np.power(qpmag, 2)
                - 8 * pmag * ppmag * np.cos(thetapp)
                + 2 * qmag * qpmag * np.cos(thetaq) * np.cos(thetaqp)
                + 2 * qmag * qpmag * np.cos(phiqp) * np.sin(thetaq) * np.sin(thetaqp)
            )
        ) / (64.0 * np.power(np.pi, 2))
        return mtx * np.sin(thetaq) * np.sin(thetapp) * np.sin(thetaqp) * 8 * np.pi**2
    else:
        raise ValueError("Invalid spin configuration in function GG!")


@lru_cache(maxsize=None)
def GG_bench(
    thetaq,
    thetapp,
    phipp,
    thetaqp,
    phiqp,
    lp,
    lamp,
    Lp,
    sp,
    twoSp,
    l,
    lam,
    L,
    s,
    twoS,
    twoJ,
    pmag,
    qmag,
    ppmag,
    qpmag,
):
    F1 = get_F1(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    # F2 = get_F2(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    q2q3 = get_q2q3(pmag, qmag, ppmag, qpmag, thetaq, thetapp, phipp, thetaqp, phiqp)
    spin_tuple = (lp, lamp, Lp, sp, twoSp, l, lam, L, s, twoS, twoJ)
    if spin_tuple == (0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1):
        mtx = -q2q3 * F1 / (16.0 * np.power(np.pi, 2))
        return mtx * np.sin(thetaq) * np.sin(thetapp) * np.sin(thetaqp) * 8 * np.pi**2
    else:
        raise ValueError("Invalid spin configuration in function GG!")


def G(lp, lamp, Lp, sp, twoSp, l, lam, L, s, twoS, twoJ, pmag, qmag, ppmag, qpmag):
    global hbarc
    thetaq_limits = [0, np.pi]
    thetapp_limits = [0, np.pi]
    phipp_limits = [0, 2 * np.pi]
    thetaqp_limits = [0, np.pi]
    phiqp_limits = [0, 2 * np.pi]
    eps = 1e-12
    opts1 = {"epsrel": eps}
    opts2 = {"epsrel": eps}
    opts3 = {"epsrel": eps}
    opts4 = {"epsrel": eps}
    opts5 = {"epsrel": eps}
    result, error = nquad(
        GG,
        [thetaq_limits, thetapp_limits, phipp_limits, thetaqp_limits, phiqp_limits],
        args=(
            lp,
            lamp,
            Lp,
            sp,
            twoSp,
            l,
            lam,
            L,
            s,
            twoS,
            twoJ,
            pmag,
            qmag,
            ppmag,
            qpmag,
        ),
        opts=[opts1, opts2, opts3, opts4, opts5],
    )
    return result * hbarc**5


lp, lamp, Lp, sp, twoSp = 0, 0, 0, 0, 1
l, lam, L, s, twoS = 0, 0, 0, 0, 1
twoJ = 1
pmag, qmag, ppmag, qpmag = 1 * hbarc, 2 * hbarc, 3 * hbarc, 4 * hbarc

t1 = time.time()
re = G(lp, lamp, Lp, sp, twoSp, l, lam, L, s, twoS, twoJ, pmag, qmag, ppmag, qpmag)
t2 = time.time()

print(f"result: {re}")
print(f"time: {t2-t1}")
