import numpy as np

# generates 3N partial-waves in JJ-coupling scheme.


jmax = 5


def gen_pw_channels(twoJ, P, twoT):
    pw_channels = []
    for j in range(0, jmax + 1, 1):
        for s in [0, 1]:
            lmin = abs(j - s)
            lmax = abs(j + s)
            for l in range(lmin, lmax + 1, 1):
                twoj3min = int(abs(twoJ - 2 * j))
                twoj3max = int(abs(twoJ + 2 * j))
                for twoj3 in range(twoj3min, twoj3max + 1, 2):
                    lammin = int(abs(twoj3 - 1) / 2)
                    lammax = int(abs(twoj3 + 1) / 2)
                    for lam in range(lammin, lammax + 1, 1):
                        if P * (-1) ** (l + lam) < 0:
                            continue
                        t = 0
                        if (l + s) % 2 == 0:
                            t = 1
                        if twoT == 3 and t == 0:
                            continue
                        pw_channels.append([l, s, j, lam, twoj3, t])

    return pw_channels


def gen_J_channels(twoJmax):
    J_channels = []
    for twoT in [1, 3]:
        for twoJ in range(1, twoJmax + 1, 2):
            for P in [1, -1]:
                J_channels.append([twoJ, P, twoT])
    return J_channels


def main():
    alpha_channels = []
    (twoJ, P, twoT) = (9, -1, 3)
    pw_channels = gen_pw_channels(twoJ, P, twoT)
    for pw_channel in pw_channels:
        pw_channel.append(twoJ)
        pw_channel.append(P)
        pw_channel.append(twoT)
    alpha_channels.append(pw_channels)

    print("{l, s, j, lambda, twoj3, t, twoJ, P, twoT}")
    file_name = f"./tool/alpha_channels_twoJ{twoJ}_P{P}_twoT{twoT}.dat"
    with open(file_name, "w") as fp:
        for pw_channel in alpha_channels:
            for chan_bra in pw_channel:
                for chan_ket in pw_channel:
                    lp, sp, jp, lamp, twoj3p, tp = (
                        chan_bra[0],
                        chan_bra[1],
                        chan_bra[2],
                        chan_bra[3],
                        chan_bra[4],
                        chan_bra[5],
                    )
                    l, s, j, lam, twoj3, t = (
                        chan_ket[0],
                        chan_ket[1],
                        chan_ket[2],
                        chan_ket[3],
                        chan_ket[4],
                        chan_ket[5],
                    )
                    twoJ, P, twoT = chan_bra[6], chan_bra[7], chan_bra[8]
                    this_channel = (
                        "{"
                        + f"{lp},{sp},{jp},{lamp},{twoj3p},{tp},{l},{s},{j},{lam},{twoj3},{t},{twoJ},{twoT}"
                        + "},\n"
                    )
                    fp.write(this_channel)


if __name__ == "__main__":
    main()
