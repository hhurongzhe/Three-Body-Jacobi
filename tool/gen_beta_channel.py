import numpy as np


def gen_pw_channels(twoJ, P, twoT, lmax, lammax):
    pw_channels = []
    for twoS in [1, 3]:
        if twoS == 1:
            s_set = [0, 1]
        else:
            s_set = [1]
        for s in s_set:
            Lmin = int(abs(twoJ - twoS) / 2)
            Lmax = int(abs(twoJ + twoS) / 2)
            for L in range(Lmin, Lmax + 1, 1):
                for l in range(0, lmax + 1, 1):
                    for lam in range(0, lammax + 1, 1):
                        if (
                            abs(l - lam) <= L
                            and (l + lam) >= L
                            and L <= abs(l + lam)
                            and P * (-1) ** (l + lam) > 0
                        ):
                            if (l + s) % 2 == 0:
                                t = 1
                            else:
                                t = 0
                            pw_channels.append([l, lam, L, s, twoS, t])
    pw_channels_sorted = sorted(pw_channels, key=lambda x: x[0])
    return pw_channels_sorted


def gen_J_channels(twoJmax):
    J_channels = []
    for twoT in [1, 3]:
        if twoJmax == 1:
            twoJ_set = [1]
        else:
            twoJ_set = range(1, twoJmax + 1, 2)
        for twoJ in twoJ_set:
            for P in [1, -1]:
                J_channels.append([twoJ, P, twoT])
    return J_channels


def main():
    # twoJmax = 9
    # lmax, lammax = 6, 6

    # beta_channels = []
    # J_channels = gen_J_channels(twoJmax)
    # for J_channel in J_channels:
    #     [twoJ, P, twoT] = J_channel
    #     pw_channels = gen_pw_channels(twoJ, P, twoT, lmax, lammax)
    #     for pw_channel in pw_channels:
    #         beta_channels.append(pw_channel + [twoJ, P, twoT])

    # for chan in beta_channels:
    #     print(chan)
    # print(len(beta_channels))

    [twoJ, P, twoT] = [1, 1, 1]
    [lmax, lammax] = [int((twoJ + 3) / 2), int((twoJ + 3) / 2)]
    pw_channels = gen_pw_channels(twoJ, P, twoT, lmax, lammax)
    for chan in pw_channels:
        chan_str = str(chan).replace("[", "{").replace("]", "}") + ","
        print(chan_str)
    print(len(pw_channels))


if __name__ == "__main__":
    main()
