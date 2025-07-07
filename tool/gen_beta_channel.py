def gen_beta_pw_channels(twoJ, P, twoT, lmax, lammax):
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
                        if abs(l - lam) <= L and (l + lam) >= L and L <= abs(l + lam) and P * (-1) ** (l + lam) > 0:
                            if (l + s) % 2 == 0:
                                t = 1
                            else:
                                t = 0
                            if not (t == 0 and twoT == 3):
                                pw_channels.append([l, lam, L, s, twoS, t, twoT, twoJ, P])
    pw_channels_sorted = sorted(pw_channels, key=lambda x: x[0])
    return pw_channels_sorted


def print_for_mma(pw_channels):
    print("{l,lam,L,s,2S,t,2T,2J}")
    for idx_chan, chan in enumerate(pw_channels):
        [l, lam, L, s, twoS, t, twoT, twoJ, P] = chan
        chan_str = "{" + str(l) + "," + str(lam) + "," + str(L) + "," + str(s) + "," + str(twoS) + "," + str(t) + "," + str(twoT) + "," + str(twoJ) + "},"
        print(idx_chan + 1, chan_str)
    print(len(pw_channels))


def save_beta_channel_info(pw_channels, filename):
    with open(filename, "w") as f:
        f.write("index l lam L s 2S t 2T 2J\n")
        for idx_chan, chan in enumerate(pw_channels):
            [l, lam, L, s, twoS, t, twoT, twoJ, P] = chan
            f.write(f"{idx_chan+1} {l} {lam} {L} {s} {twoS} {t} {twoT} {twoJ}\n")
    print(f"Channel information saved to {filename}")


def main():
    [twoJ, P, twoT] = [1, 1, 1]
    [lmax, lammax] = [int((twoJ + 3) / 2), int((twoJ + 3) / 2)]
    pw_channels = gen_beta_pw_channels(twoJ, P, twoT, lmax, lammax)
    print_for_mma(pw_channels)
    channel_file = f"data/channel_beta_info_twoJ{twoJ}_P{P}_twoT{twoT}.txt"
    save_beta_channel_info(pw_channels, channel_file)


if __name__ == "__main__":
    main()
