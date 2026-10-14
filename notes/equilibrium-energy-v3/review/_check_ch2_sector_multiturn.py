"""Sector-dipole FODO: Lambda from the turn-averaged path vs -2*pi*J*xi/C.

A single particle stays on its invariant torus, so the mean path per turn over
many turns is the torus average, free of launch-coordinate artefacts. The
per-turn path oscillates at order sqrt(J); a Hann window suppresses it.
"""
import math

from _check_ch2_chrom import track_sector, solve4, tune_unfolded, twiss_unfolded

RHO, THETA, K, N_CELL = 8.0, 0.25, 0.30, 2
C = N_CELL * 2.0 * RHO * THETA


def track_ring(x, xp, y, yp, delta):
    path = 0.0
    for _ in range(N_CELL):
        for k in (K, -K):
            xp -= k / (1.0 + delta) * x
            yp += k / (1.0 + delta) * y
            x, xp, y, yp, p = track_sector(x, xp, y, yp, delta, RHO, THETA)
            path += p
    return x, xp, y, yp, path


def fixed(delta):
    v = [0.0] * 4
    for _ in range(12):
        q = track_ring(*v, delta)
        r = [q[i] - v[i] for i in range(4)]
        if sum(map(abs, r)) < 1e-14:
            break
        e = 1e-8
        cols = []
        for j in range(4):
            w = v[:]
            w[j] += e
            qj = track_ring(*w, delta)
            cols.append([(qj[i] - w[i] - r[i]) / e for i in range(4)])
        step = solve4([[cols[j][i] for j in range(4)] for i in range(4)], [-t for t in r])
        v = [v[i] + step[i] for i in range(4)]
    return v


def matrices(delta, v):
    e = 1e-7
    cols = []
    for j in range(4):
        wp, wm = v[:], v[:]
        wp[j] += e
        wm[j] -= e
        qp, qm = track_ring(*wp, delta), track_ring(*wm, delta)
        cols.append([(qp[i] - qm[i]) / (2 * e) for i in range(4)])
    Mx = ((cols[0][0], cols[1][0]), (cols[0][1], cols[1][1]))
    My = ((cols[2][2], cols[3][2]), (cols[2][3], cols[3][3]))
    return Mx, My


def main():
    v0 = fixed(0.0)
    Mx, My = matrices(0.0, v0)
    bx, ax, _ = twiss_unfolded(Mx)
    by, ay, _ = twiss_unfolded(My)
    eps = 1e-4
    vp, vm = fixed(eps), fixed(-eps)
    xi = []
    for i in range(2):
        tp = tune_unfolded(matrices(eps, vp)[i])
        tm = tune_unfolded(matrices(-eps, vm)[i])
        xi.append((tp - tm) / (2 * eps))
    p0 = track_ring(*v0, 0.0)[4]
    print("xi", xi, "C", C)
    n_turns = 4000
    for label, Jx, Jy in (("Jx", 1e-6, 0.0), ("Jy", 0.0, 1.5e-6)):
        lams = []
        for scale in (1.0, 0.25):
            jx, jy = Jx * scale, Jy * scale
            psi, psy = 0.3, 1.1
            x = v0[0] + math.sqrt(2 * jx * bx) * math.cos(psi)
            xp = v0[1] - math.sqrt(2 * jx / bx) * (math.sin(psi) + ax * math.cos(psi))
            y = math.sqrt(2 * jy * by) * math.cos(psy)
            yp = -math.sqrt(2 * jy / by) * (math.sin(psy) + ay * math.cos(psy))
            acc = wsum = 0.0
            for n in range(n_turns):
                x, xp, y, yp, p = track_ring(x, xp, y, yp, 0.0)
                w = math.sin(math.pi * (n + 0.5) / n_turns) ** 2
                acc += w * (p - p0)
                wsum += w
            lams.append(acc / wsum / C)
        lin = (16.0 * lams[1] - lams[0]) / 3.0
        pred = -(2 * math.pi / C) * (Jx * xi[0] + Jy * xi[1])
        print(label, "direct", lin, "xi-pred", pred, "rel", (lin - pred) / pred)


if __name__ == "__main__":
    main()
