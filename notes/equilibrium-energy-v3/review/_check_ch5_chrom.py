"""Electrostatic cell: turn time, chromaticity, and the working integral for Lambda^E.

Periodic cell: hard-edge cylindrical deflectors, drifts, thin magnetic
quadrupoles, optional thin sextupoles. gamma mc^2 + q Phi is constant
between cavities. The edge is a radial plane. Transverse momenta p_r and
p_y are continuous across it; the azimuthal momentum jumps so that the
kinetic momentum matches Phi = 0 outside and Phi(r, y) inside.

delta is the relative momentum outside the deflectors. One particle is
tracked for many turns; the turn time is averaged with a Hann window, and
the piece linear in the action is taken from amplitudes J and J/4.
"""
import math
import sys

from _check_ch2_chrom import solve4, tune_unfolded, twiss_unfolded

GAMMA_S = 1.25
BETA_S = math.sqrt(1.0 - 1.0 / GAMMA_S**2)
P_S = BETA_S * GAMMA_S
V_S = BETA_S
PVS = BETA_S**2 * GAMMA_S

RHO = 8.0
THETA = 0.40
DRIFT = 0.8
K_QUAD = 0.15
N_CELL = 2
N_STEP = 200
N_TURNS = 1200

# Reference length of one period. Four drifts and two sectors per cell.
C_REF = N_CELL * 2.0 * (2.0 * DRIFT + RHO * THETA)


def field_coeffs(n):
    # psi = f0(r) + f2(r) y^2 + f4(r) y^4, Laplace through y^4,
    # median field E_r = E_0 (rho/r)^n.
    c2 = 0.5 * (n - 1.0)
    c4 = -(n - 1.0) * (n + 1.0) ** 2 / 24.0
    return c2, c4


def psi_forces(r, y, n, c2, c4):
    rho_n1 = RHO ** (n - 1.0)
    inv_r = 1.0 / r
    rn = (RHO * inv_r) ** n
    f0 = (1.0 - (RHO * inv_r) ** (n - 1.0)) / (n - 1.0)
    f0p = rn / RHO
    rm = rho_n1 * inv_r ** (n + 1.0)
    f2 = c2 * rm
    f2p = c2 * (-(n + 1.0)) * rm * inv_r
    rm4 = rm * inv_r * inv_r
    f4 = c4 * rm4
    f4p = c4 * (-(n + 3.0)) * rm4 * inv_r
    y2 = y * y
    psi = f0 + (f2 + f4 * y2) * y2
    dpsi_dr = f0p + (f2p + f4p * y2) * y2
    dpsi_dy = (2.0 * f2 + 4.0 * f4 * y2) * y
    return psi, -PVS * dpsi_dr, -PVS * dpsi_dy


def deriv(r, pr, y, py, L, g_out, n, c2, c4):
    psi, qer, qey = psi_forces(r, y, n, c2, c4)
    g = g_out - PVS * psi
    pt = L / r
    inv = g * r / pt
    return (
        pr * r / pt,
        pt + qer * inv,
        py * r / pt,
        qey * inv,
        inv,
    )


def rk4(r, pr, y, py, L, g_out, n, c2, c4, dth):
    k1 = deriv(r, pr, y, py, L, g_out, n, c2, c4)
    k2 = deriv(
        r + 0.5 * dth * k1[0], pr + 0.5 * dth * k1[1],
        y + 0.5 * dth * k1[2], py + 0.5 * dth * k1[3],
        L, g_out, n, c2, c4,
    )
    k3 = deriv(
        r + 0.5 * dth * k2[0], pr + 0.5 * dth * k2[1],
        y + 0.5 * dth * k2[2], py + 0.5 * dth * k2[3],
        L, g_out, n, c2, c4,
    )
    k4 = deriv(
        r + dth * k3[0], pr + dth * k3[1],
        y + dth * k3[2], py + dth * k3[3],
        L, g_out, n, c2, c4,
    )
    out = []
    for i in range(5):
        out.append((k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0 * dth)
    return r + out[0], pr + out[1], y + out[2], py + out[3], out[4]


def track_sector(x, xp, y, yp, p_out, g_out, n, c2, c4, nstep):
    r = RHO + x
    norm = math.sqrt(1.0 + xp * xp + yp * yp)
    pr = p_out * xp / norm
    py = p_out * yp / norm
    psi, _, _ = psi_forces(r, y, n, c2, c4)
    g = g_out - PVS * psi
    pt2 = g * g - 1.0 - pr * pr - py * py
    if pt2 <= 0.0:
        raise RuntimeError("particle reflected at deflector entrance")
    pt = math.sqrt(pt2)
    L = pt * r
    dth = THETA / nstep
    t = 0.0
    xint = 0.0
    err = 0.0
    for _ in range(nstep):
        x0 = r - RHO
        r, pr, y, py, dt = rk4(r, pr, y, py, L, g_out, n, c2, c4, dth)
        t += dt
        xint += 0.5 * (x0 + (r - RHO)) * dth
    psi, _, _ = psi_forces(r, y, n, c2, c4)
    g = g_out - PVS * psi
    p2 = pr * pr + py * py + (L / r) ** 2
    scale = max(g * g - 1.0, 1e-30)
    err = abs(p2 - (g * g - 1.0)) / scale
    pt2 = p_out * p_out - pr * pr - py * py
    if pt2 <= 0.0:
        raise RuntimeError("particle reflected at deflector exit")
    pt = math.sqrt(pt2)
    return r - RHO, pr / pt, y, py / pt, t, xint, err


def track_period(x, xp, y, yp, delta, n, c2, c4, k_quad, m_sext, nstep):
    p_out = P_S * (1.0 + delta)
    g_out = math.sqrt(1.0 + p_out * p_out)
    v = p_out / g_out
    scale = 1.0 / (1.0 + delta)
    t = 0.0
    xint = 0.0
    err = 0.0
    for _c in range(N_CELL):
        for kq, ms in ((k_quad, m_sext), (-k_quad, -m_sext)):
            xp -= kq * scale * x
            yp += kq * scale * y
            if ms != 0.0:
                xp -= 0.5 * ms * scale * (x * x - y * y)
                yp += ms * scale * x * y
            path = DRIFT * math.sqrt(1.0 + xp * xp + yp * yp)
            t += path / v
            x += DRIFT * xp
            y += DRIFT * yp
            x, xp, y, yp, ts, xi, er = track_sector(
                x, xp, y, yp, p_out, g_out, n, c2, c4, nstep
            )
            t += ts
            xint += xi
            if er > err:
                err = er
            path = DRIFT * math.sqrt(1.0 + xp * xp + yp * yp)
            t += path / v
            x += DRIFT * xp
            y += DRIFT * yp
    return x, xp, y, yp, t, xint, err


def fixed_point(delta, n, c2, c4, k_quad, m_sext, nstep):
    v = [0.0, 0.0, 0.0, 0.0]
    for _ in range(14):
        q = track_period(*v, delta, n, c2, c4, k_quad, m_sext, nstep)
        r = [q[i] - v[i] for i in range(4)]
        if sum(abs(t) for t in r) < 1e-13:
            break
        e = 1e-8
        cols = []
        for j in range(4):
            w = v[:]
            w[j] += e
            qj = track_period(*w, delta, n, c2, c4, k_quad, m_sext, nstep)
            cols.append([(qj[i] - w[i] - r[i]) / e for i in range(4)])
        step = solve4(
            [[cols[j][i] for j in range(4)] for i in range(4)],
            [-t for t in r],
        )
        v = [v[i] + step[i] for i in range(4)]
    return v


def one_turn_matrix(delta, v, n, c2, c4, k_quad, m_sext, nstep):
    e = 1e-6
    cols = []
    for j in range(4):
        wp, wm = v[:], v[:]
        wp[j] += e
        wm[j] -= e
        qp = track_period(*wp, delta, n, c2, c4, k_quad, m_sext, nstep)
        qm = track_period(*wm, delta, n, c2, c4, k_quad, m_sext, nstep)
        cols.append([(qp[i] - qm[i]) / (2.0 * e) for i in range(4)])
    mx = ((cols[0][0], cols[1][0]), (cols[0][1], cols[1][1]))
    my = ((cols[2][2], cols[3][2]), (cols[2][3], cols[3][3]))
    return mx, my


def m_mul(A, B):
    return (
        (A[0][0] * B[0][0] + A[0][1] * B[1][0],
         A[0][0] * B[0][1] + A[0][1] * B[1][1]),
        (A[1][0] * B[0][0] + A[1][1] * B[1][0],
         A[1][0] * B[0][1] + A[1][1] * B[1][1]),
    )


def twiss_through(M, beta, alpha):
    g = (1.0 + alpha * alpha) / beta
    b2 = M[0][0] ** 2 * beta - 2.0 * M[0][0] * M[0][1] * alpha + M[0][1] ** 2 * g
    a2 = (
        -M[0][0] * M[1][0] * beta
        + (M[0][0] * M[1][1] + M[0][1] * M[1][0]) * alpha
        - M[0][1] * M[1][1] * g
    )
    return b2, a2


def optics_integrals(n, c2, c4, k_quad, beta_x, alpha_x, beta_y, alpha_y, nstep):
    """Integrate gamma and beta/rho^2 along the linear optics of the tracker.

    Four trajectories carry the transfer matrix from the period entrance.
    Edges, drifts and thin quadrupoles are the same maps as in track_period.
    Sextupoles are omitted: on the reference orbit their feed-down is zero.
    """
    e = 1e-6
    p_out = P_S
    g_out = GAMMA_S
    # state: x, xp, y, yp for four particles (columns of a small identity)
    st = [[0.0, 0.0, 0.0, 0.0] for _ in range(4)]
    for j in range(4):
        st[j][j] = e

    samples = []  # (s, Mx, My) at the entrance, s=0
    s_cur = [0.0]

    def record():
        # st[j][i] = M[i][j] * e, particle j carries column j.
        mx = (
            (st[0][0] / e, st[1][0] / e),
            (st[0][1] / e, st[1][1] / e),
        )
        my = (
            (st[2][2] / e, st[3][2] / e),
            (st[2][3] / e, st[3][3] / e),
        )
        samples.append((s_cur[0], mx, my))

    def apply_quad(kq):
        for p in st:
            p[1] -= kq * p[0]
            p[3] += kq * p[2]
        record()

    def apply_drift(length):
        for p in st:
            p[0] += length * p[1]
            p[2] += length * p[3]
        s_cur[0] += length
        record()

    def apply_sector():
        dth = THETA / nstep
        # Enter, integrate, sample, exit. Same edge model as track_sector.
        inside = []
        for p in st:
            r = RHO + p[0]
            norm = math.sqrt(1.0 + p[1] * p[1] + p[3] * p[3])
            pr = p_out * p[1] / norm
            py = p_out * p[3] / norm
            psi, _, _ = psi_forces(r, p[2], n, c2, c4)
            g = g_out - PVS * psi
            pt = math.sqrt(max(g * g - 1.0 - pr * pr - py * py, 0.0))
            inside.append([r, pr, p[2], py, pt * r])
        for j, q in enumerate(inside):
            pt = q[4] / q[0]
            st[j][0] = q[0] - RHO
            st[j][1] = q[1] / pt
            st[j][2] = q[2]
            st[j][3] = q[3] / pt
        record()
        record_every = max(nstep // 40, 1)
        for step in range(nstep):
            for q in inside:
                q[0], q[1], q[2], q[3], _dt = rk4(
                    q[0], q[1], q[2], q[3], q[4], g_out, n, c2, c4, dth
                )
            s_cur[0] += RHO * dth
            if (step + 1) % record_every == 0 or step + 1 == nstep:
                for j, q in enumerate(inside):
                    pt = q[4] / q[0]
                    st[j][0] = q[0] - RHO
                    st[j][1] = q[1] / pt
                    st[j][2] = q[2]
                    st[j][3] = q[3] / pt
                record()
        for j, q in enumerate(inside):
            pr, py = q[1], q[3]
            pt = math.sqrt(max(p_out * p_out - pr * pr - py * py, 0.0))
            st[j][0] = q[0] - RHO
            st[j][1] = pr / pt
            st[j][2] = q[2]
            st[j][3] = py / pt
        record()

    record()
    for _c in range(N_CELL):
        apply_quad(k_quad)
        apply_drift(DRIFT)
        apply_sector()
        apply_drift(DRIFT)
        apply_quad(-k_quad)
        apply_drift(DRIFT)
        apply_sector()
        apply_drift(DRIFT)

    igx = igy = ibx = iby = 0.0
    prev = None
    for s, mx, my in samples:
        bx, ax = twiss_through(mx, beta_x, alpha_x)
        by, ay = twiss_through(my, beta_y, alpha_y)
        if bx <= 0.0 or by <= 0.0:
            raise RuntimeError("beta not positive along the cell")
        gx = (1.0 + ax * ax) / bx
        gy = (1.0 + ay * ay) / by
        if prev is not None:
            ds = s - prev[0]
            igx += 0.5 * (prev[1] + gx) * ds
            igy += 0.5 * (prev[2] + gy) * ds
            # 1/rho is nonzero only inside a sector. Drifts change s by DRIFT
            # exactly; sector samples change s by a multiple of RHO*dth.
            in_sector = abs(ds - DRIFT) > 1e-9
            if in_sector:
                ibx += 0.5 * (prev[3] + bx) / RHO**2 * ds
                iby += 0.5 * (prev[4] + by) / RHO**2 * ds
        prev = (s, gx, gy, bx, by)
    # Twiss returned to the entrance by the one-turn matrix.
    bx_f, ax_f = twiss_through(samples[-1][1], beta_x, alpha_x)
    by_f, ay_f = twiss_through(samples[-1][2], beta_y, alpha_y)
    return igx, igy, ibx, iby, bx_f, ax_f, by_f, ay_f, samples[-1][1], samples[-1][2]


def hann_average(n, c2, c4, k_quad, m_sext, nstep, v0, bx, ax, by, ay, Jx, Jy, n_turns):
    psi, psy = 0.3, 1.1
    x = v0[0] + math.sqrt(max(2.0 * Jx * bx, 0.0)) * math.cos(psi)
    xp = v0[1] - (
        math.sqrt(2.0 * Jx / bx) * (math.sin(psi) + ax * math.cos(psi)) if Jx else 0.0
    )
    y = v0[2] + math.sqrt(max(2.0 * Jy * by, 0.0)) * math.cos(psy)
    yp = v0[3] - (
        math.sqrt(2.0 * Jy / by) * (math.sin(psy) + ay * math.cos(psy)) if Jy else 0.0
    )
    t0 = track_period(*v0, 0.0, n, c2, c4, k_quad, m_sext, nstep)[4]
    acc_t = acc_x = wsum = 0.0
    err = 0.0
    for i in range(n_turns):
        x, xp, y, yp, t, xint, er = track_period(
            x, xp, y, yp, 0.0, n, c2, c4, k_quad, m_sext, nstep
        )
        w = math.sin(math.pi * (i + 0.5) / n_turns) ** 2
        acc_t += w * (t - t0)
        acc_x += w * xint
        wsum += w
        if er > err:
            err = er
    Ts = C_REF / V_S
    # (T-T_closed)/T_s . T_closed is t0; T_s is the reference period.
    # The definition uses the synchronous period. At delta=0 the closed
    # orbit is the reference, so t0 equals T_s up to the integrator error.
    lam = (acc_t / wsum) / Ts
    xint_avg = acc_x / wsum
    return lam, xint_avg, err, t0, Ts


def run_case(n, m_sext, nstep, n_turns):
    c2, c4 = field_coeffs(n)
    v0 = fixed_point(0.0, n, c2, c4, K_QUAD, m_sext, nstep)
    t_co, err0 = track_period(*v0, 0.0, n, c2, c4, K_QUAD, m_sext, nstep)[4:6]
    # [4] is time, [5] is xint; err is [6]. Recheck.
    packed = track_period(*v0, 0.0, n, c2, c4, K_QUAD, m_sext, nstep)
    t_co, err0 = packed[4], packed[6]
    Ts = C_REF / V_S
    Mx, My = one_turn_matrix(0.0, v0, n, c2, c4, K_QUAD, m_sext, nstep)
    bx, ax, nux = twiss_unfolded(Mx)
    by, ay, nuy = twiss_unfolded(My)
    # Central difference. h=1e-4 keeps the tune change above the
    # finite-difference noise of the one-turn matrix and below the
    # quadratic chromaticity.
    eps = 1e-4
    vp = fixed_point(eps, n, c2, c4, K_QUAD, m_sext, nstep)
    vm = fixed_point(-eps, n, c2, c4, K_QUAD, m_sext, nstep)
    mxp, myp = one_turn_matrix(eps, vp, n, c2, c4, K_QUAD, m_sext, nstep)
    mxm, mym = one_turn_matrix(-eps, vm, n, c2, c4, K_QUAD, m_sext, nstep)
    xi_x = (tune_unfolded(mxp) - tune_unfolded(mxm)) / (2.0 * eps)
    xi_y = (tune_unfolded(myp) - tune_unfolded(mym)) / (2.0 * eps)
    igx, igy, ibx, iby, bx2, ax2, by2, ay2, mx_s, my_s = optics_integrals(
        n, c2, c4, K_QUAD, bx, ax, by, ay, nstep
    )
    g2 = GAMMA_S**2
    # Per unit action, dimensionless Lambda.
    opt_x = 0.5 * igx / C_REF + (5.0 - n) / (2.0 * g2) * ibx / C_REF
    opt_y = 0.5 * igy / C_REF + (n - 1.0) / (2.0 * g2) * iby / C_REF
    rows = []
    for label, Jx, Jy in (("x", 4.0e-5, 0.0), ("y", 0.0, 6.0e-5)):
        lam1, x1, er1, _, _ = hann_average(
            n, c2, c4, K_QUAD, m_sext, nstep, v0, bx, ax, by, ay, Jx, Jy, n_turns
        )
        lam2, x2, er2, _, _ = hann_average(
            n, c2, c4, K_QUAD, m_sext, nstep, v0, bx, ax, by, ay, Jx / 4.0, Jy / 4.0, n_turns
        )
        lam = (16.0 * lam2 - lam1) / 3.0
        xint = (16.0 * x2 - x1) / 3.0
        J = Jx + Jy
        work = (opt_x * Jx + opt_y * Jy) + (1.0 + 1.0 / g2) / C_REF * xint
        pred = -(2.0 * math.pi / C_REF) * (Jx * xi_x + Jy * xi_y)
        rows.append((label, lam, pred, work, max(er1, er2)))
    return {
        "n": n,
        "m": m_sext,
        "nu": (nux, nuy),
        "xi": (xi_x, xi_y),
        "co": v0,
        "t_co": t_co,
        "Ts": Ts,
        "twiss_back": (bx2, ax2, by2, ay2),
        "twiss0": (bx, ax, by, ay),
        "err0": err0,
        "rows": rows,
        "opt": (opt_x, opt_y, igx, igy, ibx, iby),
    }


def main():
    nstep = N_STEP
    n_turns = N_TURNS
    if len(sys.argv) > 1 and sys.argv[1] == "quick":
        nstep = 80
        n_turns = 400
    cases = [(1.20, 0.0), (1.80, 0.0), (1.50, 0.05)]
    if len(sys.argv) > 1 and sys.argv[1] == "quick":
        cases = [(1.50, 0.0)]
    print("C", C_REF, "gamma", GAMMA_S, "steps", nstep, "turns", n_turns)
    for n, m in cases:
        out = run_case(n, m, nstep, n_turns)
        print(
            "n", out["n"], "m", out["m"],
            "nu", out["nu"], "xi", out["xi"],
            "co", [round(z, 12) for z in out["co"]],
            "t_co/Ts-1", out["t_co"] / out["Ts"] - 1.0,
            "err", out["err0"],
        )
        b0 = out["twiss0"]
        b1 = out["twiss_back"]
        print(
            "twiss closure dx",
            b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2], b1[3] - b0[3],
        )
        for label, lam, pred, work, er in out["rows"]:
            def rel(a, b):
                if abs(b) < 1e-18:
                    return float("nan")
                return (a - b) / b
            print(
                " plane", label,
                " direct", lam,
                " xi", pred,
                " work", work,
                " rel(d-xi)", rel(lam, pred),
                " rel(w-xi)", rel(work, pred),
                " rel(d-w)", rel(lam, work),
                " cons", er,
            )


if __name__ == "__main__":
    main()
