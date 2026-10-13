"""Symplectic identity for a one-turn map of thin Hamiltonian kicks.

The pair (tau, delta) is canonical: d(tau)/dJ = d(mu)/d(delta) = 2*pi*d(nu)/d(delta).
A second copy of the same ring drops the longitudinal half of each kick and is
not symplectic. Prints the ratio of the two sides, a sextupole scan, the
closed-orbit time against the expansion at x=0, and the Twiss quadratic form
at two observation points.
"""

from __future__ import annotations

import math

import numpy as np

# One plane. State (x, p, tau, delta). delta is the momentum conjugate to tau.
N_CELL = 6
L_DRIFT = 0.8
L_BEND = 0.6
H_BEND = 0.15
K_QUAD = 0.55
J_PROBE = 1.0e-5
N_PHASE = 32


def bend_or_drift(L, h, x, p, tau, delta):
    """Flow of H = p^2/(2(1+delta)) - h*x*delta over length L."""
    d = delta
    den = 1.0 + d
    pL = p + h * d * L
    xL = x + p * L / den + h * d * L * L / (2.0 * den)
    int_p2 = p * p * L + p * h * d * L * L + (h * d) ** 2 * L**3 / 3.0
    int_x = x * L + p * L * L / (2.0 * den) + h * d * L**3 / (6.0 * den)
    tauL = tau - int_p2 / (2.0 * den * den) - h * int_x
    return xL, pL, tauL, d


def kick(strength, power, x, p, tau, delta, keep_tau):
    """Time-1 map of f = strength * x^power / (power! * (1+delta)).

    power=2 is a quadrupole, power=3 a sextupole. The transverse kick is
    Delta p = -df/dx. The partner Delta tau = df/ddelta is required for
    symplecticity and is omitted when keep_tau is false.
    """
    den = 1.0 + delta
    fact = math.factorial(power)
    f_scale = strength / (fact * den)
    # df/dx = strength * x^(power-1) / ((power-1)! * den)
    dp = -strength * x ** (power - 1) / (math.factorial(power - 1) * den)
    dt = 0.0
    if keep_tau:
        # df/ddelta = - strength * x^power / (power! * den^2)
        dt = -strength * x**power / (fact * den * den)
    return x, p + dp, tau + dt, delta


def elements(sextupole):
    """FODO cell, repeated. Sextupole sits next to the focusing quad."""
    cell = [
        ("kick", 2, K_QUAD),
        ("kick", 3, sextupole),
        ("line", L_DRIFT, 0.0),
        ("line", L_BEND, H_BEND),
        ("kick", 2, -K_QUAD),
        ("line", L_DRIFT, 0.0),
        ("line", L_BEND, H_BEND),
    ]
    ring = []
    for _ in range(N_CELL):
        ring.extend(cell)
    return ring


def propagate(x, p, tau, delta, ring, keep_tau, start):
    n = len(ring)
    for k in range(n):
        kind, a, b = ring[(start + k) % n]
        if kind == "line":
            x, p, tau, delta = bend_or_drift(a, b, x, p, tau, delta)
        else:
            x, p, tau, delta = kick(b, a, x, p, tau, delta, keep_tau)
    return x, p, tau, delta


def jacobian(delta, ring, keep_tau, start):
    eps = 1e-7
    base = np.array(propagate(0.0, 0.0, 0.0, delta, ring, keep_tau, start))
    cols = []
    for i, step in enumerate((eps, eps, eps, eps)):
        z = [0.0, 0.0, 0.0, delta]
        z[i] += step
        out = np.array(propagate(z[0], z[1], z[2], z[3], ring, keep_tau, start))
        # delta is invariant; the fourth column is a derivative at fixed input delta
        cols.append((out - base) / step)
    return np.column_stack(cols)


def symplectic_defect(mat):
    s = np.array(
        [
            [0.0, 1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0, 0.0],
        ]
    )
    return float(np.max(np.abs(mat.T @ s @ mat - s)))


def fixed_point(delta, ring, keep_tau, start):
    x, p = 0.0, 0.0
    eps = 1e-8
    for _ in range(12):
        x0, p0, _, _ = propagate(x, p, 0.0, delta, ring, keep_tau, start)
        rx, rp = x0 - x, p0 - p
        xx, px, _, _ = propagate(x + eps, p, 0.0, delta, ring, keep_tau, start)
        xy, py, _, _ = propagate(x, p + eps, 0.0, delta, ring, keep_tau, start)
        jac = np.array([[(xx - x0) / eps, (xy - x0) / eps], [(px - p0) / eps, (py - p0) / eps]])
        step = np.linalg.solve(jac - np.eye(2), -np.array([rx, rp]))
        x += float(step[0])
        p += float(step[1])
        if abs(rx) + abs(rp) < 1e-14:
            break
    return x, p


def transverse_matrix(delta, ring, keep_tau, start, x_fp, p_fp):
    eps = 1e-7
    x0, p0, _, _ = propagate(x_fp, p_fp, 0.0, delta, ring, keep_tau, start)
    xx, px, _, _ = propagate(x_fp + eps, p_fp, 0.0, delta, ring, keep_tau, start)
    xy, py, _, _ = propagate(x_fp, p_fp + eps, 0.0, delta, ring, keep_tau, start)
    return np.array([[(xx - x0) / eps, (xy - x0) / eps], [(px - p0) / eps, (py - p0) / eps]])


def tune_twiss(mat):
    cos_mu = float(np.trace(mat) / 2.0)
    if abs(cos_mu) > 1.0 - 1e-12:
        raise RuntimeError(f"unstable or integer, cos={cos_mu}")
    sin_sign = 1.0 if mat[0, 1] >= 0.0 else -1.0
    sin_mu = sin_sign * math.sqrt(max(0.0, 1.0 - cos_mu * cos_mu))
    mu = math.atan2(sin_mu, cos_mu)
    if mu < 0.0:
        mu += 2.0 * math.pi
    beta = mat[0, 1] / sin_mu
    alpha = (mat[0, 0] - mat[1, 1]) / (2.0 * sin_mu)
    gamma = (1.0 + alpha * alpha) / beta
    return mu / (2.0 * math.pi), beta, alpha, gamma


def ellipse_average(delta, ring, keep_tau, start, action):
    x_fp, p_fp = fixed_point(delta, ring, keep_tau, start)
    nu, beta, alpha, gamma = tune_twiss(
        transverse_matrix(delta, ring, keep_tau, start, x_fp, p_fp)
    )
    amp_x = math.sqrt(2.0 * action * beta)
    amp_p = math.sqrt(2.0 * action / beta)
    dtaus = []
    xs = []
    ps = []
    for phase in np.linspace(0.0, 2.0 * math.pi, N_PHASE, endpoint=False):
        xx = x_fp + amp_x * math.cos(phase)
        pp = p_fp - amp_p * (math.sin(phase) + alpha * math.cos(phase))
        xf, pf, tau_f, _ = propagate(xx, pp, 0.0, delta, ring, keep_tau, start)
        dtaus.append(tau_f)
        xs.append(xx)
        ps.append(pp)
    return {
        "nu": nu,
        "beta": beta,
        "alpha": alpha,
        "gamma": gamma,
        "dtau": float(np.mean(dtaus)),
        "x_mean": float(np.mean(xs)),
        "p_mean": float(np.mean(ps)),
        "x_fp": x_fp,
        "p_fp": p_fp,
    }


def quadratic_form(ring, keep_tau, start):
    """Mononial coefficients of Delta tau at the origin, times Twiss."""
    eps = 2e-5

    def dtau(x, p):
        return propagate(x, p, 0.0, 0.0, ring, keep_tau, start)[2]

    c11 = (dtau(eps, 0.0) + dtau(-eps, 0.0) - 2.0 * dtau(0.0, 0.0)) / (2.0 * eps * eps)
    c22 = (dtau(0.0, eps) + dtau(0.0, -eps) - 2.0 * dtau(0.0, 0.0)) / (2.0 * eps * eps)
    c12 = (
        dtau(eps, eps) - dtau(eps, -eps) - dtau(-eps, eps) + dtau(-eps, -eps)
    ) / (4.0 * eps * eps)
    info = ellipse_average(0.0, ring, keep_tau, start, J_PROBE)
    form = c11 * info["beta"] - c12 * info["alpha"] + c22 * info["gamma"]
    # Linear time response, for the centre shift.
    r51 = (dtau(eps, 0.0) - dtau(-eps, 0.0)) / (2.0 * eps)
    r52 = (dtau(0.0, eps) - dtau(0.0, -eps)) / (2.0 * eps)
    return {
        "c11": c11,
        "c12": c12,
        "c22": c22,
        "form": form,
        "r51": r51,
        "r52": r52,
        "info": info,
    }


def chromaticity(ring, keep_tau, start):
    deltas = np.array([-2e-4, -1e-4, 1e-4, 2e-4])
    nus = []
    for d in deltas:
        nus.append(ellipse_average(float(d), ring, keep_tau, start, J_PROBE)["nu"])
    design = np.column_stack([np.ones(4), deltas, deltas**2])
    coef = np.linalg.lstsq(design, np.array(nus), rcond=None)[0]
    return float(coef[1])


def origin_vs_closed(ring, keep_tau):
    """d(Delta tau)/d(delta) at x=p=0 and on the closed orbit."""
    eps = 1e-5

    def at_origin(d):
        return propagate(0.0, 0.0, 0.0, d, ring, keep_tau, 0)[2]

    def on_orbit(d):
        x, p = fixed_point(d, ring, keep_tau, 0)
        return propagate(x, p, 0.0, d, ring, keep_tau, 0)[2]

    d_origin = (at_origin(eps) - at_origin(-eps)) / (2.0 * eps)
    d_closed = (on_orbit(eps) - on_orbit(-eps)) / (2.0 * eps)
    return d_origin, d_closed


def analyse(sextupole, keep_tau, start=0):
    ring = elements(sextupole)
    info = ellipse_average(0.0, ring, keep_tau, start, J_PROBE)
    xi = chromaticity(ring, keep_tau, start)
    pred = 2.0 * math.pi * xi
    slope = info["dtau"] / J_PROBE
    quad = quadratic_form(ring, keep_tau, start)
    # Invariant-curve centre from the one-turn transverse kick averaged on the
    # linear ellipse: (I - R) (x_c, p_c) = <M(z) - z>.
    mat = transverse_matrix(0.0, ring, keep_tau, start, info["x_fp"], info["p_fp"])
    amp_x = math.sqrt(2.0 * J_PROBE * info["beta"])
    amp_p = math.sqrt(2.0 * J_PROBE / info["beta"])
    dxs, dps = [], []
    for phase in np.linspace(0.0, 2.0 * math.pi, N_PHASE, endpoint=False):
        xx = info["x_fp"] + amp_x * math.cos(phase)
        pp = info["p_fp"] - amp_p * (math.sin(phase) + info["alpha"] * math.cos(phase))
        xf, pf, _, _ = propagate(xx, pp, 0.0, 0.0, ring, keep_tau, start)
        dxs.append(xf - xx)
        dps.append(pf - pp)
    kick = np.array([float(np.mean(dxs)), float(np.mean(dps))])
    centre_vec = np.linalg.solve(np.eye(2) - mat, kick)
    centre = quad["r51"] * centre_vec[0] + quad["r52"] * centre_vec[1]
    return {
        "xi": xi,
        "pred": pred,
        "slope": slope,
        "ratio": slope / pred if abs(pred) > 1e-14 else float("nan"),
        "form": quad["form"],
        "form_ratio": quad["form"] / pred if abs(pred) > 1e-14 else float("nan"),
        "centre": centre / J_PROBE,
        "form_plus_centre": quad["form"] + centre / J_PROBE,
        "defect": symplectic_defect(jacobian(0.0, ring, keep_tau, start)),
        "nu": info["nu"],
        "beta": info["beta"],
        "alpha": info["alpha"],
    }


def scan(keep_tau):
    strengths = [0.0, 0.05, 0.1, 0.2, 0.35]
    rows = []
    for s in strengths:
        row = analyse(s, keep_tau, start=0)
        row["s"] = s
        rows.append(row)
    pred = np.array([r["pred"] for r in rows])
    slope = np.array([r["form_plus_centre"] for r in rows])
    bare = np.array([r["form"] for r in rows])
    design = np.column_stack([pred, np.ones(len(rows))])
    kappa, b = np.linalg.lstsq(design, slope, rcond=None)[0]
    kappa_b, b_b = np.linalg.lstsq(design, bare, rcond=None)[0]
    return rows, float(kappa), float(b), float(kappa_b), float(b_b)


def main():
    ring = elements(0.0)
    print(f"symplectic defect, kicks kept: {symplectic_defect(jacobian(0.0, ring, True, 0)):.3e}")
    print(f"symplectic defect, tau dropped: {symplectic_defect(jacobian(0.0, ring, False, 0)):.3e}")
    d_ori, d_clo = origin_vs_closed(ring, True)
    print(f"d(tau)/d(delta) at x=0: {d_ori:.8e}")
    print(f"d(tau)/d(delta) on closed orbit: {d_clo:.8e}")
    print(f"relative difference: {(d_ori - d_clo) / d_clo:.6e}")

    print("--- symplectic, two observation points, S=0.2 ---")
    for start, name in ((0, "focus"), (4, "defocus")):
        row = analyse(0.2, True, start=start)
        print(
            f"{name}: nu={row['nu']:.6f} xi={row['xi']:.6f} "
            f"form/pred={row['form_ratio']:.8f} "
            f"sum/pred={row['form_plus_centre'] / row['pred']:.8f} "
            f"beta={row['beta']:.4f}"
        )

    for label, keep in (("symplectic", True), ("tau_dropped", False)):
        rows, kappa, b, kappa_b, b_b = scan(keep)
        print(
            f"--- scan {label}: with centre kappa={kappa:.6f} b={b:.3e}; "
            f"ellipse only kappa={kappa_b:.6f} b={b_b:.3e} ---"
        )
        for r in rows:
            print(
                f"S={r['s']:.2f} xi={r['xi']:.6f} pred={r['pred']:.6e} "
                f"slope={r['slope']:.6e} ratio={r['ratio']:.6f} "
                f"form={r['form']:.6e}"
            )


if __name__ == "__main__":
    main()
