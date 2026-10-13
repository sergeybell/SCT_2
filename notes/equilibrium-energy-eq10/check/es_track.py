# -*- coding: utf-8 -*-
"""Numerical check of the electrostatic slip factor on a smooth deflector.

Units c = 1, m = 1, rho = 1. A smooth ring is one deflector closed on
itself: a model of the orbit and the revolution time, not a lattice.
The slow momentum delta is defined by the total energy at zero potential.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

RTOL = 1e-11
ATOL = 1e-13


def beta2_of(gamma):
    return 1.0 - 1.0 / gamma**2


def energy_of_delta(gamma, delta):
    return np.sqrt(1.0 + (gamma**2 - 1.0) * (1.0 + delta) ** 2)


def delta_of_energy(gamma, energy):
    # p^2 = E^2 - 1, p_s^2 = gamma^2 - 1
    p2 = energy**2 - 1.0
    ps2 = gamma**2 - 1.0
    return np.sqrt(p2 / ps2) - 1.0


def qphi(kind, gamma, x, y, z):
    b2 = beta2_of(gamma)
    if kind == "cylinder":
        r = np.hypot(x, y)
        return b2 * gamma * np.log(r)
    if kind == "sphere":
        R = np.sqrt(x * x + y * y + z * z)
        return b2 * gamma * (1.0 - 1.0 / R)
    raise ValueError(kind)


def qfield(kind, gamma, x, y, z):
    b2 = beta2_of(gamma)
    if kind == "cylinder":
        r2 = x * x + y * y
        # qE = -beta^2 gamma (x, y, 0) / r^2
        return -b2 * gamma * np.array([x / r2, y / r2, 0.0])
    if kind == "sphere":
        R2 = x * x + y * y + z * z
        R3 = R2 * np.sqrt(R2)
        return -b2 * gamma * np.array([x, y, z]) / R3
    raise ValueError(kind)


def pack(x, y, z, px, py, pz):
    phi = np.arctan2(y, x)
    return np.array([x, y, z, px, py, pz, phi], dtype=float)


def deriv(kind, gamma):
    def f(_t, state):
        x, y, z, px, py, pz, _phi = state
        g = np.sqrt(1.0 + px * px + py * py + pz * pz)
        vx, vy, vz = px / g, py / g, pz / g
        ex, ey, ez = qfield(kind, gamma, x, y, z)
        r2 = x * x + y * y
        dphi = (x * vy - y * vx) / r2
        return (vx, vy, vz, ex, ey, ez, dphi)

    return f


def radial_residual(kind, gamma, energy, radius):
    """Azimuthal launch at (r,0,0): centripetal imbalance, zero on a circle."""
    phi_pot = qphi(kind, gamma, radius, 0.0, 0.0)
    gkin = energy - phi_pot
    if gkin <= 1.0 + 1e-12:
        return 1e3
    p = np.sqrt(gkin * gkin - 1.0)
    # v = (0, p/gkin, 0), qE_x = qfield(...)[0]
    ex = qfield(kind, gamma, radius, 0.0, 0.0)[0]
    ax = ex / gkin  # p_x = 0, so dv_x/dt = (dp_x/dt)/gamma
    v = p / gkin
    return ax + v * v / radius


def circular_radius(kind, gamma, delta):
    energy = energy_of_delta(gamma, delta)
    grid = np.linspace(0.55, 1.8, 40)
    vals = [radial_residual(kind, gamma, energy, float(r)) for r in grid]
    bracket = None
    for r0, r1, f0, f1 in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if f0 == 0.0:
            return float(r0), energy
        if np.isfinite(f0) and np.isfinite(f1) and f0 * f1 < 0.0:
            bracket = (float(r0), float(r1))
            break
    if bracket is None:
        raise RuntimeError(
            f"no bracket {kind} gamma={gamma} delta={delta} {list(zip(grid, vals))[:6]}"
        )
    root = root_scalar(
        lambda r: radial_residual(kind, gamma, energy, r),
        bracket=bracket,
        xtol=1e-14,
    )
    if not root.converged:
        raise RuntimeError("radius search failed")
    return float(root.root), energy


def launch_circular(kind, gamma, delta):
    radius, energy = circular_radius(kind, gamma, delta)
    gkin = energy - qphi(kind, gamma, radius, 0.0, 0.0)
    p = np.sqrt(gkin * gkin - 1.0)
    state = pack(radius, 0.0, 0.0, 0.0, p, 0.0)
    return state, radius, gkin


def one_turn(kind, gamma, state, t_hint):
    """Integrate until the cylindrical azimuth advances by 2 pi."""

    def event(_t, y):
        return y[6] - 2.0 * np.pi

    event.terminal = True
    event.direction = 1
    sol = solve_ivp(
        deriv(kind, gamma),
        (0.0, t_hint * 1.8),
        state,
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        dense_output=True,
        events=event,
    )
    if sol.status != 1 or len(sol.t_events[0]) == 0:
        raise RuntimeError(f"turn not closed, status={sol.status}")
    t_turn = float(sol.t_events[0][0])
    y_end = sol.sol(t_turn)
    return t_turn, y_end, sol


def reference_period(gamma):
    return 2.0 * np.pi / np.sqrt(beta2_of(gamma))


def fit_slip(kind, gamma, deltas):
    ts = reference_period(gamma)
    taus = []
    radii = []
    for delta in deltas:
        state, radius, _gkin = launch_circular(kind, gamma, delta)
        # hint: path ~ 2 pi r, speed ~ beta
        t_turn, _y, _sol = one_turn(kind, gamma, state, ts * max(radius, 0.5))
        taus.append(t_turn / ts - 1.0)
        radii.append(radius)
    coef = np.polyfit(deltas, taus, 4)  # c4 d^4 + ... + c1 d + c0
    # np.polyfit returns highest degree first
    eta0 = coef[-2]
    eta1 = coef[-3]
    return eta0, eta1, np.array(taus), np.array(radii)


def formula_slip(kind, gamma):
    g2 = gamma**2
    b2 = beta2_of(gamma)
    if kind == "cylinder":
        return 1.0, 0.5 + 0.5 / g2
    if kind == "sphere":
        return g2 + 2.0, 4.5 * g2 + 3.5 + 1.0 / g2
    raise ValueError(kind)


def sample_turn(sol, t_turn, nphi=360):
    """Azimuthal samples of x, y and slopes on one revolution."""
    phis = np.linspace(0.0, 2.0 * np.pi, nphi, endpoint=False)
    # phi(t) is monotone; invert on a time grid
    ts = np.linspace(0.0, t_turn, nphi * 8)
    ys = sol.sol(ts)
    phi = ys[6]
    # guard against a tiny overshoot
    phi = np.clip(phi, 0.0, 2.0 * np.pi)
    x = np.interp(phis, phi, np.hypot(ys[0], ys[1]) - 1.0)
    y = np.interp(phis, phi, ys[2])
    # slopes from the state at the same times, then interpolated
    g = np.sqrt(1.0 + ys[3] ** 2 + ys[4] ** 2 + ys[5] ** 2)
    vx, vy, vz = ys[3] / g, ys[4] / g, ys[5] / g
    r2 = ys[0] ** 2 + ys[1] ** 2
    dphi = (ys[0] * vy - ys[1] * vx) / r2
    dr = (ys[0] * vx + ys[1] * vy) / np.sqrt(r2)
    xp = np.interp(phis, phi, dr / dphi)
    yp = np.interp(phis, phi, vz / dphi)
    return x, y, xp, yp


def horizontal_tune(kind, gamma):
    if kind == "cylinder":
        return np.sqrt(1.0 + 1.0 / gamma**2)
    if kind == "sphere":
        return 1.0 / gamma
    raise ValueError(kind)


def _fit_constant(phase, values, harmonics=(1, 2)):
    """Constant term of a harmonic fit. Falls back to the mean if a column repeats it."""
    cols = [np.ones_like(phase)]
    for k in harmonics:
        cols.append(np.cos(k * phase))
        cols.append(np.sin(k * phase))
    design = np.column_stack(cols)
    # Drop columns that duplicate the constant (integer tune).
    keep = [0]
    for j in range(1, design.shape[1]):
        if np.std(design[:, j] - design[:, 0]) > 1e-8:
            keep.append(j)
    coef, *_ = np.linalg.lstsq(design[:, keep], values, rcond=None)
    return float(coef[0])


def betatron_scan(kind, gamma, mode, amplitude, nturns=40):
    """Mean tau and azimuthal moments at delta = 0.

    The turn-by-turn period still carries the betatron oscillation. Its
    constant term is read by a harmonic fit at the known tune, so a
    non-integer tune does not bias the average. Moments of a horizontal
    orbit are read the same way from x(theta).
    """
    b = np.sqrt(beta2_of(gamma))
    p = b * gamma
    if mode == "horizontal":
        pr = amplitude * p
        pt = np.sqrt(p * p - pr * pr)
        state = pack(1.0, 0.0, 0.0, pr, pt, 0.0)
        nu = horizontal_tune(kind, gamma)
    elif mode == "vertical":
        pt = p * np.cos(amplitude)
        pz = p * np.sin(amplitude)
        state = pack(1.0, 0.0, 0.0, 0.0, pt, pz)
        nu = 1.0
    else:
        raise ValueError(mode)
    ts = reference_period(gamma)
    taus = []
    chunks_x, chunks_y, chunks_xp, chunks_yp = [], [], [], []
    for _ in range(nturns):
        t_turn, y_end, sol = one_turn(kind, gamma, state, ts * 1.5)
        taus.append(t_turn / ts - 1.0)
        x, y, xp, yp = sample_turn(sol, t_turn, nphi=180)
        chunks_x.append(x)
        chunks_y.append(y)
        chunks_xp.append(xp)
        chunks_yp.append(yp)
        state = y_end.copy()
        state[6] = 0.0
    taus = np.asarray(taus)
    tau = _fit_constant(2.0 * np.pi * nu * np.arange(nturns), taus)
    x = np.concatenate(chunks_x)
    y = np.concatenate(chunks_y)
    xp = np.concatenate(chunks_xp)
    yp = np.concatenate(chunks_yp)
    theta = np.linspace(0.0, 2.0 * np.pi * nturns, x.size, endpoint=False)
    if mode == "horizontal":
        x_mean = _fit_constant(nu * theta, x)
        design = np.column_stack(
            [np.cos(nu * theta), np.sin(nu * theta), np.ones_like(theta)]
        )
        coef, *_ = np.linalg.lstsq(design, x, rcond=None)
        x2 = 0.5 * (coef[0] ** 2 + coef[1] ** 2)
        xp2 = nu**2 * x2
        y2 = float(np.mean(y * y))
        yp2 = float(np.mean(yp * yp))
        moments = {"x": x_mean, "x2": x2, "y2": y2, "xp2": xp2, "yp2": yp2}
    else:
        moments = {
            "x": float(np.mean(x)),
            "x2": float(np.mean(x * x)),
            "y2": float(np.mean(y * y)),
            "xp2": float(np.mean(xp * xp)),
            "yp2": float(np.mean(yp * yp)),
        }
    return tau, moments


def lambda_from_moments(gamma, n, moments):
    b2 = beta2_of(gamma)
    g2 = gamma**2
    return (
        (2.0 - b2) * moments["x"]
        + 0.5 * (moments["xp2"] + moments["yp2"])
        + (5.0 - n) / (2.0 * g2) * moments["x2"]
        + (n - 1.0) / (2.0 * g2) * moments["y2"]
    )


def path_from_moments(moments):
    return moments["x"] + 0.5 * (moments["xp2"] + moments["yp2"])


def main():
    # Small deltas keep the degree-4 fit free of the delta^5 tail.
    # The same fit on the algebraic sphere period shows the tail, not the integrator.
    deltas = np.array(
        [-0.002, -0.001, -0.0005, -0.00025, 0.00025, 0.0005, 0.001, 0.002]
    )
    rows = []
    for gamma in (1.25, 2.0):
        for kind in ("cylinder", "sphere"):
            eta0_f, eta1_f = formula_slip(kind, gamma)
            eta0, eta1, _taus, _radii = fit_slip(kind, gamma, deltas)
            rows.append(
                (
                    kind,
                    gamma,
                    eta0_f,
                    eta0,
                    eta0 - eta0_f,
                    eta1_f,
                    eta1,
                    eta1 - eta1_f,
                )
            )
            print(
                f"slip {kind:8s} gamma={gamma:.2f}  "
                f"eta0 {eta0:.10f} formula {eta0_f:.10f}  "
                f"d={eta0 - eta0_f:.3e}  "
                f"eta1 {eta1:.10f} formula {eta1_f:.10f}  "
                f"d={eta1 - eta1_f:.3e}"
            )
            # Residuals at this level are the delta^5 tail of the degree-4 fit.
            if abs(eta0 - eta0_f) > 1e-8 or abs(eta1 - eta1_f) > 5e-8:
                raise SystemExit(f"slip mismatch {kind} {gamma}")

    print("--- betatron ---")
    # gamma = 5/4 makes the spherical tune 4/5 and the cylindrical tune
    # irrational; the constant term is taken from the harmonic fit.
    gamma = 1.25
    g2 = gamma**2
    checks = []

    tau, mom = betatron_scan("cylinder", gamma, "horizontal", 0.01, nturns=50)
    pred = lambda_from_moments(gamma, 1, mom)
    path = path_from_moments(mom)
    x2 = mom["x2"]
    x_ratio = mom["x"] / x2
    x_target = -(g2 + 5) / (2 * (g2 + 1))
    path_ratio = path / x2
    path_target = (1 - 3 * g2) / (2 * g2 * (g2 + 1))
    tau_ratio = tau / x2
    print(
        f"cylinder horizontal  x/x2={x_ratio:.6f} ({x_target:.6f})  "
        f"tau/x2={tau_ratio:.3e}  path/x2={path_ratio:.6f} ({path_target:.6f})  "
        f"formula/x2={pred/x2:.3e}"
    )
    if abs(x_ratio - x_target) > 2e-3:
        raise SystemExit("cylinder orbit center")
    if abs(tau_ratio) > 5e-3:
        raise SystemExit("cylinder time did not cancel")
    if abs(path_ratio - path_target) > 5e-3:
        raise SystemExit("cylinder path")
    if abs((pred - tau) / x2) > 5e-3:
        raise SystemExit("cylinder formula vs time")
    checks.append(("cylinder", "horizontal", tau_ratio, path_ratio, x_ratio))

    tau, mom = betatron_scan("sphere", gamma, "horizontal", 0.01, nturns=40)
    pred = lambda_from_moments(gamma, 2, mom)
    path = path_from_moments(mom)
    x2 = mom["x2"]
    tau_ratio = tau / x2
    path_ratio = path / x2
    x_ratio = mom["x"] / x2
    print(
        f"sphere horizontal    x/x2={x_ratio:.6f} (-1)  "
        f"tau/x2={tau_ratio:.6f} ({-beta2_of(gamma):.6f})  "
        f"path/x2={path_ratio:.6f} ({-1 + 0.5/g2:.6f})  "
        f"formula/x2={pred/x2:.6f}"
    )
    if abs(x_ratio + 1) > 2e-3 or abs(tau_ratio + beta2_of(gamma)) > 2e-3:
        raise SystemExit("sphere horizontal")
    if abs(path_ratio - (-1 + 0.5 / g2)) > 2e-3:
        raise SystemExit("sphere path")
    if abs(pred / x2 - tau_ratio) > 2e-3:
        raise SystemExit("sphere formula vs time")
    checks.append(("sphere", "horizontal", tau_ratio, path_ratio, x_ratio))

    tau, mom = betatron_scan("sphere", gamma, "vertical", 0.02, nturns=20)
    y2 = mom["y2"]
    print(
        f"sphere vertical      x/y2={mom['x']/y2:.6f} (-0.5)  "
        f"tau/y2={tau/y2:.3e}  tau={tau:.3e}"
    )
    if abs(mom["x"] / y2 + 0.5) > 2e-3 or abs(tau / y2) > 1e-4:
        raise SystemExit("sphere vertical")
    checks.append(("sphere", "vertical", tau / y2, mom["x"] / y2, y2))

    print("all tracking checks passed")
    return rows, checks


if __name__ == "__main__":
    main()
