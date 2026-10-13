# -*- coding: utf-8 -*-
"""Numerical check of the longitudinal formulas against the full system.

Time is measured in units of 1/omega_s.  The integrator is DOP853.
The average is the uniform average over one synchrotron period in this time.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import ellipk

# ---------------------------------------------------------------------------
# Formulas under test
# ---------------------------------------------------------------------------

def predict(eta0, eta1, Lam, rho, a):
    """Second-order predictions.  a is the first-harmonic amplitude."""
    deq = -(eta1 * a**2 / 2.0 + Lam) / eta0
    phi_m = rho * abs(eta0) * a
    domega = (
        -(eta1 * Lam) / eta0**2
        - (5.0 / 12.0) * (eta1 * a / eta0) ** 2
        - phi_m**2 / 16.0
    )
    u2 = eta1 * a**2 / (6.0 * eta0)
    # O(a**4) piece of <delta> at cos(phi)=1, Lambda=0 (19/72 term)
    deq_a4 = deq - (19.0 / 72.0) * eta1**3 * a**4 / eta0**3
    return {
        "deq": deq,
        "domega": domega,
        "phi_m": phi_m,
        "u2": u2,
        "deq_a4": deq_a4,
    }


def pendulum_domega(phi_max):
    """Exact (Omega-1) for phi'' + sin phi = 0, amplitude phi_max."""
    if phi_max == 0.0:
        return 0.0
    k2 = np.sin(phi_max / 2.0) ** 2
    ratio = np.pi / (2.0 * ellipk(k2))
    return ratio - 1.0


# ---------------------------------------------------------------------------
# Full system
# ---------------------------------------------------------------------------

def integrate_full(eta0, eta1, Lam, rho, a_guess, n_periods=6, rtol=1e-12, atol=1e-14):
    def f(_t, y):
        phi, delta = y
        tau = eta0 * delta + eta1 * delta * delta + Lam
        return np.array([rho * tau, -np.sin(phi) / (rho * eta0)], dtype=float)

    direction = 1.0 if eta0 > 0.0 else -1.0

    def phi_cross(_t, y):
        return y[0]

    phi_cross.direction = direction
    phi_cross.terminal = False

    deq = -(eta1 * a_guess**2 / 2.0 + Lam) / eta0
    u2 = eta1 * a_guess**2 / (6.0 * eta0)
    delta0 = a_guess + deq + u2
    y0 = np.array([0.0, delta0], dtype=float)

    # rough period ~ 2*pi; integrate several of them
    t_end = n_periods * 2.0 * np.pi * 1.4
    sol = solve_ivp(
        f,
        (0.0, t_end),
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=True,
        events=phi_cross,
        max_step=0.02,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    roots = [float(t) for t in sol.t_events[0] if t > 1e-6]

    def refine(t_guess):
        half = 5e-3
        lo, hi = t_guess - half, t_guess + half
        # shrink until a sign change
        flo, fhi = sol.sol(lo)[0], sol.sol(hi)[0]
        if flo * fhi > 0.0:
            # fall back to the event time
            return t_guess
        return float(brentq(lambda t: sol.sol(t)[0], lo, hi, xtol=1e-14))

    roots = [refine(t) for t in roots]
    if len(roots) < 3:
        raise RuntimeError(f"only {len(roots)} returns, eta0={eta0}")

    # use a late period so the local truncation has not accumulated from a bad start
    t0, t1 = roots[-3], roots[-2]
    return _measure(sol, t0, t1, eta0, eta1, Lam, rho)


def integrate_reduced(eta0, eta1, Lam, a_guess, n_periods=6, rtol=1e-12, atol=1e-14):
    """cos phi = 1:  delta'' + delta = c0 + c2 delta**2, c0=-Lam/eta0, c2=-eta1/eta0."""
    c0 = -Lam / eta0
    c2 = -eta1 / eta0

    def f(_t, y):
        delta, v = y
        return np.array([v, -delta + c0 + c2 * delta * delta], dtype=float)

    def v_cross(_t, y):
        return y[1]

    # start at a maximum: v goes from 0 toward negative, next maximum has direction -1
    v_cross.direction = -1.0
    v_cross.terminal = False

    delta0 = a_guess + c0 + c2 * a_guess**2 / 3.0
    y0 = np.array([delta0, 0.0], dtype=float)
    t_end = n_periods * 2.0 * np.pi * 1.4
    sol = solve_ivp(
        f,
        (0.0, t_end),
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=True,
        events=v_cross,
        max_step=0.02,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    roots = [float(t) for t in sol.t_events[0] if t > 1e-6]

    def refine(t_guess):
        half = 5e-3
        lo, hi = t_guess - half, t_guess + half
        if sol.sol(lo)[1] * sol.sol(hi)[1] > 0.0:
            return t_guess
        return float(brentq(lambda t: sol.sol(t)[1], lo, hi, xtol=1e-14))

    roots = [refine(t) for t in roots]
    if len(roots) < 3:
        raise RuntimeError(f"reduced: only {len(roots)} returns")
    t0, t1 = roots[-3], roots[-2]
    return _measure_reduced(sol, t0, t1, eta0, eta1, Lam)


def _fourier_periodic(ts, values, n_harm):
    """Rectangle rule on a uniform periodic grid, endpoint excluded."""
    n = len(ts)
    dt_over_T = 1.0 / n
    acc = np.zeros(n_harm + 1, dtype=float)
    # values are sampled at ts, but the phase is 2*pi*i/n
    idx = np.arange(n)
    mean = float(np.sum(values) * dt_over_T)
    out = {"mean": mean}
    for k in range(1, n_harm + 1):
        ang = 2.0 * np.pi * k * idx / n
        ck = 2.0 * float(np.sum(values * np.cos(ang)) * dt_over_T)
        sk = 2.0 * float(np.sum(values * np.sin(ang)) * dt_over_T)
        out[k] = (ck, sk)
    return out


def _sample(sol, t0, t1, n, component):
    # uniform grid, exclude endpoint so the sample is periodic
    ts = t0 + (t1 - t0) * np.arange(n) / n
    ys = sol.sol(ts)
    return ts, ys[component]


def _measure(sol, t0, t1, eta0, eta1, Lam, rho, n=4096):
    _ts, phi = _sample(sol, t0, t1, n, 0)
    _ts, delta = _sample(sol, t0, t1, n, 1)
    T = t1 - t0
    fd = _fourier_periodic(_ts, delta, 4)
    a_cos, a_sin = fd[1]
    a = float(np.hypot(a_cos, a_sin))
    c2, s2 = fd[2]
    # second-harmonic cosine in the frame where the first harmonic is pure cosine.
    # phase of the first harmonic: delta ~ a cos(theta - phase) if we use atan2
    tau = eta0 * delta + eta1 * delta * delta + Lam
    tau_mean = float(np.sum(tau) / n)
    d2_mean = float(np.sum(delta * delta) / n)
    phi_max = float(np.max(np.abs(phi)))
    # c1 = a cos alpha, s1 = a sin alpha; alpha ~ 0 when t0 is a maximum of delta
    alpha = np.arctan2(a_sin, a_cos)
    # value of cos(2(theta-alpha)) coefficient
    u2 = c2 * np.cos(2 * alpha) + s2 * np.sin(2 * alpha)
    return {
        "T": T,
        "Omega": 2.0 * np.pi / T,
        "mean": fd["mean"],
        "a": a,
        "a_sin": a_sin,
        "u2": float(u2),
        "tau_mean": tau_mean,
        "d2_mean": d2_mean,
        "identity": eta0 * fd["mean"] + eta1 * d2_mean + Lam,
        "phi_max": phi_max,
        "phi_m_lin": rho * abs(eta0) * a,
    }


def _measure_reduced(sol, t0, t1, eta0, eta1, Lam, n=4096):
    _ts, delta = _sample(sol, t0, t1, n, 0)
    T = t1 - t0
    fd = _fourier_periodic(_ts, delta, 4)
    a_cos, a_sin = fd[1]
    a = float(np.hypot(a_cos, a_sin))
    c2, s2 = fd[2]
    alpha = np.arctan2(a_sin, a_cos)
    u2 = c2 * np.cos(2.0 * alpha) + s2 * np.sin(2.0 * alpha)
    d2_mean = float(np.sum(delta * delta) / n)
    return {
        "T": T,
        "Omega": 2.0 * np.pi / T,
        "mean": fd["mean"],
        "a": a,
        "u2": float(u2),
        "d2_mean": d2_mean,
        "identity": eta0 * fd["mean"] + eta1 * d2_mean + Lam,
    }


def report_case(tag, eta0, eta1, Lam, rho, a_guess):
    m = integrate_full(eta0, eta1, Lam, rho, a_guess)
    p = predict(eta0, eta1, Lam, rho, m["a"])
    dom_num = m["Omega"] - 1.0
    line = {
        "tag": tag,
        "eta0": eta0,
        "eta1": eta1,
        "Lam": Lam,
        "rho": rho,
        "a": m["a"],
        "mean": m["mean"],
        "deq": p["deq"],
        "deq_a4": p["deq_a4"],
        "err_mean": m["mean"] - p["deq"],
        "err_mean_a4": m["mean"] - p["deq_a4"],
        "domega": dom_num,
        "domega_pred": p["domega"],
        "err_omega": dom_num - p["domega"],
        "u2": m["u2"],
        "u2_pred": p["u2"],
        "identity": m["identity"],
        "phi_max": m["phi_max"],
        "phi_m": p["phi_m"],
        "tau_mean": m["tau_mean"],
    }
    if eta1 == 0.0 and Lam == 0.0:
        line["pendulum"] = pendulum_domega(m["phi_max"])
        line["err_pend"] = dom_num - line["pendulum"]
    return line


def fmt(x, sig=6):
    if x == 0 or not np.isfinite(x):
        return f"{x}"
    ax = abs(x)
    if ax != 0 and (ax < 1e-3 or ax >= 1e3):
        return f"{x:.6e}"
    return f"{x:.8g}"


def print_line(r):
    print(
        f"{r['tag']:8} eta0={r['eta0']:+.3g} eta1={r['eta1']:+.3g} "
        f"Lam={r['Lam']:.3g} rho={r['rho']:.3g} a={r['a']:.6e}"
    )
    print(
        f"         <d>={r['mean']:.8e}  formula={r['deq']:.8e}  "
        f"diff={r['err_mean']:.3e}  diff_with_a4={r['err_mean_a4']:.3e}"
    )
    print(
        f"         dOm={r['domega']:.8e}  formula={r['domega_pred']:.8e}  "
        f"diff={r['err_omega']:.3e}"
    )
    print(
        f"         u2={r['u2']:.8e}  pred={r['u2_pred']:.8e}  "
        f"identity={r['identity']:.3e}  phi_m={r['phi_m']:.6g}  phi_max={r['phi_max']:.6g}"
    )
    if "pendulum" in r:
        print(
            f"         pendulum(phi_max)={r['pendulum']:.8e}  "
            f"diff={r['err_pend']:.3e}"
        )


def scaling_study():
    """How the error of <delta> ~ a**2/2 scales with amplitude and with rho."""
    print("\n==== scaling of the mean-formula residual ====")
    eta0, eta1, Lam = -0.5, 0.5, 0.0
    amps = [0.008, 0.012, 0.016]
    rhos = [0.4, 8.0, 25.0]
    coeff = -(19.0 / 72.0) * eta1**3 / eta0**3
    print(f"19/72 coefficient of a**4: {coeff:.8e}")
    for rho in rhos:
        print(f"\n-- rho={rho}")
        for a in amps:
            m = integrate_full(eta0, eta1, Lam, rho, a, n_periods=5)
            pred = -(eta1 * m["a"] ** 2 / 2.0) / eta0
            err = m["mean"] - pred
            a4 = coeff * m["a"] ** 4
            print(
                f"  a_req={a:.4f} a={m['a']:.6e} err={err:.4e} "
                f"a4term={a4:.4e} err-a4={err-a4:.4e} "
                f"err/a**4={(err/m['a']**4):.4e} err/a**6={(err/m['a']**6):.4e} "
                f"phi_m={rho*abs(eta0)*m['a']:.4f}"
            )


def reduced_study():
    print("\n==== reduced ODE, cos phi = 1 ====")
    cases = [
        ("R1", -0.5, 0.8, 0.0, 0.02),
        ("R2", 0.4, -0.7, 0.0, 0.015),
        ("R3", -0.5, 0.5, 2e-5, 0.012),
        ("R4", 0.3, 1.1, -1e-5, 0.01),
    ]
    for tag, eta0, eta1, Lam, a in cases:
        m = integrate_reduced(eta0, eta1, Lam, a)
        deq = -(eta1 * m["a"] ** 2 / 2.0 + Lam) / eta0
        deq4 = deq - (19.0 / 72.0) * eta1**3 * m["a"] ** 4 / eta0**3
        # frequency: no phi term
        dom = (
            -(eta1 * Lam) / eta0**2
            - (5.0 / 12.0) * (eta1 * m["a"] / eta0) ** 2
        )
        # The a**4 back-reaction on frequency is not in the formula.
        # For the cross term the fixed-point shift also has O(Lam**2).
        u2 = eta1 * m["a"] ** 2 / (6.0 * eta0)
        print(
            f"{tag} a={m['a']:.6e} mean-deq={m['mean']-deq:.3e} "
            f"mean-deq4={m['mean']-deq4:.3e} "
            f"dOm-(formula)={m['Omega']-1-dom:.3e} "
            f"u2-pred={m['u2']-u2:.3e} identity={m['identity']:.3e}"
        )


def draft_spotchecks():
    """Numbers quoted by the old draft, recomputed.  Not used in the note."""
    print("\n==== old-draft spot checks ====")
    # |phi_m|=0.1, eta1=0.  Choose eta0=-0.5, rho=20, a=0.01 -> phi_m=0.1
    r = report_case("draft_phi", -0.5, 0.0, 0.0, 20.0, 0.01)
    print_line(r)
    print(f"  -phi_m**2/16 = {-r['phi_m']**2/16:.8e}")
    print(f"  series -pm2/16+pm4/3072 = {-r['phi_m']**2/16 + r['phi_m']**4/3072:.8e}")
    if "pendulum" in r:
        print(f"  exact pendulum from phi_max = {r['pendulum']:.8e}")

    # |phi_m|=0.2
    r2 = report_case("draft_phi2", -0.5, 0.0, 0.0, 40.0, 0.01)
    print_line(r2)

    # reduced oscillator quoted numbers: eta0=-0.534211, eta1=1.02, dm=1e-3
    # This is only an arithmetic check of their 5/12 evaluation, plus an integration.
    eta0, eta1, dm = -0.534211, 1.02, 1e-3
    pred_shift = -(5.0 / 12.0) * (eta1 * dm / eta0) ** 2
    pred_deq = -(eta1 * dm**2 / 2.0) / eta0
    print(f"  draft 5/12 shift {pred_shift:.8e}  (they quote -1.51902e-6)")
    print(f"  draft deq {pred_deq:.8e}  (they quote 9.54679e-7)")
    m = integrate_reduced(eta0, eta1, 0.0, dm, n_periods=8)
    print(
        f"  reduced integration dOm={m['Omega']-1:.8e} mean={m['mean']:.8e} "
        f"a={m['a']:.8e}"
    )
    # recompute formula at the MEASURED amplitude
    a = m["a"]
    print(
        f"  formula at measured a: dOm={-(5/12)*(eta1*a/eta0)**2:.8e} "
        f"deq={-(eta1*a**2/2)/eta0:.8e}"
    )


def table_cases():
    print("\n==== table cases ====")
    cases = [
        # tag, eta0, eta1, Lam, rho, a_guess
        ("A_phi", -0.50, 0.0, 0.0, 10.0, 0.010),
        ("B_phi+", 0.40, 0.0, 0.0, 10.0, 0.010),
        ("C_512", -0.50, 0.50, 0.0, 1.0, 0.010),
        ("D_512+", 0.40, -0.80, 0.0, 1.0, 0.010),
        ("E_cross", -0.50, 0.50, 2.0e-5, 1.0, 0.010),
        ("F_cross+", 0.40, -0.50, 2.0e-5, 1.0, 0.010),
        ("G_add", -0.50, 0.50, 0.0, 10.0, 0.010),
        ("H_Lam", -0.50, 0.0, 2.0e-4, 10.0, 0.010),
    ]
    rows = []
    for c in cases:
        r = report_case(*c)
        print_line(r)
        rows.append(r)
    return rows


def main():
    rows = table_cases()
    reduced_study()
    scaling_study()
    draft_spotchecks()
    print("\n==== LATEX ROWS ====")
    for r in rows:
        # columns: eta0, eta1, Lambda, rho, delta_m, mean num, mean form, dom num, dom form
        print(
            f"{r['eta0']:+.2f} & {r['eta1']:+.2f} & {r['Lam']:.1e} & {r['rho']:g} & "
            f"{r['a']:.4e} & {r['mean']:.4e} & {r['deq']:.4e} & "
            f"{r['domega']:.4e} & {r['domega_pred']:.4e} \\\\"
        )


if __name__ == "__main__":
    main()
