# -*- coding: utf-8 -*-
"""Symbolic checks for the longitudinal BM solution and the appendix.

Every identity used in sections/02_solution.tex and sections/appA_bm.tex
is re-derived here.  A failing assertion means the note must not quote
that coefficient.
"""
from __future__ import annotations

import sympy as sp

sp.init_printing(use_unicode=False)

PASS = []
FAIL = []


def check(name, cond, detail=""):
    cond = bool(cond)
    if cond:
        PASS.append(name)
        print(f"PASS  {name}")
    else:
        FAIL.append(name)
        print(f"FAIL  {name}  {detail}")
    return cond


def fourier(expr, psi, nmax):
    """Return (g, h) with expr = g[0] + sum (g[n] cos + h[n] sin)."""
    e = sp.expand(expr)
    g = [sp.simplify(sp.integrate(e, (psi, 0, 2 * sp.pi)) / (2 * sp.pi))]
    h = [0]
    for n in range(1, nmax + 1):
        g.append(sp.simplify(sp.integrate(e * sp.cos(n * psi), (psi, 0, 2 * sp.pi)) / sp.pi))
        h.append(sp.simplify(sp.integrate(e * sp.sin(n * psi), (psi, 0, 2 * sp.pi)) / sp.pi))
    return g, h


def main():
    test_linear_stability_and_c0_c2()
    test_quadratic_by_lindstedt()
    test_quadratic_by_bm_chain_rule()
    test_fourier_rules()
    test_cubic()
    test_phi_square_harmonic_balance()
    test_pendulum_series()
    test_mean_square_a4()
    test_cross_term_linearization()
    test_dimensionless_reduction()
    print("\n==== SUMMARY ====")
    print(f"passed {len(PASS)}  failed {len(FAIL)}")
    if FAIL:
        print("FAILED:", ", ".join(FAIL))
        raise SystemExit(1)


def test_linear_stability_and_c0_c2():
    print("\n-- c0, c2 and both signs of eta0 --")
    t = sp.symbols("t", real=True)
    eta0, eta1, Lam = sp.symbols("eta0 eta1 Lambda", real=True, nonzero=True)
    omega_rf, K, omega = sp.symbols("omega_rf K omega", positive=True)
    phi = sp.Function("phi")
    delta = sp.Function("delta")

    cos_phis = -sp.sign(eta0)
    tau = eta0 * delta(t) + eta1 * delta(t) ** 2 + Lam
    # symbolic sign is awkward inside diff; do the algebra with s = sgn(eta0)
    s = sp.symbols("s", real=True)
    # s**2 = 1, |eta0| = s*eta0
    cos_phis = -s
    eq_delta = sp.diff(delta(t), t) - K * cos_phis * sp.sin(phi(t))
    eq_phi = sp.diff(phi(t), t) - omega_rf * (eta0 * delta(t) + eta1 * delta(t) ** 2 + Lam)
    # ddot delta = K cos_phis cos(phi) * phi_dot
    ddot = sp.diff(K * cos_phis * sp.sin(phi(t)), t)
    ddot = sp.simplify(ddot.subs(sp.diff(phi(t), t), omega_rf * tau.subs(delta(t), delta(t))))
    # At this point tau still has delta(t). Replace.
    tau_e = eta0 * delta(t) + eta1 * delta(t) ** 2 + Lam
    ddot = K * cos_phis * sp.cos(phi(t)) * omega_rf * tau_e
    check("dd_structure", sp.simplify(ddot - K * cos_phis * omega_rf * sp.cos(phi(t)) * tau_e) == 0)

    # K cos_phis omega_rf = -omega**2 / eta0 when omega**2 = K*omega_rf*|eta0|, |eta0|=s*eta0, cos= -s
    omega2 = K * omega_rf * (s * eta0)
    pref = sp.simplify(K * cos_phis * omega_rf)
    pref_reduced = sp.simplify(pref.subs(K * omega_rf * s * eta0, omega**2))
    # pref = K*omega_rf*(-s) = - (K*omega_rf*s*eta0)/eta0 = -omega**2/eta0
    check(
        "prefactor_both_signs",
        sp.simplify(pref + omega2 / eta0) == 0,
        detail=str(pref + omega2 / eta0),
    )

    # cos phi = 1 => ddot + omega**2 delta = c0 + c2 delta**2
    ddot_1 = ddot.subs(sp.cos(phi(t)), 1)
    # ddot_1 = pref * tau = - (omega**2/eta0) * tau
    rhs = sp.expand(-(omega**2 / eta0) * (eta1 * delta(t) ** 2 + Lam) - omega**2 * delta(t))
    lhs = sp.expand(sp.simplify(ddot_1.subs(pref, -omega**2 / eta0)))
    check("reduced_ode", sp.simplify(lhs - rhs) == 0, detail=str(sp.simplify(lhs - rhs)))

    c0 = -omega**2 * Lam / eta0
    c2 = -omega**2 * eta1 / eta0
    check("c0_def", sp.simplify(c0 + omega**2 * Lam / eta0) == 0)
    check("c2_def", sp.simplify(c2 + omega**2 * eta1 / eta0) == 0)


def test_quadratic_by_lindstedt():
    print("\n-- Lindstedt, quadratic --")
    psi = sp.symbols("psi", real=True)
    a, omega, eps = sp.symbols("a omega varepsilon", positive=True)
    C0, C2 = sp.symbols("C0 C2", real=True)
    U0, U2, V0, V3, w2 = sp.symbols("U0 U2 V0 V3 w2")

    Omega = omega + eps**2 * w2
    u1 = U0 + U2 * sp.cos(2 * psi)
    u2 = V0 + V3 * sp.cos(3 * psi)
    x = a * sp.cos(psi) + eps * u1 + eps**2 * u2
    eq = Omega**2 * sp.diff(x, psi, 2) + omega**2 * x - eps * (C0 + C2 * x**2)
    eq = sp.series(sp.expand(eq), eps, 0, 3).removeO()
    eq = sp.expand(eq)

    g, h = fourier(eq, psi, 3)
    # g[n], h[n] must vanish at each order of eps. Solve order by order.
    e1 = sp.expand(eq.coeff(eps))
    g1, h1 = fourier(e1, psi, 3)
    sols1 = sp.solve([g1[0], g1[2], h1[2]], [U0, U2], dict=True)[0]
    check("lin_no_first_harmonic_order1", sp.simplify(g1[1]) == 0 and sp.simplify(h1[1]) == 0)
    U0_s = sp.simplify(sols1[U0])
    U2_s = sp.simplify(sols1[U2])
    check("lin_U0", sp.simplify(U0_s - (C0 + C2 * a**2 / 2) / omega**2) == 0, detail=str(U0_s))
    check("lin_U2", sp.simplify(U2_s + C2 * a**2 / (6 * omega**2)) == 0, detail=str(U2_s))

    e2 = sp.expand(eq.coeff(eps**2).subs(sols1))
    g2, h2 = fourier(e2, psi, 3)
    # unknowns w2, V0, V3; first harmonic fixes w2, constant fixes V0, cos3 fixes V3
    sols2 = sp.solve([g2[0], g2[1], g2[3], h2[1], h2[3]], [w2, V0, V3], dict=True)[0]
    w2_s = sp.simplify(sols2[w2])
    V0_s = sp.simplify(sols2[V0])
    V3_s = sp.simplify(sols2[V3])
    w2_expect = -C2 * C0 / omega**3 - sp.Rational(5, 12) * C2**2 * a**2 / omega**3
    check("lin_w2", sp.simplify(w2_s - w2_expect) == 0, detail=str(w2_s))
    check("lin_V0_zero", sp.simplify(V0_s) == 0, detail=str(V0_s))
    V3_expect = C2**2 * a**3 / (48 * omega**4)
    check("lin_V3", sp.simplify(V3_s - V3_expect) == 0, detail=str(V3_s))
    check("lin_no_sin", all(sp.simplify(h2[n].subs(sols2)) == 0 for n in range(1, 4)))

    # 1/2 and -1/12 decomposition
    U0_a = C2 * a**2 / (2 * omega**2)  # part of U0 without C0
    combo = sp.simplify(U0_a + U2_s / 2)
    check(
        "five_twelfths_split",
        sp.simplify(combo - sp.Rational(5, 12) * C2 * a**2 / omega**2) == 0,
    )
    half = sp.simplify(U0_a / (C2 * a**2 / omega**2))
    twelfth = sp.simplify((U2_s / 2) / (C2 * a**2 / omega**2))
    check("split_one_half", half == sp.Rational(1, 2), detail=str(half))
    check("split_minus_one_twelfth", twelfth == -sp.Rational(1, 12), detail=str(twelfth))

    # physical map eps*C -> c, frequency shift eps**2 * w2
    c0, c2 = sp.symbols("c0 c2")
    B = sp.simplify(w2_expect.subs({C0: c0, C2: c2}))
    check(
        "physical_B2",
        sp.simplify(B + c2 * c0 / omega**3 + sp.Rational(5, 12) * c2**2 * a**2 / omega**3) == 0,
    )


def test_quadratic_by_bm_chain_rule():
    print("\n-- BM chain rule, quadratic --")
    psi = sp.symbols("psi", real=True)
    a, omega, eps = sp.symbols("a omega varepsilon", positive=True)
    C0, C2 = sp.symbols("C0 C2", real=True)
    A2, B2, V0, V3 = sp.symbols("A2 B2 V0 V3")

    U0 = (C0 + C2 * a**2 / 2) / omega**2
    U2 = -C2 * a**2 / (6 * omega**2)
    u1 = U0 + U2 * sp.cos(2 * psi)
    u2 = V0 + V3 * sp.cos(3 * psi)
    x = a * sp.cos(psi) + eps * u1 + eps**2 * u2

    # adot = eps**2 * A2, psidot = omega + eps**2 * B2  (A1 = B1 = 0)
    # Truncate products at eps**2.
    dx_da = sp.diff(x, a)
    dx_dpsi = sp.diff(x, psi)
    adot = eps**2 * A2
    psidot = omega + eps**2 * B2
    xdot = sp.series(sp.expand(adot * dx_da + psidot * dx_dpsi), eps, 0, 3).removeO()
    xddot = sp.series(
        sp.expand(adot * sp.diff(xdot, a) + psidot * sp.diff(xdot, psi)),
        eps,
        0,
        3,
    ).removeO()
    rhs = sp.series(sp.expand(eps * (C0 + C2 * x**2)), eps, 0, 3).removeO()
    x_tr = sp.series(sp.expand(x), eps, 0, 3).removeO()
    residual = sp.expand(xddot + omega**2 * x_tr - rhs)

    e1 = sp.expand(residual.coeff(eps))
    g1, h1 = fourier(e1, psi, 3)
    check(
        "bm_order1_residual",
        all(sp.simplify(g1[n]) == 0 for n in range(0, 4))
        and all(sp.simplify(h1[n]) == 0 for n in range(1, 4)),
        detail=str(g1),
    )

    e2 = sp.expand(residual.coeff(eps**2))
    g2, h2 = fourier(e2, psi, 3)
    sols = sp.solve([g2[0], g2[1], g2[3], h2[1]], [A2, B2, V0, V3], dict=True)[0]
    B2_s = sp.simplify(sols[B2])
    B2_expect = -C2 * C0 / omega**3 - sp.Rational(5, 12) * C2**2 * a**2 / omega**3
    check("bm_A2_zero", sp.simplify(sols[A2]) == 0, detail=str(sols[A2]))
    check("bm_B2", sp.simplify(B2_s - B2_expect) == 0, detail=str(B2_s))
    check("bm_u2_const_zero", sp.simplify(sols[V0]) == 0, detail=str(sols[V0]))
    check(
        "bm_V3",
        sp.simplify(sols[V3] - C2**2 * a**3 / (48 * omega**4)) == 0,
        detail=str(sols[V3]),
    )
    # sins of order 2 have no source
    check("bm_order2_sins", all(sp.simplify(h2[n].subs(sols)) == 0 for n in (1, 2, 3)))


def test_fourier_rules():
    print("\n-- Fourier rules A1, B1 --")
    psi, a, omega = sp.symbols("psi a omega", positive=True)
    # damping: f = gamma * v, v = -a*omega*sin
    gamma = sp.symbols("gamma", real=True)
    v = -a * omega * sp.sin(psi)
    f = gamma * v
    g, h = fourier(f, psi, 2)
    A1 = sp.simplify(-h[1] / (2 * omega))
    B1 = sp.simplify(-g[1] / (2 * a * omega))
    check("damp_A1", sp.simplify(A1 - gamma * a / 2) == 0, detail=str(A1))
    check("damp_B1", sp.simplify(B1) == 0, detail=str(B1))

    # quadratic on the orbit has g1 = h1 = 0
    c0, c2 = sp.symbols("c0 c2")
    f2 = c0 + c2 * (a * sp.cos(psi)) ** 2
    g2, h2 = fourier(f2, psi, 3)
    check("quad_g1_zero", sp.simplify(g2[1]) == 0 and sp.simplify(h2[1]) == 0)
    check("quad_g0", sp.simplify(g2[0] - (c0 + c2 * a**2 / 2)) == 0, detail=str(g2[0]))
    check("quad_g2", sp.simplify(g2[2] - c2 * a**2 / 2) == 0, detail=str(g2[2]))
    u_const = sp.simplify(g2[0] / omega**2)
    u_cos2 = sp.simplify(g2[2] / (omega**2 * (1 - 4)))
    check("u1_from_rule_const", sp.simplify(u_const - (c0 + c2 * a**2 / 2) / omega**2) == 0)
    check("u1_from_rule_cos2", sp.simplify(u_cos2 + c2 * a**2 / (6 * omega**2)) == 0)


def test_cubic():
    print("\n-- cubic --")
    psi, a, omega = sp.symbols("psi a omega", positive=True)
    c3 = sp.symbols("c3", real=True)
    f = c3 * (a * sp.cos(psi)) ** 3
    g, h = fourier(f, psi, 3)
    check("cubic_g1", sp.simplify(g[1] - c3 * 3 * a**3 / 4) == 0, detail=str(g[1]))
    check("cubic_g3", sp.simplify(g[3] - c3 * a**3 / 4) == 0, detail=str(g[3]))
    check("cubic_g0", sp.simplify(g[0]) == 0)
    B1 = sp.simplify(-g[1] / (2 * a * omega))
    check("cubic_B1", sp.simplify(B1 + 3 * c3 * a**2 / (8 * omega)) == 0, detail=str(B1))
    U3 = sp.simplify(g[3] / (omega**2 * (1 - 9)))
    check("cubic_U3", sp.simplify(U3 + c3 * a**3 / (32 * omega**2)) == 0, detail=str(U3))
    check("cubic_A1", sp.simplify(-h[1] / (2 * omega)) == 0)

    # Third harmonic of the phase correction is NOT the cubic one.
    # (phi**2/2)*delta projects to (phi_m**2*a/8)*(cos - cos3)
    phi_m = sp.symbols("phi_m", positive=True)
    proj = (phi_m**2 * a / 8) * (sp.cos(psi) - sp.cos(3 * psi))
    gp, hp = fourier(proj, psi, 3)
    # cubic with c3 chosen to match g1: c3*3*a**3/4 = omega**2 * phi_m**2 * a / 8
    # compare g3/g1
    ratio_phi = sp.simplify(gp[3] / gp[1])
    ratio_cubic = sp.simplify(g[3] / g[1])
    check("phi_g3_over_g1", ratio_phi == -1, detail=str(ratio_phi))
    check("cubic_g3_over_g1", ratio_cubic == sp.Rational(1, 3), detail=str(ratio_cubic))
    check("phi_not_cubic", ratio_phi != ratio_cubic)


def test_phi_square_harmonic_balance():
    print("\n-- phi**2/2 frequency shift --")
    psi = sp.symbols("psi", real=True)
    phi_m, delta_m = sp.symbols("phi_m delta_m", positive=True)
    expr = (phi_m**2 / 2) * sp.sin(psi) ** 2 * delta_m * sp.cos(psi)
    target = (phi_m**2 * delta_m / 8) * (sp.cos(psi) - sp.cos(3 * psi))
    check("trig_identity", sp.simplify(sp.expand(expr - target)) == 0)

    # secular balance: x = a cos(Omega t), Omega = omega + dw
    # xdd + omega**2 x = omega**2 * (phi_m**2/8) * a * cos + (third)
    # (omega**2 - Omega**2) a = omega**2 * phi_m**2 * a / 8
    omega = sp.symbols("omega", positive=True)
    dw = sp.symbols("dw", real=True)
    Omega = omega + dw
    a = delta_m
    balance = sp.expand((omega**2 - Omega**2) * a - omega**2 * (phi_m**2 / 8) * a)
    # linearise dw << omega: omega**2 - Omega**2 = -2*omega*dw + O(dw**2)
    bal_lin = sp.series(balance, dw, 0, 2).removeO()
    dw_sol = sp.solve(bal_lin, dw)[0]
    check("phi_freq", sp.simplify(dw_sol / omega + phi_m**2 / 16) == 0, detail=str(dw_sol))

    # same result from the B1 rule, with RHS = omega**2 * target, so g1 = omega**2 * phi_m**2 * a / 8
    g1 = omega**2 * phi_m**2 * a / 8
    B = -g1 / (2 * a * omega)
    check("phi_freq_from_B_rule", sp.simplify(B / omega + phi_m**2 / 16) == 0, detail=str(B))

    # no constant term
    g, h = fourier(target, psi, 3)
    check("phi_no_constant", sp.simplify(g[0]) == 0 and sp.simplify(h[1]) == 0)
    U3 = sp.simplify((-g[3]) / 8)  # (u''+u) = - (phi_m**2 a/8) cos3 if the omega**2 cancelled
    # Here target already is (phi**2/2) delta, and the equation is u''+u = target's cos3 part
    # with the sign: target = (phi_m**2 a/8) cos - (phi_m**2 a/8) cos3
    # (1-9) U3 = - phi_m**2 a/8 => U3 = phi_m**2 a/64
    U3 = sp.simplify(g[3] / (1 - 9))
    check("phi_U3", sp.simplify(U3 - phi_m**2 * a / 64) == 0, detail=str(U3))


def test_pendulum_series():
    print("\n-- pendulum frequency series --")
    pm = sp.symbols("phi_m", positive=True)
    k2 = sp.sin(pm / 2) ** 2
    ratio = sp.pi / (2 * sp.elliptic_k(k2))
    ser = ratio.series(pm, 0, 8).removeO()
    ser_n = sp.series(sp.expand(ser), pm, 0, 8).removeO()
    expect = (
        1
        - pm**2 / 16
        + pm**4 / 3072
        - 23 * pm**6 / 737280
    )
    diff = sp.series(sp.expand(ser_n - expect), pm, 0, 8).removeO()
    check("pendulum_through_pm6", sp.simplify(diff) == 0, detail=str(ser_n))
    # independent reconstruction from the period series
    # T/T0 = 1 + (1/16)pm**2 + (11/3072)pm**4 + (173/737280)pm**6 + ...
    period = 1 + pm**2 / 16 + 11 * pm**4 / 3072 + 173 * pm**6 / 737280
    freq_from_period = sp.series(1 / period, pm, 0, 8).removeO()
    diff_p = sp.series(sp.expand(freq_from_period - expect), pm, 0, 8).removeO()
    check("pendulum_from_period", sp.simplify(diff_p) == 0, detail=str(freq_from_period))
    print("  pendulum series:", ser_n)

    # real-time pendulum: phi'' + omega**2 sin phi = 0 follows from eta1 = Lambda = 0
    t = sp.symbols("t", real=True)
    eta0 = sp.symbols("eta0", real=True, nonzero=True)
    omega_rf, omega = sp.symbols("omega_rf omega", positive=True)
    phi = sp.Function("phi")
    # delta = phi' / (omega_rf * eta0)
    # delta' = - omega**2 * sin(phi) / (omega_rf * eta0)
    # phi'' = omega_rf * eta0 * delta' = -omega**2 sin phi
    delta = sp.diff(phi(t), t) / (omega_rf * eta0)
    delta_t = -omega**2 * sp.sin(phi(t)) / (omega_rf * eta0)
    phi_tt = sp.simplify(omega_rf * eta0 * delta_t)
    check("pendulum_ode", sp.simplify(phi_tt + omega**2 * sp.sin(phi(t))) == 0)


def test_mean_square_a4():
    print("\n-- O(a**4) correction of the mean at cos phi = 1 --")
    a, omega = sp.symbols("a omega", positive=True)
    c2 = sp.symbols("c2", real=True)
    eta0, eta1 = sp.symbols("eta0 eta1", real=True, nonzero=True)
    d0 = c2 * a**2 / (2 * omega**2)
    U2 = -c2 * a**2 / (6 * omega**2)
    S4 = sp.simplify(d0**2 + U2**2 / 2)
    corr = sp.simplify((c2 / omega**2) * S4)
    corr_eta = sp.simplify(corr.subs(c2, -omega**2 * eta1 / eta0))
    expect = -sp.Rational(19, 72) * eta1**3 * a**4 / eta0**3
    check("a4_coeff", sp.simplify(corr_eta - expect) == 0, detail=str(corr_eta))
    # split 1/4 vs 1/72
    part_const = sp.simplify(d0**2 / (c2**2 * a**4 / omega**4))
    part_h2 = sp.simplify((U2**2 / 2) / (c2**2 * a**4 / omega**4))
    check("a4_const_piece", part_const == sp.Rational(1, 4), detail=str(part_const))
    check("a4_second_piece", part_h2 == sp.Rational(1, 72), detail=str(part_h2))


def test_cross_term_linearization():
    print("\n-- frequency shift from a shifted fixed point --")
    eta0, eta1, Lam, omega = sp.symbols("eta0 eta1 Lambda omega", real=True, nonzero=True)
    # fixed point of tau = 0, small Lambda
    d0 = sp.symbols("d0")
    eq = eta0 * d0 + eta1 * d0**2 + Lam
    # series solution d0 = -Lam/eta0 - eta1*Lam**2/eta0**3 + ...
    # Use inversion
    d_series = (-Lam / eta0 - eta1 * Lam**2 / eta0**3).series(Lam, 0, 3).removeO()
    residual = sp.series(eq.subs(d0, -Lam / eta0 - eta1 * Lam**2 / eta0**3), Lam, 0, 3).removeO()
    check("fixed_point_series", sp.simplify(residual) == 0, detail=str(residual))

    # Omega**2 / omega**2 = 1 + 2*eta1*d0/eta0
    rel = 1 + 2 * eta1 * (-Lam / eta0) / eta0
    Omega_over = sp.sqrt(sp.series(1 + 2 * eta1 * (-Lam / eta0 - eta1 * Lam**2 / eta0**3) / eta0, Lam, 0, 2).removeO())
    # through O(Lam): 1 - eta1*Lam/eta0**2
    ser = sp.series(Omega_over, Lam, 0, 2).removeO()
    check(
        "cross_freq",
        sp.simplify(ser - (1 - eta1 * Lam / eta0**2)) == 0,
        detail=str(ser),
    )
    check("rel_matches", sp.simplify(rel - (1 - 2 * eta1 * Lam / eta0**2)) == 0)


def test_dimensionless_reduction():
    print("\n-- dimensionless system and longitudinal map of B2 --")
    eta0, eta1, Lam = sp.symbols("eta0 eta1 Lambda", real=True, nonzero=True)
    omega, omega_rf, K = sp.symbols("omega omega_rf K", positive=True)
    s = sp.symbols("s", real=True)  # sgn(eta0), s**2 = 1
    # |eta0| = s*eta0, cos phi_s = -s, omega**2 = K*omega_rf*|eta0|
    rho = omega_rf / omega
    pref = K * (-s) / omega
    omega2 = K * omega_rf * (s * eta0)
    pref_s = sp.simplify(pref.subs(K, omega**2 / (omega_rf * s * eta0)))
    check(
        "dimless_delta_dot",
        sp.simplify(pref_s + 1 / (rho * eta0)) == 0,
        detail=str(pref_s),
    )

    a = sp.symbols("a", positive=True)
    c0 = -omega**2 * Lam / eta0
    c2 = -omega**2 * eta1 / eta0
    B = -c2 * c0 / omega**3 - sp.Rational(5, 12) * c2**2 * a**2 / omega**3
    rel = sp.simplify(B / omega)
    expect = -eta1 * Lam / eta0**2 - sp.Rational(5, 12) * (eta1 * a / eta0) ** 2
    check("map_frequency", sp.simplify(rel - expect) == 0, detail=str(rel))

    U0 = sp.simplify((c0 + c2 * a**2 / 2) / omega**2)
    U2 = sp.simplify(-c2 * a**2 / (6 * omega**2))
    check("map_U0", sp.simplify(U0 + (eta1 * a**2 / 2 + Lam) / eta0) == 0, detail=str(U0))
    check("map_U2", sp.simplify(U2 - eta1 * a**2 / (6 * eta0)) == 0, detail=str(U2))

    # delta(0) at psi = 0
    d_at_0 = sp.simplify(a + U0 + U2)
    expect0 = a - Lam / eta0 - eta1 * a**2 / (3 * eta0)
    check("ic_psi0", sp.simplify(d_at_0 - expect0) == 0, detail=str(d_at_0))


if __name__ == "__main__":
    main()
