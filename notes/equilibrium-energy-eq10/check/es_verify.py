# -*- coding: utf-8 -*-
"""Symbolic checks for the electrostatic slip factor.

Units m = c = 1. Each check is an identity used in
sections/05_electrostatic.tex. Exit status is nonzero if any identity fails.
"""
from __future__ import annotations

import sympy as sp

gamma, delta, n = sp.symbols("gamma delta n", positive=True)
d1, d2, s2 = sp.symbols("d1 d2 s2")
beta2 = 1 - 1 / gamma**2

FAILURES = []


def series_to(expr, var, order):
    return sp.series(expr, var, 0, order).removeO()


def check(name, expr):
    simp = sp.simplify(sp.together(expr))
    ok = simp == 0
    print(f"[{'OK' if ok else 'FAIL'}] {name}")
    if not ok:
        print("   ", simp)
        FAILURES.append(name)
    return ok


def gamma_of_delta(dd):
    return sp.sqrt(1 + (gamma**2 - 1) * (1 + dd) ** 2)


def speed_ratio(g_loc):
    return sp.sqrt(beta2) / sp.sqrt(1 - 1 / g_loc**2)


def electric_integrand():
    """(v_s/v)(dl/ds)-1 through delta**2 on a pure-electric dispersion orbit."""
    u = d1 * delta + d2 * delta**2
    psi = series_to(u - (n / 2) * u**2, delta, 3)
    g_loc = gamma_of_delta(delta) - psi * beta2 * gamma
    ratio = series_to(speed_ratio(g_loc), delta, 3)
    dl = 1 + u + s2 * delta**2
    return series_to(ratio * dl - 1, delta, 3)


def magnetic_kappa2():
    ratio = series_to(speed_ratio(gamma_of_delta(delta)), delta, 3)
    return sp.simplify(ratio.coeff(delta, 2))


def averaged_eta1():
    """Replace local moments by ring averages.

    <d1> = alpha0, <d1**2> = m11, <d2> = md2, <s2> = alpha1 - md2.
    """
    a0, a1, m11, md2 = sp.symbols("alpha0 alpha1 m11 md2")
    loc = sp.expand(electric_integrand().coeff(delta, 2))
    loc = loc.subs(d1**2, m11)
    loc = loc.subs(d1, a0)
    loc = loc.subs(s2, a1 - md2)
    loc = loc.subs(d2, md2)
    return sp.simplify(loc), (a0, a1, m11, md2)


def main():
    potentials()
    kinematics()
    slip()
    smooth_cylinder_exact()
    smooth_sphere_orbits()
    betatron_density()
    smooth_betatron()
    wien()
    if FAILURES:
        print(f"\n{len(FAILURES)} FAILED")
        raise SystemExit(1)
    print("\nall checks passed")


def potentials():
    print("=== potential ===")
    rho, x, y = sp.symbols("rho x y", positive=True)
    # Laplace residual of the quadratic polynomial, evaluated at the origin.
    # psi = x/rho - (n/2)(x/rho)**2 + ((n-1)/2)(y/rho)**2
    u = x / rho
    w = y / rho
    psi = u - (n / 2) * u**2 + ((n - 1) / 2) * w**2
    # operators at x=y=0: d2/dx2 + (1/rho) d/dx + d2/dy2
    pxx = sp.diff(psi, x, 2)
    px = sp.diff(psi, x)
    pyy = sp.diff(psi, y, 2)
    lap = sp.simplify((pxx + px / rho + pyy).subs({x: 0, y: 0}))
    check("Laplace at the reference point", lap)

    cyl = series_to(sp.log(1 + u), u, 3)
    claim_cyl = u - u**2 / 2
    check("cylinder Taylor is n=1", (cyl - claim_cyl).removeO())

    eps, X, Y = sp.symbols("eps X Y")
    R = sp.sqrt((rho + eps * X) ** 2 + (eps * Y) ** 2)
    sph = series_to(1 - rho / R, eps, 3)
    claim_sph = (eps * X) / rho - (eps * X) ** 2 / rho**2 + (eps * Y) ** 2 / (2 * rho**2)
    check("sphere Taylor is n=2", sp.expand(sph - claim_sph))


def kinematics():
    print("=== kinematics ===")
    ratio = series_to(speed_ratio(gamma_of_delta(delta)), delta, 3)
    check("magnetic v_s/v linear", ratio.coeff(delta, 1) + 1 / gamma**2)
    kappa2 = sp.simplify(ratio.coeff(delta, 2))
    claim = 1 / gamma**4 + sp.Rational(3, 2) * beta2 / gamma**2
    check("magnetic kappa2", kappa2 - claim)
    relE = series_to((gamma_of_delta(delta) - gamma) / (beta2 * gamma), delta, 3)
    check("(E-Es)/(pv) linear", relE.coeff(delta, 1) - 1)
    check("(E-Es)/(pv) quadratic", relE.coeff(delta, 2) - 1 / (2 * gamma**2))

    # bivariate derivatives of v_s/v at delta=psi=0
    psi = sp.symbols("psi")
    g_loc = gamma_of_delta(delta) - psi * beta2 * gamma
    R = speed_ratio(g_loc)
    specs = {
        (1, 0): -1 / gamma**2,
        (0, 1): 1 / gamma**2,
        (2, 0): (3 * gamma**2 - 1) / gamma**4,
        (1, 1): -3 / gamma**2,
        (0, 2): 3 / gamma**2,
    }
    for (kd, kp), claim_d in specs.items():
        e = R
        if kd:
            e = sp.diff(e, delta, kd)
        if kp:
            e = sp.diff(e, psi, kp)
        got = sp.simplify(e.subs({delta: 0, psi: 0}))
        check(f"speed derivative d{kd} p{kp}", got - claim_d)

    g = sp.symbols("g", positive=True)
    factored = sp.factor(sp.together((g - 1 / g) - (gamma - 1 / gamma)))
    check(
        "equal pv fixes gamma",
        factored - (g - gamma) * (g * gamma + 1) / (g * gamma),
    )


def slip():
    print("=== slip factor ===")
    ing = electric_integrand()
    eta0_loc = sp.simplify(ing.coeff(delta, 1))
    check("local eta0", eta0_loc - ((2 - beta2) * d1 - 1 / gamma**2))
    eta1, (a0, a1, m11, md2) = averaged_eta1()
    kappa2 = magnetic_kappa2()
    eta1_mag = sp.simplify(a1 - a0 / gamma**2 + kappa2)
    extra = sp.simplify(eta1 - eta1_mag)
    claim = (
        md2 + sp.Rational(1, 2) * (5 - n) * m11 - 3 * a0
    ) / gamma**2
    check("eta1 extra, general n", extra - claim)
    draft = (md2 + 2 * m11 - 3 * a0) / gamma**2
    check("draft cylindrical eta1", extra.subs(n, 1) - draft)

    # pure electric eta0 and transition
    eta0 = (2 - beta2) * a0 - 1 / gamma**2
    check(
        "pure-electric eta0 rewritten",
        sp.simplify(eta0 - (a0 - (1 - a0) / gamma**2)),
    )
    # gamma**2 = 1/alpha0 - 1 solves eta0=0
    g2 = (1 - a0) / a0
    check(
        "electrostatic transition",
        sp.simplify(eta0.subs(gamma**2, g2)),
    )


def smooth_cylinder_exact():
    print("=== smooth cylinder ===")
    arg = delta + delta**2 / (2 * gamma**2)
    tau = series_to(sp.exp(arg) - 1, delta, 3)
    check("cylinder eta0", tau.coeff(delta, 1) - 1)
    eta1_c = sp.Rational(1, 2) + 1 / (2 * gamma**2)
    check("cylinder eta1", tau.coeff(delta, 2) - eta1_c)
    eta1, (a0, a1, m11, md2) = averaged_eta1()
    sub = {
        a0: 1,
        m11: 1,
        md2: eta1_c,
        a1: eta1_c,
        n: 1,
    }
    check("formula on smooth cylinder", eta1.subs(sub) - eta1_c)
    check(
        "eta0 formula on smooth cylinder",
        ((2 - beta2) * 1 - 1 / gamma**2) - 1,
    )


def smooth_sphere_orbits():
    print("=== smooth sphere ===")
    b2 = beta2
    g1, g2 = sp.symbols("g1 g2")
    uu = sp.symbols("uu")
    gser = gamma + g1 * uu + g2 * uu**2
    lhs = series_to(gser - 1 / gser, uu, 3)
    rhs = series_to(gamma * b2 / (1 + uu), uu, 3)
    eq = sp.expand(lhs - rhs)
    sol = sp.solve([eq.coeff(uu, 1), eq.coeff(uu, 2)], [g1, g2], dict=True)[0]
    g_u = gser.subs(sol)
    gE = series_to(g_u + gamma * b2 * (1 - 1 / (1 + uu)), uu, 3)
    u1, u2 = sp.symbols("u1 u2")
    user = u1 * delta + u2 * delta**2
    gE_d = series_to(gE.subs(uu, user), delta, 3)
    target = series_to(gamma_of_delta(delta), delta, 3)
    eqs = sp.expand(gE_d - target)
    solu = sp.solve(
        [eqs.coeff(delta, 1), eqs.coeff(delta, 2)], [u1, u2], dict=True
    )[0]
    g_of_d = series_to(g_u.subs(uu, user.subs(solu)), delta, 3)
    ratio = series_to(speed_ratio(g_of_d), delta, 3)
    tau = series_to((1 + user.subs(solu)) * ratio - 1, delta, 3)
    eta0_s = sp.simplify(tau.coeff(delta, 1))
    eta1_s = sp.simplify(tau.coeff(delta, 2))
    check("sphere eta0", eta0_s - (gamma**2 + 2))
    check(
        "sphere eta1",
        eta1_s - (sp.Rational(9, 2) * gamma**2 + sp.Rational(7, 2) + 1 / gamma**2),
    )
    check("sphere dispersion", sp.factor(solu[u1] - (gamma**2 + 1)))
    eta0_f = sp.simplify((2 - b2) * solu[u1] - 1 / gamma**2)
    check("sphere eta0 vs general formula", eta0_s - eta0_f)
    eta1, (a0, a1, m11, md2) = averaged_eta1()
    sub = {
        a0: solu[u1],
        a1: solu[u2],
        m11: solu[u1] ** 2,
        md2: solu[u2],
        n: 2,
    }
    check("sphere eta1 vs general formula", sp.simplify(eta1.subs(sub) - eta1_s))
    # nonrelativistic limit, independent of the series machinery above
    check("sphere eta0 at gamma=1", sp.simplify(eta0_s.subs(gamma, 1) - 3))
    check("sphere eta1 at gamma=1", sp.simplify(eta1_s.subs(gamma, 1) - 9))


def betatron_density():
    print("=== betatron density ===")
    eps, xb, yb, xnl, ang = sp.symbols("eps xb yb xnl ang")
    xb2, yb2 = sp.symbols("xb2 yb2")
    uu = eps * xb + eps**2 * xnl
    ww = eps * yb
    psi = series_to(
        uu - (n / 2) * uu**2 + ((n - 1) / 2) * ww**2, eps, 3
    )
    g_loc = gamma - psi * beta2 * gamma
    ratio = series_to(speed_ratio(g_loc), eps, 3)
    dl = 1 + uu + eps**2 * ang
    ing = sp.expand(series_to(ratio * dl - 1, eps, 3))
    lam = sp.expand(ing.coeff(eps, 2))
    lam = lam.subs({xb**2: xb2, yb**2: yb2}).subs({xb: 0, yb: 0})
    lam = sp.together(sp.expand(lam))

    def coef_of(expr, sym):
        return sp.simplify(sp.diff(expr, sym))

    check("coef x_nl", coef_of(lam, xnl) - (2 - beta2))
    check("coef angles", coef_of(lam, ang) - 1)
    cx = coef_of(lam, xb2)
    cy = coef_of(lam, yb2)
    check("coef (x/rho)^2", cx - (5 - n) / (2 * gamma**2))
    check("coef (y/rho)^2", cy - (n - 1) / (2 * gamma**2))
    check("cylinder (x/rho)^2", cx.subs(n, 1) - 2 / gamma**2)
    check("cylinder (y/rho)^2", cy.subs(n, 1))

    # local kinetic deviation on the dispersion orbit
    u = d1 * delta + d2 * delta**2
    psi_d = series_to(u - (n / 2) * u**2, delta, 3)
    g_loc_d = gamma_of_delta(delta) - psi_d * beta2 * gamma
    dkin = series_to((g_loc_d - gamma) / (beta2 * gamma), delta, 3)
    check("delta_kin linear", sp.simplify(dkin.coeff(delta, 1) - (1 - d1)))
    d2_cyl = sp.Rational(1, 2) + 1 / (2 * gamma**2)
    check(
        "cylinder local kinetic vanishes",
        sp.simplify(dkin.subs({d1: 1, d2: d2_cyl, n: 1})),
    )


def _orbit_tau(kind):
    """Second-order horizontal tau on a smooth deflector, from 1/theta_dot.

    kind is 'cylinder' or 'sphere'. Returns (coef of <x^2> in tau, x_nl/<x^2>).
    """
    eps, A, B, C, sigma = sp.symbols("eps A B C sigma")
    psi = sp.symbols("psi")
    if kind == "cylinder":
        nu = sp.sqrt(1 + 1 / gamma**2)

        def gkin_of(rr):
            return gamma * (1 - beta2 * sp.log(rr))

        def qer_of(rr):
            return -beta2 * gamma / rr

    elif kind == "sphere":
        nu = 1 / gamma

        def gkin_of(rr):
            return gamma * (1 - beta2 + beta2 / rr)

        def qer_of(rr):
            return -beta2 * gamma / rr**2

    else:
        raise ValueError(kind)

    xphys = eps * A * sp.cos(psi) + eps**2 * (B + C * sp.cos(2 * psi))
    r = 1 + xphys
    gkin = gkin_of(r)
    pth2 = beta2 * gamma**2 + eps**2 * sigma
    pth = sp.sqrt(pth2)
    dx = sp.diff(xphys, psi) * nu
    d2x = sp.diff(xphys, psi, 2) * nu**2
    d_inv = sp.diff(1 / r**2, psi) * nu
    dpr = pth * (d2x / r**2 + dx * d_inv)
    rhs = pth / r + qer_of(r) * gkin * r**2 / pth
    eq = sp.series(dpr - rhs, eps, 0, 3).removeO()
    lin = sp.expand(eq.coeff(eps, 1))
    check(f"{kind} linear tune residual", sp.simplify(lin))
    quad = sp.expand(eq.coeff(eps, 2))
    const = sp.simplify(sp.integrate(quad, (psi, 0, 2 * sp.pi)) / (2 * sp.pi))
    c2 = sp.simplify(
        sp.integrate(quad * sp.cos(2 * psi), (psi, 0, 2 * sp.pi)) / sp.pi
    )
    r_turn = 1 + eps * A + eps**2 * (B + C)
    gkin_t = gkin_of(r_turn)
    cond = sp.series(
        gkin_t**2 - 1 - pth2 / r_turn**2, eps, 0, 3
    ).removeO()
    sol = sp.solve(
        [const, c2, sp.expand(cond).coeff(eps, 2)],
        [B, C, sigma],
        dict=True,
    )[0]
    ratio_x = sp.simplify(2 * sol[B] / A**2)  # <x> / <x^2>
    inv = gkin * r**2 / pth
    inv_s = sp.series(inv, eps, 0, 3).removeO()
    inv_avg = sp.simplify(
        sp.integrate(sp.expand(inv_s), (psi, 0, 2 * sp.pi)) / (2 * sp.pi)
    )
    tau = sp.series(inv_avg * sp.sqrt(beta2) - 1, eps, 0, 3).removeO()
    tau = sp.expand(tau.subs(sol))
    # <x^2> = eps**2 * A**2 / 2
    coef = sp.factor(sp.simplify(tau / (eps**2 * A**2 / 2)))
    return coef, ratio_x


def smooth_betatron():
    print("=== smooth betatron ===")
    coef_c, xnl_c = _orbit_tau("cylinder")
    check("cylinder horizontal Lambda", coef_c)
    claim_x = -(gamma**2 + 5) / (2 * (gamma**2 + 1))
    check("cylinder x_nl / <x^2>", sp.simplify(xnl_c - claim_x))
    # path lengthening does not vanish
    nu2 = 1 + 1 / gamma**2
    path = sp.simplify(xnl_c + sp.Rational(1, 2) * nu2)
    check(
        "cylinder path is not zero",
        sp.simplify(path * 2 * gamma**2 * (gamma**2 + 1) - (1 - 3 * gamma**2)),
    )
    # the kinematic Lambda built from x_nl, angles and (x/rho)^2 vanishes
    lam = sp.simplify((2 - beta2) * xnl_c + sp.Rational(1, 2) * nu2 + 2 / gamma**2)
    check("cylinder kinematic Lambda cancels", lam)

    coef_s, xnl_s = _orbit_tau("sphere")
    check("sphere x_nl / <x^2>", sp.simplify(xnl_s + 1))
    check("sphere horizontal Lambda", sp.simplify(coef_s + beta2))
    nu2_s = 1 / gamma**2
    # n=2: (5-n)/(2 gamma**2) = 3/(2 gamma**2)
    lam_s = sp.simplify(
        (2 - beta2) * xnl_s + sp.Rational(1, 2) * nu2_s + sp.Rational(3, 2) / gamma**2
    )
    check("sphere kinematic Lambda matches 1/theta_dot", sp.simplify(lam_s - coef_s))
    path_s = sp.simplify(xnl_s + sp.Rational(1, 2) * nu2_s)
    check(
        "sphere path differs from time",
        sp.simplify(path_s - (-1 + 1 / (2 * gamma**2))),
    )

    # tilted circular orbit in the spherical field: exact tau = 0,
    # and the second-order formula reproduces it.
    alpha, psi = sp.symbols("alpha psi", positive=True)
    X = sp.cos(psi)
    Y = sp.cos(alpha) * sp.sin(psi)
    Z = sp.sin(alpha) * sp.sin(psi)
    phi = sp.atan2(Y, X)
    r = sp.sqrt(X**2 + Y**2)
    dphi = sp.diff(phi, psi)
    xp = sp.diff(r, psi) / dphi
    yp = sp.diff(Z, psi) / dphi

    def avg(expr):
        e = sp.series(expr * dphi, alpha, 0, 4).removeO()
        return sp.simplify(sp.integrate(e, (psi, 0, 2 * sp.pi)) / (2 * sp.pi))

    ax = avg(r - 1)
    ay2 = avg(Z**2)
    ayp2 = avg(yp**2)
    axp2 = avg(xp**2)
    # n=2, rho=1; <x^2> is O(alpha**4)
    tau = sp.simplify(
        (2 - beta2) * ax + (axp2 + ayp2) / 2 + ay2 / (2 * gamma**2)
    )
    check("tilted sphere formula tau", sp.series(tau, alpha, 0, 4).removeO())
    check("tilted sphere path", avg(sp.sqrt((r) ** 2 + xp**2 + yp**2) - 1))


def wien():
    print("=== Wien / mixed linear ===")
    dE = sp.symbols("dE")
    psi = dE * delta  # straight reference, geometric curvature absent
    g_loc = gamma_of_delta(delta) - psi * beta2 * gamma
    ratio = series_to(speed_ratio(g_loc), delta, 2)
    eta0_loc = sp.simplify((ratio - 1).coeff(delta, 1))
    check("Wien density", eta0_loc - (dE - 1) / gamma**2)


if __name__ == "__main__":
    main()
