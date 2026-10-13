# -*- coding: utf-8 -*-
"""Symbolic check of Senichev's formula (25) against the smooth electric ring.

Units m = c = 1. gamma is the synchronous Lorentz factor.
Every printed identity is an algebraic check of a formula in
sections/06_formula25.tex.

Exit status is 0 only if every identity holds.
"""
from __future__ import annotations

import sys

import sympy as sp

gamma, delta = sp.symbols("gamma delta", positive=True)
beta2 = 1 - 1 / gamma**2


def gamma_of_delta(dd):
    """Lorentz factor at Phi = 0 for momentum p_s(1+delta)."""
    return sp.sqrt(1 + (gamma**2 - 1) * (1 + dd) ** 2)


def series_to(expr, var, order):
    return sp.series(expr, var, 0, order).removeO()


def check(name, expr):
    simp = sp.simplify(expr)
    ok = simp == 0
    print(f"[{'OK' if ok else 'FAIL'}] {name}")
    if not ok:
        print("   ", simp)
    return ok


def magnetic_speed():
    """(v-v_s)/v_s off the field, through delta**2."""
    g = gamma_of_delta(delta)
    ratio = sp.sqrt(1 - 1 / g**2) / sp.sqrt(beta2)
    ser = series_to(ratio - 1, delta, 3)
    b1 = sp.simplify(ser.coeff(delta, 1))
    b2 = sp.simplify(ser.coeff(delta, 2))
    return b1, b2


def smooth_ring(k):
    """Circular orbits in E proportional to r**(-k).

    Returns alpha0, alpha1, eta0, eta1 and the orbit speed
    coefficients b1, b2 in (v-v_s)/v_s = b1*delta + b2*delta**2.
    """
    u = sp.symbols("u")
    g1, g2, g3 = sp.symbols("g1 g2 g3")
    gser = gamma + g1 * u + g2 * u**2 + g3 * u**3
    lhs = series_to(gser - 1 / gser, u, 4)
    rhs = series_to(gamma * beta2 * (1 + u) ** (1 - k), u, 4)
    eq = sp.expand(lhs - rhs)
    sol = sp.solve([eq.coeff(u, i) for i in (1, 2, 3)], [g1, g2, g3], dict=True)[0]
    g_kin = gser.subs(sol)
    if k == 1:
        phi = gamma * beta2 * sp.log(1 + u)
    else:
        phi = gamma * beta2 * (1 - (1 + u) ** (1 - k)) / (k - 1)
    g_tot = series_to(g_kin + phi, u, 4)

    u1, u2, u3 = sp.symbols("u1 u2 u3")
    orbit = u1 * delta + u2 * delta**2 + u3 * delta**3
    mismatch = sp.expand(
        series_to(g_tot.subs(u, orbit), delta, 4) - series_to(gamma_of_delta(delta), delta, 4)
    )
    solu = sp.solve(
        [mismatch.coeff(delta, i) for i in (1, 2, 3)], [u1, u2, u3], dict=True
    )[0]
    a0 = sp.factor(sp.simplify(solu[u1]))
    a1 = sp.factor(sp.simplify(solu[u2]))

    g_of_d = series_to(g_kin.subs(u, orbit.subs(solu)), delta, 4)
    beta_k = sp.sqrt(1 - 1 / g_of_d**2)
    beta_s = sp.sqrt(beta2)
    dv = series_to(beta_k / beta_s - 1, delta, 3)
    b1 = sp.factor(sp.simplify(dv.coeff(delta, 1)))
    b2 = sp.factor(sp.simplify(dv.coeff(delta, 2)))
    # tau = (C/C_s)*(v_s/v) - 1, C/C_s = r/rho = 1+u
    tau = series_to((1 + orbit.subs(solu)) * (beta_s / beta_k) - 1, delta, 3)
    eta0 = sp.factor(sp.simplify(tau.coeff(delta, 1)))
    eta1 = sp.factor(sp.simplify(tau.coeff(delta, 2)))
    return a0, a1, eta0, eta1, b1, b2


def senichev_from_factor(a0, a1):
    """Slip coefficients obtained by writing Delta v/v_s = delta/(F_v gamma**2).

    F_v = 1/(1-alpha0), substituted before any division, so alpha0 = 1 is safe.
    cross_power = 1 is the uniform factor 1/F_v inside expansion (5).
    cross_power = 2 is the superscript printed in formula (25).
    The square (Delta v/v)^2 always carries 1/F_v**2.
    """
    # 1/F_v = 1-alpha0, 1/(F_v gamma**2) = (1-alpha0)/gamma**2
    lin = (1 - a0) / gamma**2
    square = (1 - a0) ** 2 / gamma**4
    cross_logic = a0 * (1 - a0) / gamma**2
    cross_print = a0 * (1 - a0) ** 2 / gamma**2
    eta0 = sp.factor(sp.simplify(a0 - lin))
    eta1_logic = sp.factor(sp.simplify(a1 - cross_logic + square))
    eta1_print = sp.factor(sp.simplify(a1 - cross_print + square))
    return eta0, eta1_logic, eta1_print


def expansion_five(a0, a1, b1, b2):
    """Formula (5) kept through delta**2, with a quadratic piece in Delta v."""
    # Delta C/C = a0 d + a1 d**2, Delta v/v_s = b1 d + b2 d**2
    # Delta T/T = dC - dv - dC*dv + dv**2
    d = delta
    dC = a0 * d + a1 * d**2
    dv = b1 * d + b2 * d**2
    expr = series_to(dC - dv - dC * dv + dv**2, d, 3)
    return sp.simplify(expr.coeff(d, 1)), sp.simplify(expr.coeff(d, 2))


def main():
    ok = True
    print("=== off-field speed, magnetic kinematics ===")
    b1_m, b2_m = magnetic_speed()
    print("  b1 =", b1_m)
    print("  b2 =", b2_m)
    ok &= check("b1 = 1/gamma**2", b1_m - 1 / gamma**2)
    missing_mag = sp.Rational(3, 2) * beta2 / gamma**2
    ok &= check("missing magnetic term is -b2 = 3 beta**2/(2 gamma**2)", -b2_m - missing_mag)

    # Senichev drops b2. Then eta1 = a1 - a0/gamma**2 + 1/gamma**4,
    # short of the true eta1 by -b2.
    a0, a1 = sp.symbols("alpha0 alpha1")
    e0, e1 = expansion_five(a0, a1, b1_m, 0)
    e0_full, e1_full = expansion_five(a0, a1, b1_m, b2_m)
    ok &= check("magnetic eta0 from (5)", e0 - (a0 - 1 / gamma**2))
    ok &= check(
        "bracket (5) misses 3 beta**2/(2 gamma**2)",
        (e1_full - e1) - missing_mag,
    )

    expected = {
        1: {
            "a0": 1,
            "a1": (gamma**2 + 1) / (2 * gamma**2),
            "eta0": 1,
            "eta1": (gamma**2 + 1) / (2 * gamma**2),
            "b1": 0,
            "b2": 0,
        },
        2: {
            "a0": gamma**2 + 1,
            "a1": 3 * gamma**2 + sp.Rational(3, 2) + 1 / (2 * gamma**2),
            "eta0": gamma**2 + 2,
            "eta1": (9 * gamma**4 + 7 * gamma**2 + 2) / (2 * gamma**2),
            "b1": -1,
            "b2": -(gamma**4 + 1) / (2 * gamma**2),
        },
    }

    print("=== smooth electric ring ===")
    for k, exp in expected.items():
        print(f"--- k = {k} ---")
        a0_k, a1_k, eta0_k, eta1_k, b1_k, b2_k = smooth_ring(k)
        print("  alpha0 =", a0_k)
        print("  alpha1 =", a1_k)
        print("  eta0   =", eta0_k)
        print("  eta1   =", eta1_k)
        print("  b1     =", b1_k)
        print("  b2     =", b2_k)
        ok &= check(f"k={k} alpha0", a0_k - exp["a0"])
        ok &= check(f"k={k} alpha1", sp.together(a1_k - exp["a1"]))
        ok &= check(f"k={k} eta0", eta0_k - exp["eta0"])
        ok &= check(f"k={k} eta1", sp.together(eta1_k - exp["eta1"]))
        ok &= check(f"k={k} b1", b1_k - exp["b1"])
        ok &= check(f"k={k} b2", sp.together(b2_k - exp["b2"]))
        ok &= check(
            f"k={k} eta0 = alpha0(1+1/gamma**2)-1/gamma**2",
            eta0_k - (a0_k * (1 + 1 / gamma**2) - 1 / gamma**2),
        )
        # orbit speed over off-field speed
        if b1_k == 0:
            ok &= check(f"k={k} alpha0 is 1, so F_v diverges", a0_k - 1)
        else:
            fv = sp.simplify((1 / gamma**2) / b1_k)
            print("  F_v =", fv)
            ok &= check(f"k={k} F_v = 1/(1-alpha0)", fv - 1 / (1 - a0_k))

        s0, s1_logic, s1_print = senichev_from_factor(a0_k, a1_k)
        print("  senichev eta0        =", s0)
        print("  senichev eta1 logic  =", s1_logic)
        print("  senichev eta1 printed=", s1_print)
        ok &= check(f"k={k} logical eta0 matches", s0 - eta0_k)
        # full (5) with the true orbit speed reproduces the period
        r0, r1 = expansion_five(a0_k, a1_k, b1_k, b2_k)
        ok &= check(f"k={k} (5) with full speed, eta0", sp.simplify(r0 - eta0_k))
        ok &= check(f"k={k} (5) with full speed, eta1", sp.together(r1 - eta1_k))
        # dropping b2 is the whole logical-minus-exact gap
        ok &= check(
            f"k={k} exact eta1 - logical eta1 = -b2",
            sp.together((eta1_k - s1_logic) + b2_k),
        )

    print("=== sphere: scaled magnetic term is not the gap ===")
    a0_s, a1_s, eta0_s, eta1_s, b1_s, b2_s = smooth_ring(2)
    _, s1_logic_s, s1_print_s = senichev_from_factor(a0_s, a1_s)
    gap = sp.factor(sp.together(eta1_s - s1_logic_s))
    print("  eta1 exact - logical =", gap)
    fv_s = -1 / gamma**2
    candidates = {
        "3 beta**2/(2 gamma**2)": missing_mag,
        "same over F_v": missing_mag / fv_s,
        "same over F_v**2": missing_mag / fv_s**2,
        "3 beta**2/(2 F_v**2)": sp.Rational(3, 2) * beta2 / fv_s**2,
    }
    for name, expr in candidates.items():
        same = sp.simplify(sp.together(gap - expr)).equals(0)
        # Pass only when the candidate is a different expression.
        ok &= check(f"gap is not {name}", sp.Integer(0 if not same else 1))
    ok &= check("gap equals -b2", sp.together(gap + b2_s))
    printed_same = sp.simplify(sp.together(s1_print_s - eta1_s)).equals(0)
    ok &= check(
        "printed eta1 differs from exact on the sphere",
        sp.Integer(0 if not printed_same else 1),
    )
    # closed forms quoted in the section
    ok &= check(
        "logical sphere eta1",
        sp.together(s1_logic_s - (8 * gamma**4 + 7 * gamma**2 + 1) / (2 * gamma**2)),
    )
    printed_closed = -gamma**4 + 2 * gamma**2 + sp.Rational(5, 2) + 1 / (2 * gamma**2)
    ok &= check("printed sphere eta1", sp.together(s1_print_s - printed_closed))

    print("=== transition condition ===")
    # F_v = 1/(1-alpha0) > 0, eta0 = 0 at gamma**2 = 1/alpha0 - 1,
    # which is the same as F_v = 1/(alpha0 gamma**2).
    aa = sp.symbols("alpha0", positive=True)
    fv = 1 / (1 - aa)
    eta0_fv = sp.simplify(aa - 1 / (fv * gamma**2))
    root = sp.simplify(eta0_fv.subs(gamma**2, 1 / aa - 1))
    ok &= check("eta0 vanishes at gamma**2 = 1/alpha0 - 1", root)
    boundary = sp.simplify(fv - 1 / (aa * (1 / aa - 1)))
    ok &= check("paper inequality is the same root", boundary)

    print("=== prefactor of formula (10) ===")
    # gamma**2 / (gamma**2 alpha0 - 1) = 1 / (alpha0 - 1/gamma**2) = +1/eta0
    pref = gamma**2 / (gamma**2 * aa - 1)
    eta0_mag = aa - 1 / gamma**2
    ok &= check("prefactor (10) equals +1/eta0", sp.simplify(pref - 1 / eta0_mag))

    print("=== summary ===")
    print("  cylinder: F_v -> infinity; eta0 and eta1 of (25) match the exact series")
    print("  sphere:   F_v = -1/gamma**2; eta0 matches; eta1 logical is short by -b2")
    print("            -b2 = (gamma**4+1)/(2 gamma**2)")
    print("  printed F_v**2 on the alpha0 term matches the cylinder and not the sphere")
    if ok:
        print("ALL CHECKS PASSED")
        return 0
    print("SOME CHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
