"""Analytic slip factor of an electrostatic dispersion orbit.

One relativistic variable gamma. beta**2 = 1 - 1/gamma**2.
The local curvature term is d1 = (D/rho), the second-order orbit
x/rho = d1*delta + d2*delta**2, and x' = xp*delta.
Potential of a cylindrical deflector: phi = beta**2*gamma*log(1 + x/rho)
in units of the rest energy, zero on the design arc.
A switch `electric` turns that potential off (magnetic limit).
"""

from __future__ import annotations

import math

import sympy as sp

GAMMA_NUM = 1.0 + 233.0 / 938.2720813


def series_coefficients(electric: bool):
    g = sp.symbols("gamma", positive=True)
    d1, d2, xp, delta = sp.symbols("d1 d2 xp delta")
    b2 = 1 - 1 / g**2
    xi = d1 * delta + d2 * delta**2
    gpath = 1 + xi + sp.Rational(1, 2) * (xp * delta) ** 2
    p_s = sp.sqrt(g**2 - 1)
    gamma_gap = sp.sqrt(1 + p_s**2 * (1 + delta) ** 2)
    phi = b2 * g * sp.log(1 + xi) if electric else 0
    gamma_k = gamma_gap - phi
    beta_s = sp.sqrt(b2)
    vs_over_v = beta_s * gamma_k / sp.sqrt(gamma_k**2 - 1)
    vs = sp.series(vs_over_v, delta, 0, 3).removeO()
    gp = sp.series(gpath, delta, 0, 3).removeO()
    prod = sp.series(sp.expand(vs * gp), delta, 0, 3).removeO()
    rel = sp.expand(prod - 1)
    c1 = sp.factor(sp.simplify(rel.coeff(delta, 1)))
    c2 = sp.factor(sp.simplify(rel.coeff(delta, 2)))
    # geometric second-order piece that defines alpha1 locally
    g2 = d2 + sp.Rational(1, 2) * xp**2
    c2_sub = sp.factor(sp.simplify(c2.subs(d2, g2 - sp.Rational(1, 2) * xp**2)))
    return g, d1, g2, xp, c1, c2_sub


def report_case(name, electric):
    g, d1, g2, xp, c1, c2 = series_coefficients(electric)
    print(f"== {name} ==")
    print("c1 =", c1)
    print("c2 =", c2)
    # moments: alpha0 = <d1>, alpha1 = <g2>, q = <d1**2>
    a0, a1, q = sp.symbols("alpha0 alpha1 q")
    # replace monomials by their ring averages, assuming <d1> = a0, <g2> = a1, <d1**2> = q
    # c2 may contain d1, d1**2, g2, xp**2. xp should have cancelled into g2.
    free = [s for s in c2.free_symbols if s != g]
    print("symbols in c2:", free)
    return g, d1, g2, c1, c2


def average_eta(c1, c2, g, d1, g2):
    """Replace local monomials by ring averages alpha0, alpha1, q=<d1**2>."""
    a0, a1, q = sp.symbols("alpha0 alpha1 q")
    # Expand into polynomial in d1, g2.
    poly = sp.expand(c1)
    # c1 should be linear in d1 plus a constant
    eta0 = sp.expand(poly.subs({d1: a0}))
    poly2 = sp.expand(c2)
    # Any xp left is an error in the reduction.
    expr = poly2
    expr = expr.subs(d1**2, q)
    # after that, leftover d1 -> a0, g2 -> a1, and d1*g2 should not exist at this order
    expr = sp.expand(expr)
    expr = expr.subs({d1: a0, g2: a1})
    expr = sp.factor(sp.simplify(expr))
    return sp.factor(sp.simplify(eta0)), expr


def smooth_check(g, d1, g2, c1, c2):
    """Full cylindrical ring: d1=1, alpha1 = 1/2 + 1/(2*gamma**2), no angle."""
    a1 = sp.Rational(1, 2) + 1 / (2 * g**2)
    e0 = sp.simplify(c1.subs({d1: 1}))
    e1 = sp.simplify(c2.subs({d1: 1, g2: a1}))
    print("smooth eta0 =", e0)
    print("smooth eta1 =", sp.factor(e1))
    print("smooth target eta1 =", sp.factor(a1))
    print("smooth eta1 - target =", sp.factor(sp.simplify(e1 - a1)))


def magnetic_target(g):
    a0, a1 = sp.symbols("alpha0 alpha1")
    b2 = 1 - 1 / g**2
    return sp.factor(a1 - a0 / g**2 + 1 / g**4 + sp.Rational(3, 2) * b2 / g**2)


def drift_cosy_check():
    """COSY time coordinate on a field-free drift, expanded in delta_K and in delta."""
    g, dK, delta = sp.symbols("gamma delta_K delta", positive=True)
    # local kinetic gamma from delta_K: gamma = 1 + (gamma_s-1)*(1+delta_K)
    gamma_k = 1 + (g - 1) * (1 + dK)
    p = sp.sqrt(gamma_k**2 - 1)
    p0 = sp.sqrt(g**2 - 1)
    xp5 = (gamma_k / g) * (p0 / p)  # a=0, h=0
    dl_ds = -(xp5 - 1) * g / (g + 1)
    series = sp.series(dl_ds, dK, 0, 3).removeO()
    c1 = sp.simplify(series.coeff(dK, 1))
    c2 = sp.simplify(series.coeff(dK, 2))
    print("== COSY drift, dl/ds vs delta_K ==")
    print(" (5|6)/L =", sp.factor(c1))
    print(" (5|66)/L =", sp.factor(c2))
    # convert to (T-Ts)/Ts = -(gamma+1)/gamma * l/L, with delta_K(delta)
    dK_of_d = (g + 1) / g * delta + (g + 1) / (2 * g**3) * delta**2
    l_over_L = c1 * dK_of_d + c2 * dK_of_d**2
    rel = sp.series(sp.expand(-(g + 1) / g * l_over_L), delta, 0, 3).removeO()
    b2 = 1 - 1 / g**2
    target1 = -1 / g**2
    target2 = 1 / g**4 + sp.Rational(3, 2) * b2 / g**2
    print(" eta0 =", sp.factor(rel.coeff(delta, 1)), " target", sp.factor(target1))
    print(
        " eta1 =",
        sp.factor(sp.simplify(rel.coeff(delta, 2))),
        " target",
        sp.factor(target2),
    )
    print(
        " eta1 residual",
        sp.factor(sp.simplify(rel.coeff(delta, 2) - target2)),
    )


def velocity_amplitude():
    """v_s/v - 1 at delta=0 as a series in xi = x/rho, cylindrical potential."""
    g, xi = sp.symbols("gamma xi", positive=True)
    b2 = 1 - 1 / g**2
    phi = b2 * g * sp.log(1 + xi)
    gamma_k = g - phi
    beta_s = sp.sqrt(b2)
    vs = beta_s * gamma_k / sp.sqrt(gamma_k**2 - 1)
    rel = sp.series(vs - 1, xi, 0, 3).removeO()
    print("== v_s/v - 1 at delta=0 ==")
    print(sp.factor(sp.simplify(rel)))
    print(" linear", sp.factor(rel.coeff(xi, 1)))
    print(" quadratic", sp.factor(sp.simplify(rel.coeff(xi, 2))))


def numbers():
    g = GAMMA_NUM
    b2 = 1.0 - 1.0 / g**2
    print("== numbers, proton 233 MeV ==")
    print(f"gamma={g:.10f} beta2={b2:.10f} 1/gamma2={1/g**2:.10f} 2-beta2={2-b2:.10f}")
    a0 = 0.0964
    eta_mag = a0 - 1 / g**2
    eta_es = (2 - b2) * a0 - 1 / g**2
    print(f"mh alpha0={a0}")
    print(f"eta0_mag={eta_mag:.6f} eta0_es={eta_es:.6f}")
    print(f"gamma_tr_mag={1/math.sqrt(a0):.6f} gamma_tr_es={math.sqrt(1/a0 - 1):.6f}")
    print(f"F_v=1/(1-alpha0)={1/(1-a0):.6f}")
    # model ring from the summary
    for tag, a0m, a1m, eta1 in (
        ("model S=0", 0.63745000, 0.52577761, 0.75898477),
        ("model S=0.004", 0.63745000, 0.25782781, 0.31908628),
    ):
        f = 1.0 / (1.0 - a0m)
        sen = a1m - a0m / (f * g**2) + 1.0 / (f**2 * g**4)
        sen3 = sen + 1.5 * b2 / g**2
        print(
            f"{tag}: F={f:.6f} sen25={sen:.6f} sen25+3b2/2g2={sen3:.6f} "
            f"track={eta1:.6f} mag={a1m - a0m/g**2 + 1/g**4 + 1.5*b2/g**2:.6f}"
        )


def main():
    g, d1, g2, c1, c2 = report_case("electrostatic cylinder", True)
    eta0, eta1 = average_eta(c1, c2, g, d1, g2)
    print("eta0 =", eta0)
    print("eta1 =", eta1)
    smooth_check(g, d1, g2, c1, c2)
    print("magnetic formula target", magnetic_target(g))
    # difference eta1_es - eta1_mag
    a0, a1, q = sp.symbols("alpha0 alpha1 q")
    diff = sp.factor(sp.simplify(eta1 - magnetic_target(g)))
    print("eta1_es - eta1_mag =", diff)
    print()
    g2s, d1m, g2m, c1m, c2m = report_case("magnetic, phi=0", False)
    e0m, e1m = average_eta(c1m, c2m, g2s, d1m, g2m)
    print("mag eta0 =", e0m)
    print("mag eta1 =", e1m)
    print("mag eta1 - target", sp.factor(sp.simplify(e1m - magnetic_target(g2s))))
    print()
    drift_cosy_check()
    print()
    velocity_amplitude()
    print()
    numbers()


if __name__ == "__main__":
    main()
