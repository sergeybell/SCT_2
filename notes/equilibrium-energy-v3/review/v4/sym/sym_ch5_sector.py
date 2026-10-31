"""
Sector magnetic dipole: Hamiltonian through cubic order and the
chromaticity / path-length identity, same method as sym_ch5_closure.py.

h = -(1+x/rho)*sqrt((1+delta)**2 - px**2 - py**2) + x/rho + x**2/(2*rho**2)

Uniform vertical field, perpendicular (sector) edges. No field index, no sextupole.
"""
from __future__ import annotations

import sympy as sp

rho = sp.symbols("rho", positive=True)
x, y, px, py, delta = sp.symbols("x y p_x p_y delta", real=True)
D, Dp = sp.symbols("D Dprime", real=True)
beta_x, beta_y = sp.symbols("beta_x beta_y", real=True)
alpha_x, alpha_y = sp.symbols("alpha_x alpha_y", real=True)
gamma_x, gamma_y = sp.symbols("gamma_x gamma_y", real=True)
Jx, Jy = sp.symbols("J_x J_y", positive=True)
apx, apy = sp.symbols("alpha_x_prime alpha_y_prime", real=True)


def hamiltonian_poly() -> sp.Expr:
    h = -(1 + x / rho) * sp.sqrt((1 + delta) ** 2 - px**2 - py**2) + x / rho + x**2 / (2 * rho**2)
    eps = sp.symbols("eps")
    sub = {x: eps * x, y: eps * y, px: eps * px, py: eps * py, delta: eps * delta}
    s = sp.expand(h.subs(sub).series(eps, 0, 4).removeO().subs(eps, 1))
    s = sp.expand(s.series(delta, 0, 2).removeO())
    return s


def truncate(expr: sp.Expr, vars_, deg: int) -> sp.Expr:
    expr = sp.expand(expr)
    acc = 0
    for t in sp.Add.make_args(expr):
        if t == 0:
            continue
        d = sp.total_degree(sp.Poly(sp.expand(t), *vars_))
        if d <= deg:
            acc += t
    return sp.expand(acc)


def forces(h: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    xp, yp = sp.symbols("xp yp")
    h0 = sp.expand(h.subs(delta, 0))
    x_of_p = sp.expand(sp.diff(h0, px))
    y_of_p = sp.expand(sp.diff(h0, py))

    def poisson(f: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(f, x) * sp.diff(h0, px)
            - sp.diff(f, px) * sp.diff(h0, x)
            + sp.diff(f, y) * sp.diff(h0, py)
            - sp.diff(f, py) * sp.diff(h0, y)
        )

    quad_x = sp.expand(x_of_p - px)
    quad_y = sp.expand(y_of_p - py)
    px_s = truncate(xp - quad_x.subs({px: xp, py: yp}), {x, y, xp, yp}, 2)
    py_s = truncate(yp - quad_y.subs({px: xp, py: yp}), {x, y, xp, yp}, 2)
    xpp = truncate(poisson(x_of_p).subs({px: px_s, py: py_s}), {x, y, xp, yp}, 2)
    ypp = truncate(poisson(y_of_p).subs({px: px_s, py: py_s}), {x, y, xp, yp}, 2)
    return xpp, ypp, px_s, py_s


def qcoeff(q: sp.Expr, mon: sp.Expr, gens) -> sp.Expr:
    return sp.factor(
        sp.Poly(sp.expand(q), *gens).coeff_monomial(sp.Poly(mon, *gens).monoms()[0])
    )


def main() -> None:
    h = hamiltonian_poly()
    print("=== cubic Hamiltonian at delta=0, deviation from -1 ===")
    hdev = sp.expand(h.subs(delta, 0) + 1)
    print(hdev)

    xp, yp = sp.symbols("xp yp")
    xpp, ypp, px_s, py_s = forces(h)
    print("\nx'' =", sp.factor(xpp))
    print("y'' =", sp.factor(ypp))
    print("px(xp) =", px_s)
    Kx = sp.simplify(-xpp.coeff(x).subs({x: 0, y: 0, xp: 0, yp: 0}))
    print("Kx =", Kx)
    Fx = sp.expand(xpp + Kx * x)
    print("Fx =", sp.factor(Fx))
    for mon in [x**2, y**2, xp**2, yp**2, x * xp]:
        print(f"  Fx[{mon}] =", sp.factor(qcoeff(Fx, mon, (x, y, xp, yp))))

    # averages
    Favg = (
        qcoeff(Fx, x**2, (x, y, xp, yp)) * Jx * beta_x
        + qcoeff(Fx, y**2, (x, y, xp, yp)) * Jy * beta_y
        + qcoeff(Fx, xp**2, (x, y, xp, yp)) * Jx * gamma_x
        + qcoeff(Fx, yp**2, (x, y, xp, yp)) * Jy * gamma_y
    )
    print("<Fx> =", sp.factor(sp.expand(Favg)))

    # Chromaticity around dispersed orbit
    xb, yb, pxb, pyb = sp.symbols("x_b y_b px_b py_b")
    h_shift = sp.expand(
        h.subs({x: xb + D * delta, px: pxb + Dp * delta, y: yb, py: pyb}).series(delta, 0, 2).removeO()
    )
    poly = sp.Poly(h_shift, xb, yb, pxb, pyb, delta)
    terms = []
    for mon, c in poly.terms():
        if mon[4] == 1 and sum(mon[:4]) == 2:
            terms.append(c * xb**mon[0] * yb**mon[1] * pxb**mon[2] * pyb**mon[3])
    q = sp.expand(sum(terms))
    print("\nquad*delta =", sp.factor(q))
    gens = (xb, yb, pxb, pyb)
    c_x2 = qcoeff(q, xb**2, gens)
    c_px2 = qcoeff(q, pxb**2, gens)
    c_xpx = qcoeff(q, xb * pxb, gens)
    c_y2 = qcoeff(q, yb**2, gens)
    c_py2 = qcoeff(q, pyb**2, gens)
    c_ypy = qcoeff(q, yb * pyb, gens)
    print("c_x2", c_x2)
    print("c_y2", c_y2)
    print("c_px2", c_px2)
    print("c_py2", c_py2)
    print("c_xpx", c_xpx)
    print("c_ypy", c_ypy)

    Ix = sp.expand(2 * c_x2 * beta_x + 2 * c_px2 * gamma_x - 2 * c_xpx * alpha_x)
    Iy = sp.expand(2 * c_y2 * beta_y + 2 * c_py2 * gamma_y - 2 * c_ypy * alpha_y)
    print("\nIx =", sp.factor(Ix))
    print("Iy =", sp.factor(Iy))

    # Linear dispersion source: x'' + Kx x - S delta at degree 1
    # From Hamiltonian, coeff of delta*x in the quadratic form gives B, etc.
    # Check explicit Ix at D=0.
    print("Ix(D=0) =", sp.factor(Ix.subs({D: 0, Dp: 0})))
    print("Iy(D=0) =", sp.factor(Iy.subs({D: 0, Dp: 0})))

    # Path-length Lambda density at delta=0 is x/rho + (xp**2+yp**2)/2,
    # so after average and x = x_beta + x_nl:
    # x_nl/rho + Jx gamma_x/2 + Jy gamma_y/2
    # Explicit match to -J/2 * I:
    target = sp.expand(Jx * gamma_x / 2 + Jy * gamma_y / 2)
    from_xi = sp.expand(
        -(Jx / 2) * Ix.subs({D: 0, Dp: 0}) - (Jy / 2) * Iy.subs({D: 0, Dp: 0})
    )
    print("\nexplicit target", target)
    print("from xi", sp.factor(from_xi))
    print("diff", sp.factor(sp.expand(target - from_xi)))

    Ix_feed = sp.expand(Ix - Ix.subs({D: 0, Dp: 0}))
    Iy_feed = sp.expand(Iy - Iy.subs({D: 0, Dp: 0}))
    from_feed = sp.expand(-(Jx / 2) * Ix_feed - (Jy / 2) * Iy_feed)
    green = sp.expand(Favg * D)
    mismatch = sp.expand(from_feed - green)
    print("\nfeed from xi =", sp.factor(from_feed))
    print("<F>D =", sp.factor(green))
    print("mismatch =", sp.factor(mismatch))
    # gamma = K beta - alpha', Kx=1/rho**2, Ky=0
    mismatch2 = sp.expand(
        mismatch.subs({gamma_x: beta_x / rho**2 - apx, gamma_y: 0 * beta_y - apy})
    )
    # Ky=0 so gamma_y = - alpha_y'
    print("after Twiss =", sp.factor(mismatch2))

    # ansatz Q
    a1, a2, a3, a4 = sp.symbols("a1 a2 a3 a4")
    S = 1 / rho
    Kx = 1 / rho**2
    Dpp = S - Kx * D
    Q = (
        a1 * Jx * D * alpha_x / rho
        + a2 * Jy * D * alpha_y / rho
        + a3 * Jx * Dp * beta_x / rho
        + a4 * Jy * Dp * beta_y / rho
    )
    dQ = (
        a1 * Jx * (Dp * alpha_x + D * apx) / rho
        + a2 * Jy * (Dp * alpha_y + D * apy) / rho
        + a3 * Jx * (Dpp * beta_x + Dp * (-2 * alpha_x)) / rho
        + a4 * Jy * (Dpp * beta_y + Dp * (-2 * alpha_y)) / rho
    )
    resid = sp.expand(mismatch2 - dQ)
    poly = sp.Poly(sp.expand(resid), Jx, Jy, D, Dp, alpha_x, alpha_y, apx, apy, beta_x, beta_y)
    eqs = [sp.simplify(c) for _, c in poly.terms() if sp.simplify(c) != 0]
    print("n eqs", len(eqs))
    sol = sp.solve(eqs, [a1, a2, a3, a4], dict=True)
    print("sol", sol)

    # Also print the dispersion linear equation coefficient of delta from x''.
    # Rebuild x'' keeping delta*linear.
    print("\n=== linear chromatic Hill derivative ===")
    # At x=0 the Hill K from x'' + Kx/(1+delta) approximately.
    # coeff: expand px' and x'.
    h1 = sp.expand(h.series(delta, 0, 2).removeO())
    # Use only linear transverse, first order delta: already in h.
    # x' = dh/dpx, at px small, x=0: 
    xp_lin = sp.series(sp.diff(h, px).subs({y: 0, py: 0, x: 0}), px, 0, 2).removeO()
    print("x'(x=0) =", sp.factor(sp.series(xp_lin, delta, 0, 2).removeO()))


if __name__ == "__main__":
    main()
