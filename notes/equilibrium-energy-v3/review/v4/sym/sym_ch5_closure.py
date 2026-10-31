"""
Term-by-term closure of Lambda^E and -2*pi*(Jx*xi_x + Jy*xi_y)/C_s
for the body of a cylindrical electrostatic deflector.

Uses the cubic Hamiltonian produced by the same expansion as
sym_ch5_hamiltonian.py. All steps are symbolic identities.
No lattice, no floats.
"""
from __future__ import annotations

import sympy as sp

n, rho, gamma_s = sp.symbols("n rho gamma_s", positive=True)
x, y, px, py, delta = sp.symbols("x y p_x p_y delta", real=True)
D, Dp, beta_x, beta_y = sp.symbols("D Dprime beta_x beta_y", real=True)
alpha_x, alpha_y = sp.symbols("alpha_x alpha_y", real=True)
gamma_x, gamma_y = sp.symbols("gamma_x gamma_y", real=True)
Jx, Jy = sp.symbols("J_x J_y", positive=True)
# Twiss derivative alpha' = d alpha/ds. K*beta = gamma + alpha'.
apx, apy = sp.symbols("alpha_x_prime alpha_y_prime", real=True)

ig2 = 1 / gamma_s**2
beta2 = 1 - ig2
xi = x / rho
eta = y / rho

a3 = n * (n + 1) / 6
c3 = (1 - n**2) / 2
psi = (
    xi
    - n * xi**2 / 2
    + (n - 1) * eta**2 / 2
    + a3 * xi**3
    + c3 * xi * eta**2
)


def hamiltonian_poly() -> sp.Expr:
    pi2 = (1 + 2 * delta) - 2 * psi * (1 + beta2 * delta) + beta2 * psi**2
    h = -(1 + xi) * sp.sqrt(pi2 - px**2 - py**2)
    eps = sp.symbols("eps")
    sub = {x: eps * x, y: eps * y, px: eps * px, py: eps * py, delta: eps * delta}
    s = h.subs(sub).series(eps, 0, 4).removeO()
    s = sp.expand(s.subs(eps, 1))
    s = sp.expand(s.series(delta, 0, 2).removeO())
    return s


def coeff(expr: sp.Expr, mon: sp.Expr) -> sp.Expr:
    poly = sp.Poly(sp.expand(expr), x, y, px, py, delta)
    mpoly = sp.Poly(sp.expand(mon), x, y, px, py, delta)
    ((exponents, cmon),) = mpoly.terms()
    return sp.simplify(poly.coeff_monomial(exponents) / cmon)


def quadratic_force(h: sp.Expr) -> dict[str, sp.Expr]:
    """x'' + Kx x and y'' + Ky y through quadratic order at delta=0.

    Slopes xp, yp are restored so that the source is written in x, y, xp, yp.
    """
    xp, yp = sp.symbols("xp yp")
    # Work at delta = 0, transverse degree <= 2 in the final second derivative.
    h0 = sp.expand(h.subs(delta, 0))
    # Keep degree <= 3 in h (already).
    x_of_p = sp.expand(sp.diff(h0, px))
    y_of_p = sp.expand(sp.diff(h0, py))
    px_dot = sp.expand(-sp.diff(h0, x))
    py_dot = sp.expand(-sp.diff(h0, y))

    # Series reversion px = xp + quadratic corrections.
    # x_of_p = px + quadratic form in (x,y,px,py).
    # Iterate once: px = xp - (x_of_p - px)|_{px=xp, py=yp}
    quad_x = sp.expand(x_of_p - px)
    quad_y = sp.expand(y_of_p - py)
    px_s = sp.expand(xp - quad_x.subs({px: xp, py: yp}))
    py_s = sp.expand(yp - quad_y.subs({px: xp, py: yp}))
    # Drop cubic leftovers created by the substitution into quadratic pieces.
    px_s = truncate(px_s, {x, y, xp, yp}, 2)
    py_s = truncate(py_s, {x, y, xp, yp}, 2)

    # x'' = Poisson(x_of_p, h) = dx_of_p/d(phase)
    # Evaluate with momenta, then substitute slopes, then truncate to degree 2.
    def poisson_dt(f: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(f, x) * sp.diff(h0, px)
            - sp.diff(f, px) * sp.diff(h0, x)
            + sp.diff(f, y) * sp.diff(h0, py)
            - sp.diff(f, py) * sp.diff(h0, y)
        )

    xpp = poisson_dt(x_of_p)
    ypp = poisson_dt(y_of_p)
    xpp = truncate(xpp.subs({px: px_s, py: py_s}), {x, y, xp, yp}, 2)
    ypp = truncate(ypp.subs({px: px_s, py: py_s}), {x, y, xp, yp}, 2)
    xpp = sp.expand(xpp)
    ypp = sp.expand(ypp)

    Kx = sp.simplify(-xpp.coeff(x))
    Ky = sp.simplify(-ypp.coeff(y))
    # Remove the linear piece; what remains is the quadratic source F = x'' + K x.
    Fx = sp.expand(xpp + Kx * x)
    Fy = sp.expand(ypp + Ky * y)
    # Kx as read off must not depend on coordinates.
    if any(Fx.has(s) and False for s in ()):
        pass
    return {
        "Kx": sp.factor(Kx.subs({x: 0, y: 0, xp: 0, yp: 0})),
        "Ky": sp.factor(Ky.subs({x: 0, y: 0, xp: 0, yp: 0})),
        "Fx": sp.factor(Fx),
        "Fy": sp.factor(Fy),
        "px_of_xp": px_s,
        "py_of_yp": py_s,
    }


def truncate(expr: sp.Expr, vars_: set, deg: int) -> sp.Expr:
    expr = sp.expand(expr)
    acc = 0
    for t in sp.Add.make_args(expr):
        if t == 0:
            continue
        d = sp.total_degree(sp.Poly(sp.expand(t), *vars_))
        if d <= deg:
            acc += t
    return sp.expand(acc)


def phase_average_horizontal(Fx: sp.Expr) -> sp.Expr:
    """Replace quadratic monomials by Twiss averages. Jx, Jy factors kept.

    <x**2> = Jx beta_x, <xp**2> = Jx gamma_x, <x*xp> = -Jx alpha_x,
    and the same for y. Cross-plane products of odd factors vanish.
    """
    xp, yp = sp.symbols("xp yp")
    repl = {
        x**2: Jx * beta_x,
        y**2: Jy * beta_y,
        xp**2: Jx * gamma_x,
        yp**2: Jy * gamma_y,
        x * xp: -Jx * alpha_x,
        y * yp: -Jy * alpha_y,
        x * y: 0,
        x * yp: 0,
        y * xp: 0,
        xp * yp: 0,
    }
    # Expand and replace each quadratic monomial. Linear and cubic absent.
    e = sp.expand(Fx)
    acc = 0
    monoms = {
        x**2: Jx * beta_x,
        y**2: Jy * beta_y,
        xp**2: Jx * gamma_x,
        yp**2: Jy * gamma_y,
        x * xp: -Jx * alpha_x,
        y * yp: -Jy * alpha_y,
    }
    for mon, avg in monoms.items():
        c = sp.Poly(e, x, y, xp, yp).coeff_monomial(
            sp.Poly(mon, x, y, xp, yp).monoms()[0]
        )
        acc += sp.simplify(c) * avg
    return sp.expand(acc)


def main() -> None:
    h = hamiltonian_poly()
    forces = quadratic_force(h)
    print("Kx =", forces["Kx"])
    print("Ky =", forces["Ky"])
    print("Fx =", forces["Fx"])
    print("Fy =", forces["Fy"])
    print("px(xp) =", forces["px_of_xp"])
    print("py(yp) =", forces["py_of_yp"])

    xp, yp = sp.symbols("xp yp")
    Fx = forces["Fx"]
    Fy = forces["Fy"]
    print("\nFx coefficients:")
    for mon in [x**2, y**2, xp**2, yp**2, x * xp, y * yp, x * y, xp * yp]:
        c = sp.Poly(sp.expand(Fx), x, y, xp, yp).coeff_monomial(
            sp.Poly(mon, x, y, xp, yp).monoms()[0]
        )
        print(f"  Fx[{mon}] = {sp.factor(c)}")
    print("Fy coefficients:")
    for mon in [x * y, xp * yp, x * yp, y * xp, y**2, x**2]:
        c = sp.Poly(sp.expand(Fy), x, y, xp, yp).coeff_monomial(
            sp.Poly(mon, x, y, xp, yp).monoms()[0]
        )
        print(f"  Fy[{mon}] = {sp.factor(c)}")

    Favg = phase_average_horizontal(Fx)
    print("\n<Fx> =", sp.factor(Favg))

    # --- Hamiltonian chromaticity integrand ---
    # Quadratic form of h around x = D*delta, px = Dp*delta, at order delta^1.
    xb, yb, pxb, pyb = sp.symbols("x_b y_b px_b py_b")
    h_shift = h.subs({
        x: xb + D * delta,
        px: pxb + Dp * delta,
        y: yb,
        py: pyb,
    })
    h_shift = sp.expand(h_shift.series(delta, 0, 2).removeO())
    # Keep the part quadratic in betatron variables and linear in delta.
    quad_delta = []
    poly = sp.Poly(sp.expand(h_shift), xb, yb, pxb, pyb, delta)
    for mon, c in poly.terms():
        # mon = (ex, ey, epx, epy, ed)
        if mon[4] == 1 and sum(mon[:4]) == 2:
            quad_delta.append(c * xb**mon[0] * yb**mon[1] * pxb**mon[2] * pyb**mon[3])
    q = sp.expand(sum(quad_delta))
    print("\nQuadratic*delta Hamiltonian around the dispersed orbit:")
    print(sp.factor(q))

    def qcoeff(mon: sp.Expr) -> sp.Expr:
        return sp.factor(
            sp.Poly(q, xb, yb, pxb, pyb).coeff_monomial(
                sp.Poly(mon, xb, yb, pxb, pyb).monoms()[0]
            )
        )

    print("  coeff xb**2 =", qcoeff(xb**2))
    print("  coeff yb**2 =", qcoeff(yb**2))
    print("  coeff pxb**2 =", qcoeff(pxb**2))
    print("  coeff pyb**2 =", qcoeff(pyb**2))
    print("  coeff xb*pxb =", qcoeff(xb * pxb))
    print("  coeff yb*pyb =", qcoeff(yb * pyb))

    # <d h / d delta> = q with averages. xi = 1/(2*pi*J) * oint < >
    # Build Ix such that xi_x = 1/(4*pi) * oint Ix ds, similarly Iy.
    # <xb**2>=Jx beta_x, <pxb**2>=Jx gamma_x, <xb*pxb>=-Jx alpha_x.
    c_x2 = qcoeff(xb**2)
    c_px2 = qcoeff(pxb**2)
    c_xpx = qcoeff(xb * pxb)
    c_y2 = qcoeff(yb**2)
    c_py2 = qcoeff(pyb**2)
    c_ypy = qcoeff(yb * pyb)

    # <q> = c_x2*Jx*beta + c_px2*Jx*gamma + c_xpx*(-Jx*alpha) + vertical
    # xi_x = 1/(2*pi) * oint (<q>_x / Jx) ds
    #      = 1/(4*pi) * oint Ix, so Ix = 2*(c_x2*beta + c_px2*gamma - c_xpx*alpha)*2? 
    # 1/(2*pi) * (c_x2*beta) = 1/(4*pi) * (2 c_x2 beta)
    # Ix = 2*c_x2*beta_x + 2*c_px2*gamma_x - 2*c_xpx*alpha_x
    Ix = sp.expand(
        2 * c_x2 * beta_x + 2 * c_px2 * gamma_x - 2 * c_xpx * alpha_x
    )
    Iy = sp.expand(
        2 * c_y2 * beta_y + 2 * c_py2 * gamma_y - 2 * c_ypy * alpha_y
    )
    print("\nIx (xi_x = 1/(4*pi) oint Ix) =")
    print(sp.factor(Ix))
    print("Iy =")
    print(sp.factor(Iy))

    # Separate explicit (D=Dp=0) and feed-down.
    Ix0 = sp.factor(Ix.subs({D: 0, Dp: 0}))
    Iy0 = sp.factor(Iy.subs({D: 0, Dp: 0}))
    print("\nIx explicit =", Ix0)
    print("Iy explicit =", Iy0)
    print("Ix feed =", sp.factor(Ix - Ix.subs({D: 0, Dp: 0})))
    print("Iy feed =", sp.factor(Iy - Iy.subs({D: 0, Dp: 0})))

    # Hill derivative check at x=0: K_hill = A*B, A=1-delta, B from x**2 coeff *2
    # Compare text gap.
    Kx = (2 - n + ig2) / rho**2
    Ky = (n - 1) / rho**2
    # From hamiltonian, coeff of delta*x**2 is d(B/2)/ddelta
    Bdelta_x = 2 * coeff(h - h.subs({x: 0, y: 0, px: 0, py: 0, delta: 0}), delta * x**2)
    Adelta = 2 * coeff(h - h.subs({x: 0, y: 0, px: 0, py: 0, delta: 0}), delta * px**2)
    dKhill_x = sp.simplify(Adelta * Kx + Bdelta_x)  # A0=1
    gap = sp.simplify(dKhill_x + Kx + (5 - n) / (gamma_s**2 * rho**2))
    print("\nBdelta_x =", sp.factor(Bdelta_x))
    print("Adelta =", sp.factor(Adelta))
    print("dK_hill_x/ddelta at x=0 =", sp.factor(dKhill_x))
    print("GAP (should be 0 if text is wrong) =", gap)

    Bdelta_y = 2 * coeff(h - h.subs({x: 0, y: 0, px: 0, py: 0, delta: 0}), delta * y**2)
    dKhill_y = sp.simplify(Adelta * Ky + Bdelta_y)
    print("dK_hill_y/ddelta at x=0 =", sp.factor(dKhill_y))
    print("text dKy =", sp.factor(-Ky * (1 + ig2)))
    print("match y", sp.simplify(dKhill_y + Ky * (1 + ig2)))

    # ---------- closure identity ----------
    # C*Lambda_x density integrated equals -2*pi*Jx*xi = -Jx/2 * oint Ix
    # Working integrand contribution for plane x, times C (the oint-and):
    #   S*x_nl + Jx/2 * gamma_x + Jx/2 * (5-n)*beta_x/(gamma**2 rho**2)
    #   (+ vertical)
    S = (1 + ig2) / rho
    # Claim: oint S x_nl = oint <Fx> D   + edge, and edges handled separately.
    # First check BULK identity assuming continuous x' (no edge), i.e. whether
    # <Fx>*D  equals  - (Jx/2) * (Ix - Ix_explicit_matched)
    # Explicit part of -Jx/2 * Ix should match Jx/2*gamma + Jx/2*(5-n)beta/...
    explicit_target = sp.expand(
        Jx / 2 * gamma_x + Jx / 2 * (5 - n) * beta_x / (gamma_s**2 * rho**2)
        + Jy / 2 * gamma_y + Jy / 2 * (n - 1) * beta_y / (gamma_s**2 * rho**2)
    )
    explicit_from_xi = sp.expand(-(Jx / 2) * Ix0 - (Jy / 2) * Iy0)
    # Replace gamma = K*beta - alpha' so the difference is a derivative.
    print("\n=== Explicit (no orbit shift) Lambda density vs -J/2 * I ===")
    print("target =", explicit_target)
    print("from xi =", sp.factor(explicit_from_xi))
    diff_expl = sp.expand(explicit_target - explicit_from_xi)
    print("difference =", sp.factor(diff_expl))
    # Express gamma = K beta - alpha_prime and simplify difference.
    diff2 = sp.expand(
        diff_expl.subs({
            gamma_x: Kx * beta_x - apx,
            gamma_y: Ky * beta_y - apy,
        })
    )
    print("difference after gamma = K beta - alpha' =", sp.factor(diff2))

    # Feed-down: -J/2 * I_feed should equal <Fx>*D   ? for the x_nl channel
    # plus any total derivatives.
    Ix_feed = sp.expand(Ix - Ix.subs({D: 0, Dp: 0}))
    Iy_feed = sp.expand(Iy - Iy.subs({D: 0, Dp: 0}))
    from_xi_feed = sp.expand(-(Jx / 2) * Ix_feed - (Jy / 2) * Iy_feed)
    print("\n=== Feed-down part of -J/2 I (should match S*x_nl after Green) ===")
    print(sp.factor(from_xi_feed))

    print("\n<Fx> =", sp.factor(Favg))
    # If x_nl'' + K x_nl = <Fx>, and D'' + K D = S, and no edges,
    # oint S x_nl = oint <Fx> D.
    # Then S*x_nl corresponds, after Green, to <Fx>*D.
    green_density = sp.expand(Favg * D)
    print("<Fx>*D =", sp.factor(green_density))
    mismatch = sp.expand(from_xi_feed - green_density)
    print("feed xi minus <Fx>*D =", sp.factor(mismatch))
    # Substitute gamma = K beta - alpha' and Dp = D' .
    Kx_s = (2 - n + ig2) / rho**2
    Ky_s = (n - 1) / rho**2
    mismatch2 = sp.expand(
        mismatch.subs({
            gamma_x: Kx_s * beta_x - apx,
            gamma_y: Ky_s * beta_y - apy,
        })
    )
    print("after Twiss identity =", sp.factor(mismatch2))

    # Try to write mismatch as a total derivative involving alpha, D.
    # Candidate: d/ds ( f(D, alpha) ). We look for a linear combination
    # of D'*alpha, D*alpha', alpha*D  that reproduces mismatch2.
    print("\n=== decompose mismatch into total derivatives ===")
    decompose_mismatch(mismatch2)


def decompose_mismatch(expr: sp.Expr) -> None:
    """Write expr as dQ/ds plus a remainder, with Q built from D, alpha, beta, J.

    Inside the deflector K, rho, n are constant, so
    d/ds (D * alpha) = Dp * alpha + D * alpha'
    d/ds (D**2) = 2 D Dp
    beta' = -2 alpha, so d/ds (alpha * something).
    Also d/ds (D' * beta) etc. Use the dispersion ODE to remove D''.
    """
    Kx = (2 - n + ig2) / rho**2
    Ky = (n - 1) / rho**2
    S = (1 + ig2) / rho
    # alpha' = K*beta - gamma, but gamma already eliminated, so alpha' = apx symbol.
    # beta' = -2*alpha.
    # We treat symbols (D, Dp, alpha_x, alpha_y, apx, apy, beta_x, beta_y) as fields.
    # Express a general ansatz for Q and match coefficients.
    a1, a2, a3_, a4, a5, a6 = sp.symbols("a1 a2 a3 a4 a5 a6")
    # Q = a1*Jx*D*alpha_x/rho + a2*Jx*Dp*beta_x/rho + ... let sympy solve.
    # From structure, mismatch should be linear in J and in (D, Dp).
    Q = (
        a1 * Jx * D * alpha_x / rho
        + a2 * Jx * Dp * beta_x / rho
        + a3_ * Jy * D * alpha_y / rho
        + a4 * Jy * Dp * beta_y / rho
        + a5 * Jx * D * Dp / rho
        + a6 * Jy * D * Dp / rho
    )
    # dQ/ds. Derivatives:
    # D' = Dp, Dp' = D'' = S - Kx*D, alpha' = ap, beta' = -2 alpha
    Dpp = S - Kx * D
    dQ = (
        a1 * Jx * (Dp * alpha_x + D * apx) / rho
        + a2 * Jx * (Dpp * beta_x + Dp * (-2 * alpha_x)) / rho
        + a3_ * Jy * (Dp * alpha_y + D * apy) / rho
        + a4 * Jy * (Dpp * beta_y + Dp * (-2 * alpha_y)) / rho
        + a5 * Jx * (Dp * Dp + D * Dpp) / rho
        + a6 * Jy * (Dp**2 + D * Dpp) / rho
    )
    dQ = sp.expand(dQ)
    resid = sp.expand(expr - dQ)
    # Collect monomials in the field variables and set coefficients to 0.
    fields = (D, Dp, alpha_x, alpha_y, apx, apy, beta_x, beta_y)
    # resid is polynomial in these and J, with coeffs depending on a_i, n, gamma, rho
    poly = sp.Poly(sp.expand(resid), Jx, Jy, D, Dp, alpha_x, alpha_y, apx, apy, beta_x, beta_y)
    eqs = []
    for mon, c in poly.terms():
        c_s = sp.simplify(c)
        if c_s != 0:
            eqs.append(c_s)
    print("number of coefficient equations", len(eqs))
    sol = sp.solve(eqs, [a1, a2, a3_, a4, a5, a6], dict=True)
    print("solutions", sol)
    if sol:
        print("remainder check", sp.simplify(resid.subs(sol[0])))


if __name__ == "__main__":
    main()
