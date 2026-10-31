"""
Electrostatic deflector: Hamiltonian expansion through cubic order.

Independent variable s, metric h_s = 1 + x/rho.
Scaled Hamiltonian h = K/p_s = -(1+x/rho)*sqrt(Pi**2 - px**2 - py**2),
Pi**2 = p_kin**2 / p_s**2 from total energy conservation.

Orders kept:
  transverse degree <= 3 at delta = 0,
  transverse degree <= 2 times delta^1.
No numerical lattice, no floats.
"""
from __future__ import annotations

import sympy as sp

n, rho, gamma_s = sp.symbols("n rho gamma_s", positive=True)
delta = sp.symbols("delta", real=True)
x, y, px, py = sp.symbols("x y p_x p_y", real=True)
eps = sp.symbols("eps", real=True)

ig2 = 1 / gamma_s**2
beta2 = 1 - ig2
xi = x / rho
eta = y / rho

# Cubic potential of a cylindrical field E_r ~ r**(-n), off-plane from Laplace.
# a3 fixed by the median-plane power law; c3 by Laplace at order xi.
a3 = n * (n + 1) / 6
c3 = sp.symbols("c3")
psi_cubic = (
    xi
    - n * xi**2 / 2
    + (n - 1) * eta**2 / 2
    + a3 * xi**3
    + c3 * xi * eta**2
)


def laplace_residual(psi_expr: sp.Expr) -> sp.Expr:
    """psi_xixi + psi_etaeta + psi_xi/(1+xi) in dimensionless coordinates."""
    X, Y = sp.symbols("X Y")
    f = psi_expr.subs({xi: X, eta: Y})
    lap = sp.diff(f, X, 2) + sp.diff(f, Y, 2) + sp.diff(f, X) / (1 + X)
    return sp.series(sp.series(lap, X, 0, 2).removeO(), Y, 0, 2).removeO()


def fix_c3() -> sp.Expr:
    lap = laplace_residual(psi_cubic)
    # Coefficient of X (order 1). Y**2 residual starts at quartic potential.
    coef = sp.expand(lap).coeff(sp.symbols("X"))
    sol = sp.solve(coef, c3)
    if len(sol) != 1:
        raise RuntimeError(f"c3 solution not unique: {sol}, residual={lap}")
    return sp.simplify(sol[0])


def kinetic_pi2(psi: sp.Expr) -> sp.Expr:
    """Pi**2 through O(delta) and O(psi**2), dropping delta**2."""
    # E_out/E_s = 1 + beta2*delta + O(delta**2)
    # Pi**2 = (1+delta)**2 - 2*psi*(E_out/E_s) + beta2*psi**2
    return sp.expand(
        (1 + 2 * delta) - 2 * psi * (1 + beta2 * delta) + beta2 * psi**2
    )


def scaled_hamiltonian(psi: sp.Expr) -> sp.Expr:
    pi2 = kinetic_pi2(psi)
    return -(1 + xi) * sp.sqrt(pi2 - px**2 - py**2)


def series_h(hexpr: sp.Expr, n_trans: int, n_delta: int) -> sp.Expr:
    """Joint series: transverse variables ~ eps, delta ~ eps (then drop excess)."""
    sub = {
        x: eps * x,
        y: eps * y,
        px: eps * px,
        py: eps * py,
        delta: eps * delta,
    }
    s = hexpr.subs(sub)
    s = s.series(eps, 0, n_trans + 1).removeO()
    s = sp.expand(s.subs(eps, 1))
    # Drop delta**2 and higher, and delta * (transverse degree >= 3).
    s = s.series(delta, 0, n_delta + 1).removeO()
    s = sp.expand(s)
    # Remove pure numbers' higher delta already handled. Filter monomials.
    terms = sp.Add.make_args(s)
    kept = []
    vars_t = (x, y, px, py)
    for t in terms:
        td = sp.degree(sp.together(t), delta) if t.has(delta) else 0
        # total transverse degree
        mon = sp.Poly(sp.expand(t.subs(delta, 1)), *vars_t)
        # t may be independent of transverse vars
        if t.subs({x: 0, y: 0, px: 0, py: 0}) == t:
            deg_t = 0
        else:
            deg_t = sp.total_degree(sp.Poly(sp.expand(t / delta**td), *vars_t))
        if td <= n_delta and deg_t <= n_trans and not (td >= 1 and deg_t >= 3):
            kept.append(t)
    return sp.expand(sum(kept))


def monomial_coeff(expr: sp.Expr, mon: sp.Expr) -> sp.Expr:
    """Coefficient of a multivariate monomial. mon is a product of symbols."""
    e = sp.expand(expr)
    if mon == 1:
        return e.subs({x: 0, y: 0, px: 0, py: 0, delta: 0})
    # Use Poly
    poly = sp.Poly(sp.expand(e), x, y, px, py, delta)
    # mon as expression
    mpoly = sp.Poly(sp.expand(mon), x, y, px, py, delta)
    if len(mpoly.terms()) != 1:
        raise ValueError(mon)
    ((exponents, coeff_m),) = mpoly.terms()
    return sp.simplify(poly.coeff_monomial(exponents) / coeff_m)


def main() -> None:
    c3_val = fix_c3()
    print("c3 =", c3_val)
    psi = psi_cubic.subs(c3, c3_val)
    print("psi =", sp.factor(psi))

    h = series_h(scaled_hamiltonian(psi), n_trans=3, n_delta=1)
    # Subtract the reference value -1 so the deviation Hamiltonian is visible,
    # but keep delta-dependent constants.
    print("\n=== Hamiltonian deviation h - h(0) through the kept orders ===")
    h0 = h.subs({x: 0, y: 0, px: 0, py: 0, delta: 0})
    print("h(0) =", h0)
    dh = sp.expand(h - h0)

    monomials = [
        delta,
        x,
        delta * x,
        x**2,
        y**2,
        px**2,
        py**2,
        delta * x**2,
        delta * y**2,
        delta * px**2,
        delta * py**2,
        x**3,
        x * y**2,
        x * px**2,
        x * py**2,
        y**2 * px**0,
    ]
    print("\nCoefficients (term / monomial):")
    for mon in [
        delta,
        x,
        delta * x,
        x**2,
        y**2,
        px**2,
        py**2,
        delta * x**2,
        delta * y**2,
        delta * px**2,
        delta * py**2,
        x**3,
        x * y**2,
        x * px**2,
        x * py**2,
        y * px,  # should vanish
        x * px,
        delta * x * px,
    ]:
        c = monomial_coeff(dh, mon)
        print(f"  [{mon}] = {c}")

    # Linear Hill coefficients at delta=0 from the quadratic form,
    # assuming h = px**2/2 + (Kx/2) x**2 + ...
    print("\n=== Implied linear rigidities (h contains K/2 * coord**2) ===")
    cx2 = monomial_coeff(dh, x**2)
    cy2 = monomial_coeff(dh, y**2)
    cpx = monomial_coeff(dh, px**2)
    cpy = monomial_coeff(dh, py**2)
    print("Kx =", sp.simplify(2 * cx2))
    print("Ky =", sp.simplify(2 * cy2))
    print("A_px (coeff of px**2/2) =", sp.simplify(2 * cpx))
    print("A_py =", sp.simplify(2 * cpy))

    # Explicit delta derivatives of quadratic coefficients at x=0.
    print("\n=== Explicit d(coeff)/d delta at the reference orbit ===")
    print("d(coeff x**2)/ddelta =", monomial_coeff(dh, delta * x**2))
    print("d(coeff y**2)/ddelta =", monomial_coeff(dh, delta * y**2))
    print("d(coeff px**2)/ddelta =", monomial_coeff(dh, delta * px**2))
    print("d(coeff py**2)/ddelta =", monomial_coeff(dh, delta * py**2))
    print("Bdelta_x = d(Kx)/d delta from x**2 piece only =", sp.simplify(2 * monomial_coeff(dh, delta * x**2)))
    print("Bdelta_y =", sp.simplify(2 * monomial_coeff(dh, delta * y**2)))

    # Time of flight integrand v_s dt/ds = -dh/ddelta, at delta=0, degree <= 2.
    # h = series including delta * f, so -dh/ddelta |_{delta=0} = -coeff of those terms.
    tof = sp.expand((-dh).coeff(delta))
    # -dh/ddelta at delta=0 equals minus the coefficient of delta in dh,
    # but dh = ... + delta * F(transverse) + (terms without delta).
    # derivative of (terms without delta) w.r.t delta is 0.
    # Wait: dh includes delta*x etc, so diff(dh, delta) at 0 is the coeff.
    # v_s dt/ds = - diff(h, delta). h0 is constant. diff(dh, delta)|_0.
    integrand = sp.expand(sp.diff(h, delta).subs(delta, 0))
    # This is diff, so TOF factor is -integrand, and reference should be 1.
    tof_factor = sp.expand(-integrand)
    print("\n=== -dh/ddelta at delta=0 through degree 2 (equals (dl/ds)(v_s/v)) ===")
    print(tof_factor)
    for mon in [1, x, y, x**2, y**2, px**2, py**2, x * y, px * py]:
        print(f"  TOF [{mon}] = {monomial_coeff(tof_factor, mon)}")

    # Notes' integrand at delta=0, degree 2, with x' = px, y' = py:
    # 1 + x/rho + (px**2+py**2)/2 + psi/gamma**2 + (x/rho)*psi/gamma**2 + 3*psi**2/(2*gamma**2)
    psi2 = xi - n * xi**2 / 2 + (n - 1) * eta**2 / 2  # degree <= 2 is enough
    notes = (
        1
        + xi
        + (px**2 + py**2) / 2
        + psi2 / gamma_s**2
        + xi * psi2 / gamma_s**2
        + 3 * psi2**2 / (2 * gamma_s**2)
    )
    notes = series_h(notes + delta * 0, n_trans=2, n_delta=0)
    # series_h with delta*0 still expands; notes has no delta. Do direct truncation.
    notes_e = sp.expand(notes.subs({x: eps * x, y: eps * y, px: eps * px, py: eps * py}))
    notes_e = sp.expand(notes_e.series(eps, 0, 3).removeO().subs(eps, 1))
    diff_tof = sp.simplify(sp.expand(tof_factor - notes_e))
    print("\nTOF minus notes integrand (should be 0) =", diff_tof)

    # Dispersion drive: linear EOM.
    # At quadratic+linear level, h = cpx*px**2 + cx2*x**2 + cx*x + cdx*delta*x + cd0*delta
    # plus delta*px**2 etc. For the driven Hill equation at first order in delta
    # and first order in x, use the delta-dependent quadratic form carefully.
    print("\n=== Linear driven equation via Poisson bracket, O(eps) ===")
    eom = linear_eom(h)
    for key, val in eom.items():
        print(f"  {key} = {sp.factor(val)}")


def linear_eom(h: sp.Expr) -> dict[str, sp.Expr]:
    """x'' + Kx x - S delta, and y'' + Ky y, at linear order.

    Momenta are eliminated in favour of slopes. Only degree 1 is kept,
    so px = x' holds.
    """
    # Use a fresh slope expansion with eps.
    xp, yp = sp.symbols("xprime yprime", real=True)
    sub = {
        x: eps * x,
        y: eps * y,
        px: eps * px,
        py: eps * py,
        delta: eps * delta,
    }
    hs = sp.expand(h.subs(sub))
    # h already truncated; subs eps again would double-count. Use h directly
    # with its own smallness: treat all variables as O(eps) by the series
    # already stored in h. Compute EOM and truncate products.
    xprime = sp.expand(sp.diff(h, px))
    yprime = sp.expand(sp.diff(h, py))
    pxprime = sp.expand(-sp.diff(h, x))
    pyprime = sp.expand(-sp.diff(h, y))

    def deg_filter(expr: sp.Expr, maxdeg: int) -> sp.Expr:
        expr = sp.expand(expr)
        acc = 0
        for t in sp.Add.make_args(expr):
            # degree in (x,y,px,py,delta)
            if t == 0:
                continue
            poly = sp.Poly(sp.expand(t * eps), x, y, px, py, delta)
            # hack: total degree
            d = sp.total_degree(sp.Poly(sp.expand(t), x, y, px, py, delta))
            if d <= maxdeg:
                acc += t
        return sp.expand(acc)

    # x'' = dx'/ds = d(x')/dx * x' + d(x')/dpx * px' + ...
    # Keep only degree <= 1 in the final equation. Quadratic pieces of the
    # coefficients times linear velocities produce degree 2, drop them.
    xprime1 = deg_filter(xprime, 1)
    yprime1 = deg_filter(yprime, 1)
    pxprime1 = deg_filter(pxprime, 1)
    pyprime1 = deg_filter(pyprime, 1)

    # At this order x' = A px with A = 1 + O(delta), so px = x' * (1+O(delta)).
    # Substitute px -> xp, py -> yp after solving the linear relation.
    # xprime1 = c_px * px + c_d * delta * px + ... no bare x*px (degree 2).
    # Collect: xprime1 = alpha(delta) * px, with alpha = 1 - delta + O(2) at x=0,
    # but alpha may multiply only px. Check leftover.
    sol_px = sp.solve(sp.Eq(xprime1.subs(px, px), xp), px)
    # Better: series reversion. xprime1 should be linear in px.
    c_px = sp.simplify(xprime1.coeff(px))
    rest = sp.expand(xprime1 - c_px * px)
    if rest != 0:
        raise RuntimeError(f"x' not proportional to px at linear order: {rest}")
    # px = xp / c_px, c_px = 1 - delta + O(delta**2) possibly
    px_of = sp.series(xp / c_px, delta, 0, 2).removeO()
    py_c = sp.simplify(yprime1.coeff(py))
    resty = sp.expand(yprime1 - py_c * py)
    if resty != 0:
        raise RuntimeError(f"y' not proportional to py: {resty}")
    py_of = sp.series(yp / py_c, delta, 0, 2).removeO()

    # px' at linear order, then x'' = d(px_of)/ds approximately.
    # Since px = (1+delta) xp + O(2) and delta constant,
    # x''_slope related: d(px)/ds = px' = (1+delta) x''.
    # Use: x'' = d(x')/ds and x' = c_px * px, c_px depends on delta only (not s),
    # so x'' = c_px * px' .
    xpp = sp.expand(c_px * pxprime1.subs(px, px_of).subs(py, py_of))
    ypp = sp.expand(py_c * pyprime1.subs(px, px_of).subs(py, py_of))
    xpp = sp.expand(sp.series(xpp, delta, 0, 2).removeO())
    ypp = sp.expand(sp.series(ypp, delta, 0, 2).removeO())
    xpp = deg_filter(xpp.subs(xp, 0).subs(yp, 0), 1)  # slopes don't enter the linear force
    # Wait, deg_filter after subs xp=0 removes slope. Linear EOM shouldn't
    # depend on slopes. Check the part proportional to xp before dropping.
    xpp_full = sp.expand(sp.series(
        sp.expand(c_px * pxprime1.subs({px: px_of, py: py_of})),
        delta, 0, 2,
    ).removeO())
    ypp_full = sp.expand(sp.series(
        sp.expand(py_c * pyprime1.subs({px: px_of, py: py_of})),
        delta, 0, 2,
    ).removeO())

    def lin_part(expr: sp.Expr) -> sp.Expr:
        expr = sp.expand(expr)
        acc = 0
        for t in sp.Add.make_args(expr):
            d = sp.total_degree(sp.Poly(sp.expand(t), x, y, xp, yp, delta))
            if d <= 1:
                acc += t
        return sp.expand(acc)

    xpp_l = lin_part(xpp_full)
    ypp_l = lin_part(ypp_full)
    # Write x'' = -Kx x + S delta + (slope leftovers)
    return {
        "x''": xpp_l,
        "y''": ypp_l,
        "c_px": sp.series(c_px, delta, 0, 2).removeO(),
        "Kx": sp.simplify(-xpp_l.coeff(x)),
        "Sx": sp.simplify(xpp_l.coeff(delta)),
        "Ky": sp.simplify(-ypp_l.coeff(y)),
        "Sy": sp.simplify(ypp_l.coeff(delta)),
        "x''_xp_coeff": sp.simplify(xpp_l.coeff(xp)),
        "x''_y_coeff": sp.simplify(xpp_l.coeff(y)),
    }


if __name__ == "__main__":
    main()
