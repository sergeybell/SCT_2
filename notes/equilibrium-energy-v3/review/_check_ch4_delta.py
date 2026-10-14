# Fixed-rational checks for chapter 4. One symbol (delta); gamma and n are rationals.
import sympy as sp

delta = sp.symbols("delta")
g = sp.symbols("g", positive=True)  # gamma_s**2, used only in polynomial identity


def vsv_formula(delta, psi, g):
    # g = gamma^2
    return (
        1
        - (delta - psi) / g
        + sp.Rational(3, 2) * (delta - psi) ** 2 / g
        - delta ** 2 / (2 * g ** 2)
    )


def closed_quad(d1, d2, p2, n, g):
    eta0 = (1 + 1 / g) * d1 - 1 / g
    eta1 = (
        p2 / 2
        + (1 + 1 / g) * d2
        + (5 - n) / (2 * g) * d1 ** 2
        - 4 * d1 / g
        + (3 * g - 1) / (2 * g ** 2)
    )
    return sp.expand(eta0), sp.expand(eta1)


def product_coeffs(d1, d2, p2, n, g):
    psi = d1 * delta + (d2 - n * d1 ** 2 / 2) * delta ** 2
    vsv = vsv_formula(delta, psi, g)
    path = 1 + d1 * delta + (d2 + p2 / 2) * delta ** 2
    prod = sp.series(sp.expand(path * vsv), delta, 0, 3).removeO()
    # constant, linear, quadratic
    c0 = sp.expand(prod.subs(delta, 0))
    c1 = sp.expand(prod.coeff(delta))
    c2 = sp.expand(prod.coeff(delta ** 2))
    return c0, c1, c2


def check_algebra():
    d1, d2, p2, n = sp.symbols("d1 d2 p2 n")
    c0, c1, c2 = product_coeffs(d1, d2, p2, n, g)
    e0, e1 = closed_quad(d1, d2, p2, n, g)
    print("algebra c0", sp.simplify(c0 - 1))
    print("algebra c1", sp.simplify(c1 - e0))
    print("algebra c2", sp.simplify(c2 - e1))


def check_vsv(gamma_s):
    psi = sp.symbols("psi")
    gs2 = sp.Integer(gamma_s.p) ** 2 / sp.Integer(gamma_s.q) ** 2 if False else gamma_s ** 2
    beta2 = 1 - 1 / gs2
    gamma_out = sp.sqrt(1 + (gs2 - 1) * (1 + delta) ** 2)
    gamma_loc = gamma_out - beta2 * gamma_s * psi
    beta_s = sp.sqrt(beta2)
    beta = sp.sqrt(1 - 1 / gamma_loc ** 2)
    exact = beta_s / beta
    # joint series to total order 2
    s = exact.series(delta, 0, 3).removeO()
    s = sp.series(s, psi, 0, 3).removeO()
    # drop degree >= 3
    s = sum(
        term
        for term in sp.Add.make_args(sp.expand(s))
        if sp.total_degree(term, delta, psi) <= 2
    )
    formula = vsv_formula(delta, psi, gs2)
    formula = sp.expand(formula)
    # formula may contain degree 3 from (delta-psi)**2 * nothing; it is degree 2
    diff = sp.expand(s - formula)
    diff = sum(
        term
        for term in sp.Add.make_args(diff)
        if sp.total_degree(term, delta, psi) <= 2
    )
    print("vsv gamma", gamma_s, "diff", sp.simplify(diff))
    # psi=0 vs chapter-1 coefficient
    ch1 = (
        1
        - delta / gs2
        + (3 * gs2 - 1) / (2 * gs2 ** 2) * delta ** 2
    )
    d0 = sp.series(sp.expand(formula.subs(psi, 0) - ch1), delta, 0, 3).removeO()
    print("vsv psi0", sp.simplify(d0))


def orbit_etas(n, gamma_s):
    gs2 = gamma_s ** 2
    beta2 = sp.simplify(1 - 1 / gs2)
    k = sp.simplify(beta2 * gamma_s)  # p v /(mc^2) on reference, rho=1
    a, b = sp.symbols("a b")
    r = 1 + a * delta + b * delta ** 2
    gamma_out = sp.sqrt(1 + (gs2 - 1) * (1 + delta) ** 2)
    if n == 1:
        psi = sp.log(r)
    else:
        psi = (r ** (1 - n) - 1) / (1 - n)
    gamma_loc = gamma_out - k * psi
    lhs = gamma_loc - 1 / gamma_loc
    rhs = k * r ** (1 - n)
    res = sp.series(lhs - rhs, delta, 0, 3).removeO()
    # constant term must vanish identically
    c0 = sp.simplify(res.subs(delta, 0))
    eq1 = sp.simplify(res.coeff(delta))
    eq2 = sp.simplify(res.coeff(delta ** 2))
    sol = sp.solve([eq1, eq2], [a, b], dict=True)
    if len(sol) != 1:
        # filter the physical branch a finite
        sol = [s for s in sol if s[a].is_finite]
    sol = sol[0]
    a_v = sp.simplify(sol[a])
    b_v = sp.simplify(sol[b])
    r_s = 1 + a_v * delta + b_v * delta ** 2
    if n == 1:
        psi_s = sp.log(r_s)
    else:
        psi_s = (r_s ** (1 - n) - 1) / (1 - n)
    g_out = sp.sqrt(1 + (gs2 - 1) * (1 + delta) ** 2)
    g_loc = g_out - k * psi_s
    beta = sp.sqrt(1 - 1 / g_loc ** 2)
    vs = sp.sqrt(beta2)
    tau = sp.series(r_s * vs / beta - 1, delta, 0, 3).removeO()
    eta0_ex = sp.simplify(sp.expand(tau.coeff(delta)))
    eta1_ex = sp.simplify(sp.expand(tau.coeff(delta ** 2)))
    e0, e1 = closed_quad(a_v, b_v, 0, n, gs2)
    print(
        "n",
        n,
        "gamma",
        gamma_s,
        "res0",
        c0,
        "D/rho",
        a_v,
        "D2/rho",
        b_v,
    )
    print("  exact eta", eta0_ex, eta1_ex)
    print("  formula  ", sp.simplify(e0), sp.simplify(e1))
    print(
        "  match",
        sp.simplify(eta0_ex - e0) == 0,
        sp.simplify(eta1_ex - e1) == 0,
    )


def laplace_and_known_fields():
    x, y, rho, n = sp.symbols("x y rho n", positive=True)
    psi = (
        x / rho
        - n / 2 * (x / rho) ** 2
        + (n - 1) / 2 * (y / rho) ** 2
    )
    # Laplace: dxx + dyy + 1/(rho+x) * dx  == 0 at x=y=0
    dx = sp.diff(psi, x)
    lap = sp.diff(psi, x, 2) + sp.diff(psi, y, 2) + dx / (rho + x)
    print("laplace at 0", sp.simplify(lap.subs({x: 0, y: 0})))
    # spherical n=2 expansion
    xi, zeta = sp.symbols("xi zeta")
    R = sp.sqrt((1 + xi) ** 2 + zeta ** 2)
    psi_sph = sp.series(1 - 1 / R, xi, 0, 3).removeO()
    psi_sph = sp.series(psi_sph, zeta, 0, 3).removeO()
    psi_sph = sum(
        t
        for t in sp.Add.make_args(sp.expand(psi_sph))
        if sp.total_degree(t, xi, zeta) <= 2
    )
    expect = xi - xi ** 2 + sp.Rational(1, 2) * zeta ** 2
    print("spherical", sp.expand(psi_sph - expect))
    # cylindrical n=1
    psi_cyl = sp.series(sp.log(1 + xi), xi, 0, 3).removeO()
    print("cylindrical", sp.expand(psi_cyl - (xi - xi ** 2 / 2)))


def rational_points():
    # integrand identity at several rational points (d1, d2, p2, n, gamma)
    pts = [
        (sp.Rational(1, 1), sp.Rational(41, 50), 0, 1, sp.Rational(5, 4)),
        (sp.Rational(41, 16), sp.Rational(3, 7), 0, 2, sp.Rational(5, 4)),
        (sp.Rational(2, 3), sp.Rational(-1, 5), sp.Rational(4, 9), sp.Rational(3, 2), sp.Rational(5, 4)),
        (sp.Rational(-1, 4), sp.Rational(2, 1), sp.Rational(1, 3), 2, 2),
    ]
    for d1, d2, p2, n, gamma in pts:
        gs2 = gamma ** 2
        c0, c1, c2 = product_coeffs(d1, d2, p2, n, gs2)
        e0, e1 = closed_quad(d1, d2, p2, n, gs2)
        ok = (
            sp.simplify(c0 - 1) == 0
            and sp.simplify(c1 - e0) == 0
            and sp.simplify(c2 - e1) == 0
        )
        print("point", d1, d2, p2, n, gamma, ok, sp.simplify(e0), sp.simplify(e1))


if __name__ == "__main__":
    laplace_and_known_fields()
    check_algebra()
    check_vsv(sp.Rational(5, 4))
    check_vsv(sp.Integer(2))
    for n in (1, 2, sp.Rational(3, 2)):
        orbit_etas(n, sp.Rational(5, 4))
    orbit_etas(1, sp.Integer(2))
    orbit_etas(2, sp.Integer(2))
    rational_points()
