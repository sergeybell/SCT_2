# Symbolic check of local gamma, delta_loc and v_s/v inside the deflector.
# Joint order is tracked by a bookkeeping parameter eps. No numeric substitution.
import sympy as sp

g, delta, psi, eps = sp.symbols("g delta psi eps", real=True)
# g = gamma_s**2 > 1. beta_s**2 = 1 - 1/g. gamma_s = sqrt(g).
beta2 = 1 - 1 / g
gamma_s = sp.sqrt(g)
k = beta2 * gamma_s  # = beta_s**2 * gamma_s = p_s v_s / (m c**2)


def main():
    d = eps * delta
    p = eps * psi
    gamma_out = sp.sqrt(1 + (g - 1) * (1 + d) ** 2)
    gamma_loc = gamma_out - k * p

    # --- claim 1: gamma_loc series, quadratic piece depends on delta only ---
    gamma_series = sp.expand(sp.series(gamma_loc, eps, 0, 3).removeO().subs(eps, 1))
    claimed_gamma = (
        gamma_s
        + beta2 * gamma_s * (delta - psi)
        + beta2 * delta ** 2 / (2 * gamma_s)
    )
    diff_gamma = sp.simplify(sp.together(gamma_series - claimed_gamma))
    print("gamma_loc - claimed:", diff_gamma)

    # quadratic monomials actually present
    quad = sp.expand(gamma_series - gamma_s - beta2 * gamma_s * (delta - psi))
    print("gamma quadratic piece:", sp.factor(quad))

    # --- claim 2: delta_loc from the same gamma series ---
    # gamma(delta_loc) = gamma_s + beta2*gamma_s*delta_loc + beta2*delta_loc**2/(2*gamma_s)
    # set equal to claimed_gamma and solve order by order
    dl = sp.symbols("dl")
    gamma_of = (
        gamma_s
        + beta2 * gamma_s * dl
        + beta2 * dl ** 2 / (2 * gamma_s)
    )
    # invert the series: dl = (delta - psi) + quadratic
    A, B = sp.symbols("A B")
    # general quadratic ansatz in delta, psi
    ansatz = (delta - psi) + A * delta ** 2 + B * psi ** 2 + sp.symbols("C") * delta * psi
    C = sp.symbols("C")
    ansatz = (delta - psi) + A * delta ** 2 + B * psi ** 2 + C * delta * psi
    residual = sp.expand(gamma_of.subs(dl, ansatz) - claimed_gamma)
    residual = sum(
        term
        for term in sp.Add.make_args(residual)
        if sp.total_degree(term, delta, psi) <= 2
    )
    residual = sp.expand(residual)
    eqs = [residual.coeff(delta ** 2), residual.coeff(psi ** 2), residual.coeff(delta * psi)]
    sol = sp.solve(eqs, [A, B, C], dict=True)[0]
    print("series inversion A,B,C:", {kk: sp.simplify(vv) for kk, vv in sol.items()})
    claimed_dl = delta - psi + (2 * delta * psi - psi ** 2) / (2 * g)
    series_dl = sp.expand(ansatz.subs(sol))
    print("series delta_loc - claimed:", sp.simplify(series_dl - claimed_dl))

    # exact kinetic-momentum definition, then series
    # gamma_loc**2 = 1 + (g-1)*(1+delta_loc)**2
    gamma_loc_s = sp.series(gamma_loc, eps, 0, 3).removeO()
    ratio = (gamma_loc_s ** 2 - 1) / (g - 1)
    delta_loc_exact = sp.series(sp.sqrt(ratio) - 1, eps, 0, 3).removeO().subs(eps, 1)
    delta_loc_exact = sp.expand(delta_loc_exact)
    print("exact delta_loc - claimed:", sp.simplify(sp.together(delta_loc_exact - claimed_dl)))

    # --- claim 3: magnetic vsv at delta_loc equals eq:vsv-es ---
    def vsv_mag(x):
        return 1 - x / g + (3 * g - 1) / (2 * g ** 2) * x ** 2

    vsv_at = sp.expand(vsv_mag(claimed_dl))
    # drop degree >= 3
    vsv_at = sum(
        term
        for term in sp.Add.make_args(sp.expand(vsv_at))
        if sp.total_degree(term, delta, psi) <= 2
    )
    vsv_es = (
        1
        - (delta - psi) / g
        + sp.Rational(3, 2) * (delta - psi) ** 2 / g
        - delta ** 2 / (2 * g ** 2)
    )
    print("vsv(delta_loc) - vsv_es:", sp.simplify(sp.expand(vsv_at - vsv_es)))

    # direct speed from gamma_loc
    beta_s = sp.sqrt(beta2)
    beta = sp.sqrt(1 - 1 / gamma_loc ** 2)
    vratio = sp.series(beta_s / beta, eps, 0, 3).removeO().subs(eps, 1)
    vratio = sp.expand(vratio)
    print("direct v_s/v - vsv_es:", sp.simplify(sp.together(vratio - vsv_es)))

    # psi = 0 reduces to the magnetic coefficient (3*g-1)/(2*g**2)
    mag_quad = sp.simplify(sp.expand(vsv_es.subs(psi, 0)).coeff(delta ** 2))
    print("psi=0 quad coeff - magnetic:", sp.simplify(mag_quad - (3 * g - 1) / (2 * g ** 2)))

    # gamma(delta_loc) reproduces gamma_loc through order 2 (consistency)
    gamma_from_dl = sp.expand(gamma_of.subs(dl, claimed_dl))
    gamma_from_dl = sum(
        term
        for term in sp.Add.make_args(gamma_from_dl)
        if sp.total_degree(term, delta, psi) <= 2
    )
    print("gamma(delta_loc) - gamma_loc:", sp.simplify(sp.expand(gamma_from_dl - claimed_gamma)))


if __name__ == "__main__":
    main()
