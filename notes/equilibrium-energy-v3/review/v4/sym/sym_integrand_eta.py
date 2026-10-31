# Symbolic identity: path * (v_s/v) -> eta0^E, eta1^E, and ch.6 rearrangements.
# Coefficients stay symbols. No numeric substitution.
import sympy as sp

delta = sp.symbols("delta")
d1, d2, p2, n, g = sp.symbols("d1 d2 p2 n g", real=True)
a0, a1, a0m, a1m = sp.symbols("alpha0 alpha1 alpha0_M alpha1_M", real=True)


def product_from_exact(d1, d2, p2, n, g):
    """Series of (dl/ds)*(v_s/v) from exact gamma and the quadratic potential."""
    eps = sp.symbols("eps")
    beta2 = 1 - 1 / g
    gamma_s = sp.sqrt(g)
    d = eps * delta
    x_over = d1 * d + d2 * d ** 2
    psi = x_over - (n / 2) * (d1 * d) ** 2
    gamma_out = sp.sqrt(1 + (g - 1) * (1 + d) ** 2)
    gamma_loc = gamma_out - beta2 * gamma_s * psi
    beta_s = sp.sqrt(beta2)
    beta = sp.sqrt(1 - 1 / gamma_loc ** 2)
    vsv = beta_s / beta
    # path: 1 + x/rho + x'^2/2, x' = D1' * delta + O(delta**2), x'^2 = p2 * delta**2
    path = 1 + x_over + sp.Rational(1, 2) * p2 * d ** 2
    prod = sp.series(path * vsv, eps, 0, 3).removeO().subs(eps, 1)
    prod = sp.expand(prod)
    c0 = sp.expand(prod.subs(delta, 0))
    c1 = sp.expand(prod.coeff(delta))
    c2 = sp.expand(prod.coeff(delta ** 2))
    return c0, c1, c2


def closed(d1, d2, p2, n, g):
    eta0 = (1 + 1 / g) * d1 - 1 / g
    eta1 = (
        p2 / 2
        + (1 + 1 / g) * d2
        + (5 - n) / (2 * g) * d1 ** 2
        - 4 * d1 / g
        + (3 * g - 1) / (2 * g ** 2)
    )
    return sp.expand(eta0), sp.expand(eta1)


def main():
    c0, c1, c2 = product_from_exact(d1, d2, p2, n, g)
    e0, e1 = closed(d1, d2, p2, n, g)
    print("c0 - 1:", sp.simplify(c0 - 1))
    print("c1 - eta0:", sp.simplify(sp.together(c1 - e0)))
    print("c2 - eta1:", sp.simplify(sp.together(c2 - e1)))

    # split of the linear coefficient: geometry, potential slowing, free kinematic
    geom = d1
    slowing = d1 / g
    free = -1 / g
    print("linear split:", sp.simplify(c1 - (geom + slowing + free)))

    # (5-n)/(2g) composition: -n/(2g) + 3/(2g) + 1/g
    comp = -n / (2 * g) + sp.Rational(3, 2) / g + 1 / g
    print("(5-n) composition:", sp.simplify(comp - (5 - n) / (2 * g)))

    # -4/g = -3/g from 3(delta-psi)**2 cross term, -1/g from path * (-delta/g)
    print("minus four split:", sp.simplify(-3 / g - 1 / g + 4 / g))

    # chapter 6: eta1^E = magnetic form on alpha^E + (1/g)*bracket
    # alpha0 = d1, alpha1 = d2 + p2/2   (local densities; averages are linear)
    alpha0 = d1
    alpha1 = d2 + p2 / 2
    mag_form = alpha1 - alpha0 / g + (3 * g - 1) / (2 * g ** 2)
    bracket = d2 + (5 - n) / 2 * d1 ** 2 - 3 * alpha0
    rewritten = sp.expand(mag_form + bracket / g)
    print("cmp-eta1 rewrite:", sp.simplify(sp.together(rewritten - e1)))

    # differences between rings. D-terms in the bracket are the ES ones.
    eta0_m = a0m - 1 / g
    eta0_e = (1 + 1 / g) * a0 - 1 / g
    diff0 = sp.expand(eta0_e - eta0_m)
    claimed0 = (a0 - a0m) + a0 / g
    print("cmp-eta0-diff:", sp.simplify(diff0 - claimed0))

    eta1_m = a1m - a0m / g + (3 * g - 1) / (2 * g ** 2)
    # use averaged symbols: <D2/rho> separate from alpha1
    d2_avg, d1sq = sp.symbols("D2_over_rho D1sq")
    # alpha1_E = <D2/rho> + <D'^2>/2, and e1 uses both
    # For the difference formula the slope sits inside alpha1 only.
    eta1_e = (
        a1
        - a0 / g
        + (3 * g - 1) / (2 * g ** 2)
        + (1 / g) * (d2_avg + (5 - n) / 2 * d1sq - 3 * a0)
    )
    diff1 = sp.expand(eta1_e - eta1_m)
    claimed1 = (
        (a1 - a1m)
        - (a0 - a0m) / g
        + (1 / g) * (d2_avg + (5 - n) / 2 * d1sq - 3 * a0)
    )
    print("cmp-eta1-diff:", sp.simplify(diff1 - claimed1))

    # outside the deflector d1=d2=0 reproduces magnetic vsv * (1 + x'^2/2)
    c0o, c1o, c2o = product_from_exact(0, 0, p2, n, g)
    mag_v = 1 - delta / g + (3 * g - 1) / (2 * g ** 2) * delta ** 2
    path_o = 1 + sp.Rational(1, 2) * p2 * delta ** 2
    mag_prod = sp.expand(sp.series(path_o * mag_v, delta, 0, 3).removeO())
    print(
        "outside match:",
        sp.simplify(c0o - mag_prod.subs(delta, 0)),
        sp.simplify(c1o - mag_prod.coeff(delta)),
        sp.simplify(sp.together(c2o - mag_prod.coeff(delta ** 2))),
    )


if __name__ == "__main__":
    main()
