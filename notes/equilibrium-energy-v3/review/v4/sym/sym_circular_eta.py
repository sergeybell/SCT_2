# Second-order circular orbit by Taylor polynomials, not by series() of a root.
# gamma(delta) is the already checked magnetic series. g and n stay symbols.
import sympy as sp

g, n, delta = sp.symbols("g n delta", positive=True)
a, b = sp.symbols("a b")


def main():
    gamma_s = sp.sqrt(g)
    k = (g - 1) / gamma_s  # beta**2 * gamma
    x = a * delta + b * delta ** 2
    psi = x - (n / 2) * a ** 2 * delta ** 2
    # gamma_loc = gamma_out - k*psi, magnetic series for gamma_out
    gamma_loc = (
        gamma_s
        + k * delta
        + k * delta ** 2 / (2 * g)
        - k * psi
    )
    # keep total degree <= 2
    gamma_loc = sp.expand(gamma_loc)
    dgamma = sp.expand(gamma_loc - gamma_s)
    dgamma = _deg(dgamma, 2)
    # f = gamma - 1/gamma, f' = 1+1/g, f'' = -2/gamma**3
    fp = 1 + 1 / g
    fpp = -2 / gamma_s ** 3
    f = k + fp * dgamma + sp.Rational(1, 2) * fpp * _deg(sp.expand(dgamma ** 2), 2)
    f = _deg(sp.expand(f), 2)
    # pv/(m c**2) target = k * r**(1-n)
    pv = k * (
        1
        + (1 - n) * x
        + (1 - n) * (-n) / 2 * a ** 2 * delta ** 2
    )
    pv = _deg(sp.expand(pv), 2)
    residual = _deg(sp.expand(f - pv), 2)
    c1 = sp.together(sp.simplify(residual.coeff(delta)))
    c2 = sp.together(sp.simplify(residual.coeff(delta ** 2)))
    a_val = sp.factor(sp.together(sp.simplify(sp.solve(c1, a)[0])))
    b_val = sp.factor(
        sp.together(sp.simplify(sp.solve(sp.expand(c2.subs(a, a_val)), b)[0]))
    )
    a_ode = sp.factor(sp.together((1 + 1 / g) / (2 - n + 1 / g)))
    print("a:", a_val)
    print("a - ODE:", sp.simplify(sp.together(a_val - a_ode)))
    print("b:", b_val)

    # flight time from the checked vsv series and C/C_s = 1 + x
    vsv = (
        1
        - (delta - psi) / g
        + sp.Rational(3, 2) * (delta - psi) ** 2 / g
        - delta ** 2 / (2 * g ** 2)
    )
    vsv = vsv.subs({a: a_val, b: b_val})
    # psi still contains a, b; substitute before expanding in delta
    psi_s = psi.subs({a: a_val, b: b_val})
    vsv = (
        1
        - (delta - psi_s) / g
        + sp.Rational(3, 2) * (delta - psi_s) ** 2 / g
        - delta ** 2 / (2 * g ** 2)
    )
    path = 1 + x.subs({a: a_val, b: b_val})
    prod = _deg(sp.expand(path * vsv), 2)
    tau1 = sp.factor(sp.together(prod.coeff(delta)))
    tau2 = sp.factor(sp.together(prod.coeff(delta ** 2)))
    eta0 = (1 + 1 / g) * a_val - 1 / g
    eta1 = (
        (1 + 1 / g) * b_val
        + (5 - n) / (2 * g) * a_val ** 2
        - 4 * a_val / g
        + (3 * g - 1) / (2 * g ** 2)
    )
    print("eta0 - tau:", sp.simplify(sp.together(sp.expand(eta0 - tau1))))
    print("eta1 - tau:", sp.simplify(sp.together(sp.expand(eta1 - tau2))))

    a1 = sp.simplify(a_val.subs(n, 1))
    b1 = sp.factor(sp.together(b_val.subs(n, 1)))
    print("n=1 a-1:", sp.simplify(a1 - 1))
    print("n=1 b-(1+1/g)/2:", sp.simplify(sp.together(b1 - (1 + 1 / g) / 2)))
    a2 = sp.simplify(sp.together(a_val.subs(n, 2)))
    print("n=2 a-(g+1):", sp.simplify(sp.together(a2 - (g + 1))))
    eta0_2 = sp.simplify(sp.together(((1 + 1 / g) * a2 - 1 / g)))
    print("n=2 eta0-(g+2):", sp.simplify(sp.together(eta0_2 - (g + 2))))


def _deg(expr, order):
    expr = sp.expand(expr)
    return sum(
        t
        for t in sp.Add.make_args(expr)
        if sp.total_degree(t, delta) <= order
    )


if __name__ == "__main__":
    main()
