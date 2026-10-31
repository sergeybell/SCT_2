# Laplace in the deflector frame and the quadratic potential.
import sympy as sp

rho, x, y, n = sp.symbols("rho x y n", real=True)
hs = 1 + x / rho


def lap_orthogonal(Phi):
    """Laplace operator with h_x = h_y = 1, h_s = 1 + x/rho, no s dependence."""
    term_x = sp.diff(hs * sp.diff(Phi, x), x)
    term_y = sp.diff(hs * sp.diff(Phi, y), y)
    return sp.simplify(term_x + term_y) / hs


def main():
    # div form used in the note, set to zero
    psi = (
        x / rho
        - n / 2 * (x / rho) ** 2
        + (n - 1) / 2 * (y / rho) ** 2
    )
    div = sp.diff(hs * sp.diff(psi, x), x) + sp.diff(hs * sp.diff(psi, y), y)
    # value at the origin fixes b
    b = sp.symbols("b")
    trial = x / rho - n / 2 * (x / rho) ** 2 + b * (y / rho) ** 2
    lap = (
        sp.diff(trial, x, 2)
        + sp.diff(trial, y, 2)
        + 1 / (rho + x) * sp.diff(trial, x)
    )
    at0 = sp.simplify(lap.subs({x: 0, y: 0}))
    print("Laplace at 0:", at0)
    print("b solved:", sp.solve(at0, b))

    # residual of the adopted quadratic is O(x), constant term vanishes
    lap_q = (
        sp.diff(psi, x, 2)
        + sp.diff(psi, y, 2)
        + 1 / (rho + x) * sp.diff(psi, x)
    )
    res = sp.series(lap_q.subs(y, 0), x, 0, 2).removeO()
    print("median residual through O(x):", sp.simplify(res))
    print("residual at 0:", sp.simplify(lap_q.subs({x: 0, y: 0})))

    # cylindrical n=1: psi = ln(1 + x/rho), no y
    psi_c = sp.log(1 + x / rho)
    series_c = sp.series(psi_c, x, 0, 3).removeO()
    claimed_c = x / rho - sp.Rational(1, 2) * (x / rho) ** 2
    print("cylinder minus claimed:", sp.simplify(sp.expand(series_c - claimed_c)))
    print("cylinder lap:", sp.simplify(lap_orthogonal(psi_c)))

    # spherical n=2: psi = 1 - rho/R, rho>0 removes Abs from the series
    rho_p = sp.symbols("rho_p", positive=True)
    R = sp.sqrt((rho_p + x) ** 2 + y ** 2)
    psi_s = 1 - rho_p / R
    series_s = sp.series(sp.series(psi_s, y, 0, 3).removeO(), x, 0, 3).removeO()
    claimed_s = x / rho_p - (x / rho_p) ** 2 + sp.Rational(1, 2) * (y / rho_p) ** 2
    # series_s may still contain mixed orders; subtract and drop O(3)
    diff = sp.expand(series_s - claimed_s)
    diff = sum(
        t
        for t in sp.Add.make_args(diff)
        if sp.total_degree(t, x, y) <= 2
    )
    print("sphere minus claimed:", sp.factor(sp.together(sp.simplify(diff))))


if __name__ == "__main__":
    main()
