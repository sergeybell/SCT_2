# Linear electrostatic orbit equation, analyzer limits, edge jump order.
# Exact symbols only: n=1 and gamma->1 are the classical limits, not a fit.
import sympy as sp

x, rho, delta, n, g = sp.symbols("x rho delta n g", positive=True)
xp = sp.symbols("xprime", real=True)


def main():
    # pv / (p_s v_s) = 1 + (1 + 1/g) * (delta - x/rho)     linear
    # 1/R = (1/rho) / (pv/(p_s v_s)) * (1 - n*x/rho)
    gfac = 1 + 1 / g
    pv_ratio = 1 + gfac * (delta - x / rho)
    inv_R = (1 / rho) * (1 / pv_ratio) * (1 - n * x / rho)
    eps = sp.symbols("eps")
    inv_R_s = inv_R.subs({x: eps * x, delta: eps * delta})
    inv_R_lin = sp.series(inv_R_s, eps, 0, 2).removeO().subs(eps, 1)
    inv_R_lin = sp.expand(inv_R_lin)

    # x'' + x/rho**2 = 1/rho - 1/R
    # equilibrium of the linear ODE: coeff * x = source * delta, with x'' expressed
    rhs = sp.expand(1 / rho - inv_R_lin)
    # rhs = A * delta + B * x
    A = sp.simplify(rhs.coeff(delta))
    B = sp.simplify(rhs.coeff(x))
    # move B*x to the left: x'' - B*x = A*delta, but the geometric x/rho**2 is NOT
    # yet in rhs. Full left coefficient of x is 1/rho**2 - B.
    # Equation: x'' + x/rho**2 - rhs = 0, and rhs contains B*x + A*delta
    # x'' + (1/rho**2 - B) x = A delta
    K = sp.simplify(1 / rho ** 2 - B)
    source = sp.simplify(A)
    print("K:", K)
    print("source:", source)
    K_claimed = (2 - n + 1 / g) / rho ** 2
    S_claimed = gfac / rho
    print("K diff:", sp.simplify(K - K_claimed))
    print("source diff:", sp.simplify(source - S_claimed))

    # equivalent literature form  (3 - n - beta**2)/rho**2
    beta2 = 1 - 1 / g
    K_lit = (3 - n - beta2) / rho ** 2
    print("K vs 3-n-beta^2:", sp.simplify(K_claimed - K_lit))

    # phase advance per bend angle: d^2 x / d theta^2 + (K rho**2) x = 0
    # first-order image when sqrt(K)* (rho*theta) = pi, theta = pi/sqrt(K rho**2)
    nu2 = sp.simplify(K_claimed * rho ** 2)
    theta = sp.pi / sp.sqrt(nu2)
    theta_cyl = sp.simplify(theta.subs({n: 1, g: 1}))
    theta_sph = sp.simplify(theta.subs({n: 2, g: 1}))
    print("n=1 nonrel angle:", theta_cyl, "minus pi/sqrt(2):", sp.simplify(theta_cyl - sp.pi / sp.sqrt(2)))
    print("n=2 nonrel angle:", theta_sph, "minus pi:", sp.simplify(theta_sph - sp.pi))
    # ultra-relativistic cylindrical: g -> oo, n=1 -> theta = pi
    theta_rel = sp.simplify(sp.limit(theta.subs(n, 1), g, sp.oo))
    print("n=1 ultra-rel angle minus pi:", sp.simplify(theta_rel - sp.pi))

    # circular-orbit dispersion from force balance, linear
    # D/rho = source / (K) / delta * rho wait: x = (source/K) * delta if x''=0
    D_over = sp.simplify((S_claimed / K_claimed) / rho)
    print("D1/rho:", D_over)
    print("n=1 D1/rho - 1:", sp.simplify(D_over.subs(n, 1) - 1))
    print("n=2 D1/rho - (g+1):", sp.simplify(D_over.subs(n, 2) - (g + 1)))

    # edge: p_x continuous, p_parallel jumps with delta_loc
    # x'_in = x'_out * p_out / p_in, p_out = 1+delta, p_in = 1+delta_loc
    psi = sp.symbols("psi", real=True)
    dl = delta - psi + (2 * delta * psi - psi ** 2) / (2 * g)
    p_out = 1 + delta
    p_in = 1 + dl
    eps = sp.symbols("eps")
    p_out_s = 1 + eps * delta
    dl_s = dl.subs({delta: eps * delta, psi: eps * psi})
    p_in_s = 1 + dl_s
    jump = xp * (p_out_s / p_in_s - 1)
    jump_s = sp.expand(sp.series(jump, eps, 0, 3).removeO().subs(eps, 1))
    print("angle jump through order 2:", jump_s)
    print("jump - xprime*psi:", sp.simplify(sp.expand(jump_s - xp * psi)))


if __name__ == "__main__":
    main()
