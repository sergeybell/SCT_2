"""NOT USED for the audit coefficients.

This file expands a Hamiltonian in which A_s is the straight multipole
polynomial times (1+h x). That normalization does not reproduce the
sector weak focusing x''+(K_1+h**2)x = h*delta (the curvilinear curl
was not built into A_s). The sector sources of S2 that the report quotes
are computed in sym_D2_sector.py from the polar equation of motion.
The thin-lens bookkeeping below the Hamiltonian block remains valid and
agrees with sym_D2_sector.py on the sextupole piece -K_2 D_1**2/2.
"""
import sympy as sp

h, K1, K2 = sp.symbols("h K_1 K_2", real=True)
x, px, delta = sp.symbols("x p_x delta", real=True)
xp = sp.symbols("xprime", real=True)

S = sp.sqrt((1 + delta) ** 2 - px**2)
Poly = h * x + K1 * x**2 / 2 + K2 * x**3 / 6
H = -(1 + h * x) * S + (1 + h * x) * Poly

# reference orbit: p_x' = 0 at x=px=delta=0
px_dot = sp.simplify(-sp.diff(H, x))
print("p_x' at 0:", sp.simplify(px_dot.subs({x: 0, px: 0, delta: 0})))

x_dot = sp.diff(H, px)
# x'' = d/ds (x') = Xx x' + Xp px'
Xx = sp.diff(x_dot, x)
Xp = sp.diff(x_dot, px)
xpp_exact = sp.simplify(Xx * x_dot + Xp * px_dot)

# Expand the second-order ODE in amplitudes. Order counting:
# x, px, delta are first order; keep through second order in the EOM
# after it is written as x'' = F.
xpp_series = sp.series(xpp_exact, x, 0, 3).removeO()
xpp_series = sp.series(xpp_series, px, 0, 3).removeO()
xpp_series = sp.series(xpp_series, delta, 0, 3).removeO()
xpp_series = sp.expand(xpp_series)

# Invert x' = x_dot(x, px, delta) for px, through second order.
xd = sp.series(sp.series(sp.series(x_dot, x, 0, 3).removeO(), px, 0, 3).removeO(), delta, 0, 3).removeO()
xd = sp.expand(xd)
# px = a1 xp + quadratic form. Use successive substitution.
# At linear order xd = px * (coefficient). 
lin = sp.expand(xd.subs({x: 0, delta: 0}))
# lin should be px * c
c_lin = sp.simplify(lin / px)
print("linear x' / p_x:", c_lin)

# px = xp + quadratic corrections. Linear response of x' to p_x is 1.

# Let eps scale x, px, delta, xp.
eps = sp.symbols("eps")
xd_e = sp.expand(xd.subs({x: eps * x, px: eps * px, delta: eps * delta}))
# revert: set px = b1 + eps*b2 with b1 linear in xp, x, delta.
b1 = xp  # since c_lin should be 1 at delta=x=0; check
# We'll solve coefficients of the polynomial identity xd - xp == 0.

px_ans = xp
for order in range(1, 4):
    trial = sp.expand(xd.subs(px, px_ans) - xp)
    # residual starts at the next order; collect homogeneous pieces of degree `order`
    # scale
    res = sp.expand(trial.subs({x: eps * x, xp: eps * xp, delta: eps * delta}))
    res = sp.series(res, eps, 0, order + 1).removeO()
    # piece exactly degree `order`
    piece = sp.expand(res.coeff(eps, order))
    # piece is a polynomial in x, xp, delta. It equals (dxd/dpx at 0) * correction
    # dx'/dpx at 0 is 1, but at finite delta the linear map is px |-> px/(1+delta) roughly.
    # correction to px that cancels `piece`: since xd = px/(1+delta) + higher,
    # at each new order the linear response of xd to dpx is 1 (because lower orders
    # already match and the derivative at the current point...).
    # Use: d(xd)/d(px) at px=px_ans, linearized, leading term is 1 at eps=0.
    # Simplest reliable method: subtract piece from px_ans, then re-expand,
    # because xd = px + (terms that vanish at px=x=delta=0) + quadratic,
    # and the coefficient of the newly added px correction in xd is 1 + O(eps).
    px_ans = sp.expand(px_ans - piece)
    # re-substitute once so products of corrections are included at higher order
    px_ans = sp.series(px_ans.subs({x: eps * x, xp: eps * xp, delta: eps * delta}), eps, 0, 4).removeO()
    px_ans = sp.expand(px_ans.subs(eps, 1))

residual = sp.expand(xd.subs(px, px_ans) - xp)
residual = sp.series(residual.subs({x: eps * x, xp: eps * xp, delta: eps * delta}), eps, 0, 4).removeO()
print("inversion residual through eps^3:", sp.simplify(residual))
print("px(xp,x,delta) =")
print(sp.factor(px_ans))

# Now x'' expressed in x, xp, delta
xpp_x = sp.expand(xpp_series.subs(px, px_ans))
xpp_x = sp.series(xpp_x.subs({x: eps * x, xp: eps * xp, delta: eps * delta}), eps, 0, 4).removeO()
xpp_x = sp.expand(xpp_x.subs(eps, 1))
print("\n=== x'' through second order ===")
# write x'' + (K1+h**2)*x - h*delta = quadratic sources
lin_op = sp.expand((K1 + h**2) * x - h * delta)
rhs_quad = sp.expand(xpp_x + (K1 + h**2) * x - h * delta)
# rhs_quad should contain only quadratic terms (the equation is x'' = - (K1+h**2)x + h delta + quadratic)
# Check linear part of x'' + (K1+h**2)x - h delta vanishes
lin_check = sp.expand(rhs_quad.subs({
    x: eps * x, xp: eps * xp, delta: eps * delta
}))
print("linear piece of residual (expect 0):", sp.expand(sp.series(lin_check, eps, 0, 2).removeO()))
quad = sp.expand(sp.series(lin_check, eps, 0, 3).removeO() - sp.series(lin_check, eps, 0, 2).removeO())
quad = sp.expand(quad.subs(eps, 1))
print("quadratic residual x''+(K1+h^2)x - h delta =")
for term in sp.Add.make_args(sp.expand(quad)):
    print("  ", term)

# Closed orbit x = D1*delta + D2*delta**2, xp = Dp1*delta + Dp2*delta**2
# x'' = D1''*delta + D2''*delta**2
# At order delta: D1'' + (K1+h**2) D1 - h = 0
# At order delta**2: D2'' + (K1+h**2) D2 = S2
D1, D2, Dp1 = sp.symbols("D_1 D_2 Dp_1", real=True)
# xp on orbit to the order needed in quadratic terms: xp enters quadratically,
# so xp = D1' * delta is enough ( (D1' delta)^2 is delta**2; D2' would be delta**3 in squares).
x_co = D1 * delta
xp_co = Dp1 * delta
quad_co = sp.expand(quad.subs({x: x_co, xp: xp_co}))
# S2 is minus? Equation:
# x'' + (K1+h**2) x - h delta = quad(x,xp,delta)
# D2'' delta**2 + (K1+h**2) D2 delta**2 = coefficient of delta**2 in quad
S2 = sp.simplify(sp.expand(quad_co).coeff(delta, 2))
print("\n=== S2 (D2'' + (K1+h^2) D2 = S2) ===")
print(S2)
S2e = sp.expand(S2)
# collect structures
print("S2 expanded terms:")
for term in sp.Add.make_args(S2e):
    print("  ", term)

sext = sp.simplify(S2.coeff(K2))
print("coeff of K2 in S2 (times K2 gives the piece):", sext)
print("sextupole piece:", sp.simplify(K2 * sext))

# Also split S2 = S2_sext + S2_other
S2_sext = sp.simplify(-(K2 / 2) * D1**2)
print("S2 - (-K2/2 D1^2)  [other sources]:")
other = sp.expand(S2 - S2_sext)
print(other)

# --- thin-bend model of the note: drop body factors of h in the Hamiltonian
# kinetic and potential, keep a pure curvature drive by replacing the dipole
# with the straight-coordinate equation used in eq:hill-d.
# We still RECORD the sector result above. The thin model is expanded directly
# from the note's equation, to show the sextupole piece is unchanged.

d = delta
xx = x
# eq:hill-d with h = 1/rho: x'' + K1/(1+d) x = -K2/(2(1+d)) x**2 + h*d
# Move to x'' + K1 x = S, expand to delta**2 on x = D1 d + D2 d**2.
left_extra = sp.series(K1 / (1 + d) * (D1 * d + D2 * d**2), d, 0, 3).removeO()
# x'' + K1 x = K1 x - K1/(1+d) x - K2/(2(1+d)) x**2 + h d
rhs = sp.series(
    K1 * (D1 * d + D2 * d**2) - K1 / (1 + d) * (D1 * d + D2 * d**2)
    - K2 / (2 * (1 + d)) * (D1 * d) ** 2
    + h * d,
    d, 0, 3,
).removeO()
# order d: D1'' + K1 D1 = h, so the d^1 part of rhs must be h*d
# order d**2: D2'' + K1 D2 = coeff
print("\n=== thin model (note eq:hill-d) ===")
print("order delta of RHS:", sp.expand(rhs.coeff(d, 1)))
S2_thin = sp.expand(rhs.coeff(d, 2))
print("S2_thin:", S2_thin)
print("S2_thin sextupole:", sp.simplify(S2_thin.coeff(K2) * K2))
print("S2_thin other:", sp.expand(S2_thin + (K2 / 2) * D1**2))

# If the drive is the sector-consistent h*delta/(1+delta) instead of h*delta:
rhs2 = sp.series(
    K1 * (D1 * d + D2 * d**2) - K1 / (1 + d) * (D1 * d + D2 * d**2)
    - K2 / (2 * (1 + d)) * (D1 * d) ** 2
    + h * d / (1 + d),
    d, 0, 3,
).removeO()
print("S2 if drive h*delta/(1+delta):", sp.expand(rhs2.coeff(d, 2)))

# --- path length -> Delta alpha_1
# C/Cs = 1 + a0 d + a1 d**2
# (C-Cs)/Cs = (1/Cs) * int (x/rho + xp**2/2) ds
# x = D1 d + D2 d**2, xp = D1' d + ...
# a0 = (1/Cs) int D1 h ds
# a1 = (1/Cs) int (D2 h + Dp1**2 / 2) ds
# Sextupole changes only D2. Delta a1 = (1/Cs) int (Delta D2) h ds
# Green: int (Delta D2) h ds = int S2_sext * D1 ds
# if (d^2+Kx) Delta D2 = S2_sext and (d^2+Kx) D1 = h, kernel symmetric.
print("\n=== sextupole shift ===")
print("Delta alpha_1 = Delta eta_1 = (1/Cs) * oint S2_sext * D1 ds")
print("S2_sext =", S2_sext)
print("oint S2_sext D1 = -1/2 oint K2 D1**3")
print("Delta alpha_1 = -1/(2 Cs) oint K2 D1**3 ds")

# Show the quadrupole-feeddown piece does NOT come from K2
# and the slope piece h*Dp**2 is real in the sector model.
print("\nother sector sources (no K2):")
print(sp.collect(sp.expand(other), [D1, Dp1, h, K1]))
