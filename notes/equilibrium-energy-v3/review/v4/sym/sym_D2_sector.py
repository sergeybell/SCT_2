"""S2 from the median-plane sector equation (polar geometry), through order delta**2.

Exact local equation for constant curvature h (uniform sector), x outward:
  x'' = h*xp**2/(1+h*x) + h + h**2*x
        - b(x)/(1+delta) * (1+h*x) * sqrt((1+h*x)**2 + xp**2)
  b(x) = h + K1*x + (K2/2)*x**2
with b = -q B /(p0), so b(0) = h and K1 = (1/B rho) dB/dx in the note's sign.

Checked limits:
  equilibrium x=xp=delta=0 => x''=0
  linear: x'' + (K1+h**2)*x = h*delta
"""
import sympy as sp

h, K1, K2 = sp.symbols("h K_1 K_2", real=True)
x, xp, delta = sp.symbols("x xprime delta", real=True)
D1, D2, Dp1 = sp.symbols("D_1 D_2 Dp_1", real=True)

b = h + K1 * x + K2 * x**2 / 2
expr = (
    h * xp**2 / (1 + h * x)
    + h
    + h**2 * x
    - b / (1 + delta) * (1 + h * x) * sp.sqrt((1 + h * x) ** 2 + xp**2)
)

# equilibrium
print("eq at 0:", sp.simplify(expr.subs({x: 0, xp: 0, delta: 0})))

eps = sp.symbols("eps")
ex = sp.series(expr.subs({x: eps * x, xp: eps * xp, delta: eps * delta}), eps, 0, 3).removeO()
ex = sp.expand(ex.subs(eps, 1))

# linear operator expected: x'' + (K1+h**2)*x - h*delta = quadratic
residual = sp.expand(ex + (K1 + h**2) * x - h * delta)
lin = sp.expand(sp.series(residual.subs({x: eps*x, xp: eps*xp, delta: eps*delta}), eps, 0, 2).removeO())
print("linear residual (expect 0):", sp.simplify(lin))

quad = sp.expand(
    sp.series(residual.subs({x: eps*x, xp: eps*xp, delta: eps*delta}), eps, 0, 3).removeO()
    - sp.series(residual.subs({x: eps*x, xp: eps*xp, delta: eps*delta}), eps, 0, 2).removeO()
).subs(eps, 1)
quad = sp.expand(quad)
print("quadratic form of x''+(K1+h^2)x-h*delta:")
for term in sp.Add.make_args(quad):
    print("  ", term)

# orbit x = D1*delta, xp = Dp1*delta (enough for delta**2)
S_total = sp.expand(quad.subs({x: D1 * delta, xp: Dp1 * delta})).coeff(delta, 2)
# Equation: x'' + (K1+h**2) x - h delta = quad
# so D2'' + (K1+h**2) D2 = S2, with S2 = coeff of delta**2 in quad
S2 = sp.expand(S_total)
print("\nS2 =")
for term in sp.Add.make_args(S2):
    print("  ", term)

S2_sext = -K2 * D1**2 / 2
other = sp.expand(S2 - S2_sext)
print("\nsextupole piece:", S2_sext)
print("coeff check K2:", sp.simplify(S2.coeff(K2)))
print("other sources:")
print(sp.factor(other))
print(sp.collect(other, [D1, Dp1]))

# Group other sources by mechanism
print("\n--- grouped ---")
# terms with K2 already removed
# K1*h*D1**2 : quadrupole in curved frame
# h**3 * D1**2 : pure geometry of sector field expansion / path
# K1*D1 : 1/(1+delta) feed-down of quadrupole
# h**2*D1 : 1/(1+delta) feed-down of weak focusing
# h*Dp**2 : slope geometry
# h : drive expansion of delta/(1+delta) and curvature

# Thin-bend bookkeeping of the note (eq:hill-d): no h**2 in K, drive h*delta
d = delta
rhs = sp.series(
    K1 * (D1 * d) - K1 / (1 + d) * (D1 * d + D2 * d**2)
    - (K2 / 2) / (1 + d) * (D1 * d) ** 2
    + h * d,
    d, 0, 3,
).removeO()
# This rhs is for x''+K1 x. Order d**2 coefficient is S2, and D2 term on the left
# was moved: K1*D2 is cancelled by the expansion only through the -K1/(1+d)*D2 piece.
# Explicit:
print("\nthin-model S2 (note hill-d):", sp.expand(rhs.coeff(d, 2)))

# Green identity is algebraic once S2_sext = -K2/2 D1**2 and the linear
# operator is the same self-adjoint Hill operator for D1 and D2.
print("\nDelta alpha_1 from sextupole = (1/C_s) * oint S2_sext * D1 ds")
print("= -1/(2 C_s) * oint K2 * D1**3 ds")
print("same for Delta eta_1, kinematic eta1 has no K2")
