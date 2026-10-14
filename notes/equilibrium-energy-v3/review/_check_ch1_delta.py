# Short rational checks for chapter 1 (delta). Each identity is exact arithmetic.
import sympy as sp

delta = sp.symbols("delta")
g2 = sp.Rational(9, 4)  # gamma_s**2
g_ratio_sq = (1 + (g2 - 1) * (1 + delta) ** 2) / g2
vsv = sp.series(sp.sqrt(g_ratio_sq) / (1 + delta), delta, 0, 3).removeO()
vsv_f = 1 - delta / g2 + (3 * g2 - 1) / (2 * g2**2) * delta**2
assert sp.expand(vsv - vsv_f) == 0

gs = sp.sqrt(g2)
gamma = sp.series(gs * sp.sqrt(g_ratio_sq), delta, 0, 3).removeO()
gamma_f = (
    gs
    + gs * (1 - 1 / g2) * delta
    + (g2 - 1) / (2 * g2 * gs) * delta**2
)
assert sp.expand(gamma - gamma_f) == 0

# gamma_s**2 -> 1: v_s/v -> 1/(1+delta)
v_nr = sp.series(1 / (1 + delta), delta, 0, 3).removeO()
v_formula_nr = 1 - delta + (3 * 1 - 1) / 2 * delta**2
assert sp.expand(v_nr - v_formula_nr) == 0

# BM coefficients at fixed rationals
eta0 = sp.Rational(2)
eta1 = sp.Rational(3)
w2 = sp.Rational(5)
gg = sp.Rational(4)
a = sp.Rational(1, 5)
c2 = -(w2 / eta0) * (eta1 - eta0 / gg)
cd = -sp.Rational(1) / gg
s_plus = sp.together(c2 + cd * w2)
s_minus = sp.together(c2 - cd * w2)
assert s_plus == -w2 * eta1 / eta0
assert s_minus == -w2 * eta1 / eta0 + 2 * w2 / gg
D0 = sp.together((a**2 / 2) * s_plus / w2)
D2 = sp.together(-a**2 * s_minus / (6 * w2))
assert D0 == -eta1 * a**2 / (2 * eta0)
assert D2 == (a**2 / 6) * (eta1 / eta0 - 2 / gg)

# constant of w = delta + delta**2/(2 gamma**2) matches -eta_tilde * a**2 / (2 eta0)
eta_tilde = eta1 - eta0 / (2 * gg)
w_const = D0 + (a**2 / 2) / (2 * gg)
w_harm = D2 + (a**2 / 2) / (2 * gg)
assert sp.together(w_const - (-eta_tilde * a**2 / (2 * eta0))) == 0
assert sp.together(w_harm - (eta_tilde * a**2 / (6 * eta0))) == 0

# phase: integrate tau, constant term vanishes, sin(2psi) coefficient
# tau cos2 coeff = (2/3) eta1 a**2 - eta0 a**2 / (3 gg)
# phi sin2 coeff relative to omega_rf/omega = half of that
phi_sin2 = sp.Rational(1, 2) * ((sp.Rational(2, 3)) * eta1 * a**2 - eta0 * a**2 / (3 * gg))
phi_expected = (eta1 / 3) * a**2 - eta0 * a**2 / (6 * gg)
assert sp.together(phi_sin2 - phi_expected) == 0

# <delta^2> pieces: delta = D0 + a cos + D2 cos2
# mean of cos^2 = 1/2; cross harmonics 0; D0^2 and D2^2 are O(a^4)
assert sp.together(D0**2 / (a**2 / 2)) == sp.together(D0**2 * 2 / a**2)

print("all checks passed")
print("D0", D0)
print("D2", D2)
print("vsv", vsv)
print("gamma", sp.simplify(gamma))
