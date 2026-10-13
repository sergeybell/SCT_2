"""Slip-factor definitions: Ng (18.17)--(18.21), Senichev (3)--(5), variable changes.

Series coefficients are the multipliers in
    (T - T_s)/T_s = eta0*delta + eta1*delta**2 + O(delta**3),
    C/C_s = 1 + alpha0*delta + alpha1*delta**2 + O(delta**3).
A derivative is not divided by the factorial. Run: python ng_verify.py
"""
import sympy as sp

d = sp.symbols("delta", real=True)
g = sp.symbols("gamma_s", positive=True)
a0, a1, a1ng = sp.symbols("alpha0 alpha1 alpha1_Ng", real=True)
b2 = 1 - 1 / g**2
checks = []


def check(name, expr):
    simplified = sp.simplify(expr)
    ok = simplified == 0
    checks.append(ok)
    print(f"{'OK' if ok else 'FAIL'}  {name}")
    if not ok:
        print("     ", simplified)


def coeff(expr, n):
    return sp.series(expr, d, 0, n + 1).removeO().coeff(d, n)


print("=== kinematic series of tau = (T-T_s)/T_s ===")
bg_s = sp.sqrt(g**2 - 1)
beta = bg_s * (1 + d) / sp.sqrt(1 + bg_s**2 * (1 + d) ** 2)
beta_s = bg_s / g
C_over_Cs = 1 + a0 * d + a1 * d**2
T_over_Ts = (C_over_Cs) * (beta_s / beta)
eta0 = sp.simplify(coeff(T_over_Ts - 1, 1))
eta1 = sp.simplify(coeff(T_over_Ts - 1, 2))
eta0_ref = a0 - 1 / g**2
eta1_ref = a1 - a0 / g**2 + 1 / g**4 + sp.Rational(3, 2) * b2 / g**2
print("eta0 =", eta0)
print("eta1 =", eta1)
check("eta0 matches section 01", eta0 - eta0_ref)
check("eta1 matches section 01", eta1 - eta1_ref)

print("\n=== Ng (18.15): beta derivatives at delta=0 ===")
beta_ratio = sp.simplify(beta / beta_s)
# (18.15) quotes (1/beta) d^n beta / d delta^n, which equals the same for beta/beta_s.
bp = sp.diff(beta_ratio, d)
bpp = sp.diff(beta_ratio, d, 2)
bppp = sp.diff(beta_ratio, d, 3)
check("beta'/beta = 1/gamma^2", bp.subs(d, 0) - 1 / g**2)
check("beta''/beta = -3 beta^2/gamma^2", bpp.subs(d, 0) - (-3 * b2 / g**2))
check(
    "beta'''/beta = -3 beta^2 (1-5 beta^2)/gamma^2",
    bppp.subs(d, 0) - (-3 * b2 * (1 - 5 * b2) / g**2),
)

print("\n=== Ng (18.17), (18.18) after alpha1 = alpha0 * alpha1_Ng ===")
# (18.18): eta1 = alpha0*alpha1_Ng + 3 beta^2/(2 gamma^2) - eta0/gamma^2
ng_eta1 = a0 * a1ng + sp.Rational(3, 2) * b2 / g**2 - eta0_ref / g**2
check("alpha1 = alpha0*alpha1_Ng turns (18.18) into eta1", ng_eta1.subs(a1ng, a1 / a0) - eta1_ref)

print("\n=== (18.20) and (18.21) from omega/omega_s = T_s/T ===")
e0, e1 = sp.symbols("eta0 eta1", real=True)
omega_ratio = 1 / (1 + e0 * d + e1 * d**2)  # omega/omega_s
# Delta omega/omega_s = -eta * delta, eta = eta0 + eta1_1820 * delta + ...
slip_1820 = sp.series(-(omega_ratio - 1) / d, d, 0, 2).removeO()
eta1_1820 = sp.expand(slip_1820.coeff(d, 1))
# (18.21): eta = -(1/omega_s) d omega / d delta
slip_1821 = sp.series(-sp.diff(omega_ratio, d), d, 0, 2).removeO()
eta1_1821 = sp.expand(slip_1821.coeff(d, 1))
print("(18.20) coefficient of delta in eta =", eta1_1820)
print("(18.21) coefficient of delta in eta =", eta1_1821)
check("(18.20) gives eta1 - eta0**2", eta1_1820 - (e1 - e0**2))
check("(18.21) gives 2(eta1 - eta0**2)", eta1_1821 - 2 * (e1 - e0**2))

# Same quadratic from the particle's own period in the denominator.
own = (e0 * d + e1 * d**2) / (1 + e0 * d + e1 * d**2)  # (T-T_s)/T
own_quad = sp.expand(coeff(own, 2))
check("(T-T_s)/T has quadratic eta1 - eta0**2", own_quad - (e1 - e0**2))

print("\n=== error in delta_eq if that eta1 is used in <tau>=0 ===")
dm = sp.symbols("delta_m", real=True)
deq = -e1 * dm**2 / (2 * e0)
deq_1820 = -(e1 - e0**2) * dm**2 / (2 * e0)
deq_1821 = -2 * (e1 - e0**2) * dm**2 / (2 * e0)
print("delta_eq(18.20) - delta_eq =", sp.factor(deq_1820 - deq))
print("delta_eq(18.21) - delta_eq =", sp.factor(deq_1821 - deq))
check("18.20 shifts delta_eq by eta0*delta_m**2/2", deq_1820 - deq - e0 * dm**2 / 2)
check(
    "18.21 shifts delta_eq by delta_eq + eta0*delta_m**2",
    deq_1821 - deq - (deq + e0 * dm**2),
)

print("\n=== Senichev (5): denominator T_s, velocity cut at delta/gamma^2 ===")
u = d / g**2  # Delta v/v_s kept to first order only
sen = sp.expand(
    (1 + a0 * d + a1 * d**2) * (1 - u + u**2) - 1
)  # (T-T_s)/T_s with that cut
sen_quad = sp.collect(sp.series(sen, d, 0, 3).removeO(), d).coeff(d, 2)
sen_ref = a1 - a0 / g**2 + 1 / g**4
print("Senichev quadratic =", sp.simplify(sen_quad))
check("Senichev quadratic is alpha1 - alpha0/gamma^2 + 1/gamma^4", sen_quad - sen_ref)
check("missing piece is 3 beta^2/(2 gamma^2), not eta0**2", eta1_ref - sen_ref - sp.Rational(3, 2) * b2 / g**2)
# (T-T_s)/T with the FULL velocity carries eta1 - eta0**2, not Senichev's bracket.
full_own_quad = sp.simplify(eta1_ref - eta0_ref**2)
diff_den = sp.simplify(full_own_quad - sen_ref)
print("full (T-T_s)/T quadratic - Senichev =", diff_den)
diff_den_num = sp.simplify(
    diff_den.subs({a0: sp.Rational(1, 5), a1: sp.Rational(1, 3), g: sp.Rational(5, 4)})
)
check(
    "(T-T_s)/T quadratic differs from Senichev (5)",
    sp.Integer(0 if diff_den_num != 0 else 1),
)

print("\n=== half the second derivative is the series coefficient ===")
# Generic: f = 1 + c1 d + c2 d**2, c2 = (1/2) f''(0).
c1, c2 = sp.symbols("c1 c2")
f = 1 + c1 * d + c2 * d**2
check("series coeff of delta**2 is half the second derivative", c2 - sp.diff(f, d, 2).subs(d, 0) / 2)

print("\n=== COSY delta_K and the scaled time ===")
# delta_K = (gamma-gamma_s)/(gamma_s-1), gamma^2 = 1+(gamma_s^2-1)*(1+delta)**2
gam = sp.sqrt(1 + (g**2 - 1) * (1 + d) ** 2)
dk = (gam - g) / (g - 1)
dk1 = sp.simplify(coeff(dk, 1))
dk2 = sp.simplify(coeff(dk, 2))
print("delta_K linear =", dk1)
print("delta_K quadratic =", dk2)
check("d delta_K/d delta = (gamma+1)/gamma", dk1 - (g + 1) / g)
check("delta_K quadratic = (gamma+1)/(2 gamma**3)", dk2 - (g + 1) / (2 * g**3))

# l = e1*delta_K + e2*delta_K**2, tau = -(gamma+1)/gamma * l/C_s
e1m, e2m, Cs = sp.symbols("e1 e2 C_s", real=True)
l_over_C = (e1m * dk + e2m * dk**2) / Cs
tau_cosy = -(g + 1) / g * l_over_C
eta0_map = sp.simplify(coeff(tau_cosy, 1))
eta1_map = sp.simplify(coeff(tau_cosy, 2))
print("eta0 from (l, delta_K) =", eta0_map)
print("eta1 from (l, delta_K) =", eta1_map)
# Doubling the monomial e2 (reading a derivative as a coefficient) does not double eta1.
eta1_doubled = eta1_map.subs(e2m, 2 * e2m)
double_num = sp.simplify(
    (eta1_doubled - 2 * eta1_map).subs({e1m: 1, e2m: 1, Cs: 1, g: sp.Rational(5, 4)})
)
print("eta1(2*e2) - 2*eta1 at a sample point =", double_num)
check(
    "doubling (5|66) does not double eta1",
    sp.Integer(0 if double_num != 0 else 1),
)

print("\n=== MAD-X pt = (E-E0)/(p0 c) ===")
# c=1 units of chapter 39: pt = (gamma-gamma_s)/(beta_s*gamma_s)
pt = (gam - g) / (bg_s)
pt1 = sp.simplify(coeff(pt, 1))
pt2 = sp.simplify(coeff(pt, 2))
print("pt linear =", pt1)
print("pt quadratic =", pt2)
check("d pt/d delta = beta_s", pt1 - sp.sqrt(b2))
# Re-expansion: Delta C/C = A*pt + B*pt**2 changes alpha1.
A, B = sp.symbols("A B", real=True)
path_pt = A * pt + B * pt**2
alpha0_from_pt = sp.simplify(coeff(path_pt, 1))
alpha1_from_pt = sp.simplify(coeff(path_pt, 2))
print("alpha0 from a pt-series =", alpha0_from_pt)
print("alpha1 from a pt-series =", alpha1_from_pt)
pt_quad_num = sp.simplify(
    (alpha1_from_pt - B).subs({A: 1, B: 1, g: sp.Rational(5, 4)})
)
print("alpha1(from pt) - B at a sample point =", pt_quad_num)
check(
    "pt re-expansion changes the quadratic coefficient",
    sp.Integer(0 if pt_quad_num != 0 else 1),
)
# ALPHA_C_P is the second derivative, so the series coefficient is ALPHA_C_P/2.
ac, acp = sp.symbols("ALPHA_C ALPHA_C_P", real=True)
path_deriv = ac * d + acp * d**2 / 2
check("ALPHA_C_P/2 is the series alpha1", coeff(path_deriv, 2) - acp / 2)

print("\n=== using ALPHA_C_P as alpha1 shifts eta1 by the true alpha1 ===")
eta1_if_deriv = eta1_ref.subs(a1, acp)  # wrongly
# true alpha1 = acp/2, true eta1 uses a1=acp/2
shift = sp.simplify(eta1_if_deriv.subs(acp, 2 * a1) - eta1_ref)
print("eta1(alpha1=ALPHA_C_P) - eta1 =", shift)
check("that shift equals the true alpha1", shift - a1)

print()
if all(checks):
    print(f"ALL {len(checks)} CHECKS PASSED")
else:
    print(f"{sum(checks)}/{len(checks)} passed")
    raise SystemExit(1)
