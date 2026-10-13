# -*- coding: utf-8 -*-
"""Проверка равенств раздела sections/01_magnetic.tex.

Каждая проверка — тождество sympy. Ненулевой остаток печатается и даёт
ненулевой код выхода. Числа внизу — подстановка в уже доказанные ряды,
не параметры какого-либо кольца.
"""

from __future__ import annotations

import sys

import sympy as sp

delta = sp.symbols("delta", real=True)
gs = sp.symbols("gamma_s", positive=True)
a0 = sp.symbols("alpha_0", real=True)
a1 = sp.symbols("alpha_1", real=True)
K = sp.symbols("K", positive=True)
wrf = sp.symbols("omega_rf", positive=True)
eta_abs = sp.symbols("eta_abs", positive=True)
dm = sp.symbols("delta_m", real=True)

bs2 = 1 - 1 / gs**2  # beta_s**2
half = sp.Rational(1, 2)
three_half = sp.Rational(3, 2)

gam = sp.sqrt(1 + (gs**2 - 1) * (1 + delta) ** 2)
beta_over_betas = gs * (1 + delta) / gam
vs_over_v = gam / (gs * (1 + delta))

results: list[tuple[str, bool, object]] = []


def check(name: str, expr) -> None:
    residual = sp.simplify(sp.expand(expr))
    # series objects and equations
    if isinstance(residual, sp.Expr) and residual.has(sp.O):
        residual = sp.simplify(residual.removeO())
    ok = residual == 0
    results.append((name, ok, residual))
    status = "PASS" if ok else "FAIL"
    print(f"{status}  {name}")
    if not ok:
        print(f"      residual: {residual}")


def coeff(expr, n: int):
    return sp.expand(expr).series(delta, 0, n + 1).removeO().coeff(delta, n)


# --- точные кинематические тождества ---
# beta*gamma = beta_s*gamma_s*(1+delta), gamma**2 = 1+(beta*gamma)**2
beta_gamma = sp.sqrt(gs**2 - 1) * (1 + delta)
check("gamma**2 = 1+(beta*gamma)**2", gam**2 - (1 + beta_gamma**2))
check(
    "beta/beta_s = gamma_s(1+delta)/gamma",
    beta_over_betas - beta_gamma / gam * gs / sp.sqrt(gs**2 - 1),
)

# d ln beta / d ln p = 1/gamma**2
# p = p_s(1+delta), beta/beta_s = beta_over_betas
# d ln beta / d delta = d(ln beta_over_betas)/d delta, since beta_s fixed
# d ln p / d delta = 1/(1+delta)
dlnb_dd = sp.diff(sp.log(beta_over_betas), delta)
dlnp_dd = 1 / (1 + delta)
check(
    "d ln beta / d ln p = 1/gamma**2",
    sp.simplify(dlnb_dd / dlnp_dd - 1 / gam**2),
)

lam = 1 / (gam**2 * (1 + delta))
check("lambda = (1/beta) d beta/d delta", sp.simplify(lam - dlnb_dd))
check("lambda_s = 1/gamma_s**2", sp.simplify(lam.subs(delta, 0) - 1 / gs**2))

d_gam_inv = sp.diff(gam ** (-2), delta)
check(
    "d(gamma**(-2))/d delta = -2(gamma_s**2-1)(1+delta)/gamma**4",
    sp.simplify(d_gam_inv + 2 * (gs**2 - 1) * (1 + delta) / gam**4),
)
lam_prime = sp.diff(lam, delta)
check(
    "lambda' formula",
    sp.simplify(
        lam_prime
        - (d_gam_inv / (1 + delta) - 1 / (gam**2 * (1 + delta) ** 2))
    ),
)
lam_prime_s = sp.simplify(lam_prime.subs(delta, 0))
check(
    "lambda'_s = -2 beta_s**2/gamma_s**2 - 1/gamma_s**2",
    sp.simplify(lam_prime_s - (-2 * bs2 / gs**2 - 1 / gs**2)),
)
bpp = sp.simplify((lam.subs(delta, 0)) ** 2 + lam_prime_s)
check(
    "(beta''/beta)_s raw = 1/gamma**4 - 2 beta**2/gamma**2 - 1/gamma**2",
    sp.simplify(bpp - (1 / gs**4 - 2 * bs2 / gs**2 - 1 / gs**2)),
)
check(
    "(beta''/beta)_s = -3 beta_s**2/gamma_s**2",
    sp.simplify(bpp + 3 * bs2 / gs**2),
)

# прямое дифференцирование точной beta/beta_s
b1 = sp.simplify(sp.diff(beta_over_betas, delta).subs(delta, 0))
b2 = sp.simplify(sp.diff(beta_over_betas, delta, 2).subs(delta, 0))
check("direct (beta'/beta_s)_s = 1/gamma_s**2", sp.simplify(b1 - 1 / gs**2))
check("direct (beta''/beta_s)_s = -3 beta_s**2/gamma_s**2", sp.simplify(b2 + 3 * bs2 / gs**2))

beta_series = (
    1
    + delta / gs**2
    - three_half * bs2 / gs**2 * delta**2
)
check(
    "series beta/beta_s through delta**2",
    coeff(beta_over_betas - beta_series, 0)
    + coeff(beta_over_betas - beta_series, 1) * delta
    + coeff(beta_over_betas - beta_series, 2) * delta**2,
)

a = 1 / gs**2
b = -three_half * bs2 / gs**2
# 1/(1+a*delta+b*delta**2) = 1 - a delta + (a**2 - b) delta**2 + O(delta**3)
vsv_from_inversion = 1 - a * delta + (a**2 - b) * delta**2
vsv_expected = 1 - delta / gs**2 + (1 / gs**4 + three_half * bs2 / gs**2) * delta**2
check("inversion: a**2 = 1/gamma_s**4", sp.simplify(a**2 - 1 / gs**4))
check(
    "inversion: -b = 3 beta_s**2/(2 gamma_s**2)",
    sp.simplify(-b - three_half * bs2 / gs**2),
)
check(
    "v_s/v from inversion equals expected",
    sp.simplify(sp.expand(vsv_from_inversion - vsv_expected)),
)
check(
    "series v_s/v through delta**2",
    sp.simplify(coeff(vs_over_v - vsv_expected, 0))
    + sp.simplify(coeff(vs_over_v - vsv_expected, 1))
    + sp.simplify(coeff(vs_over_v - vsv_expected, 2)),
)
# equivalent writing of the delta**2 coefficient
check(
    "1/gamma**2 + beta**2/(2 gamma**2) = 1/gamma**4 + 3 beta**2/(2 gamma**2)",
    sp.simplify(
        (1 / gs**2 + bs2 / (2 * gs**2))
        - (1 / gs**4 + three_half * bs2 / gs**2)
    ),
)

eta0 = a0 - 1 / gs**2
eta1 = a1 - a0 / gs**2 + 1 / gs**4 + three_half * bs2 / gs**2
prod = sp.expand((1 + a0 * delta + a1 * delta**2) * vsv_expected)
check("product constant term", coeff(prod, 0) - 1)
check("product -> eta0", sp.simplify(coeff(prod, 1) - eta0))
check("product -> eta1", sp.simplify(coeff(prod, 2) - eta1))

# exact product of exact v_s/v and the length series, through delta**2
length = 1 + a0 * delta + a1 * delta**2
exact_tau = length * vs_over_v - 1
check("exact tau delta^1 = eta0", sp.simplify(coeff(exact_tau, 1) - eta0))
check("exact tau delta^2 = eta1", sp.simplify(coeff(exact_tau, 2) - eta1))

# (T-T_s)/T = tau/(1+tau) carries -eta0**2 at delta**2
tau_poly = eta0 * delta + eta1 * delta**2
over_own = sp.series(tau_poly / (1 + tau_poly), delta, 0, 3).removeO()
check(
    "(T-T_s)/T delta**2 = eta1 - eta0**2",
    sp.simplify(coeff(over_own, 2) - (eta1 - eta0**2)),
)

# --- Senichev (5): Delta v / v_s = delta/gamma**2 only ---
dC = a0 * delta + a1 * delta**2
dv = delta / gs**2
sen5 = sp.expand(dC - dv - dC * dv + dv**2)
kappa_sen = a1 - a0 / gs**2 + 1 / gs**4
check("Senichev (5) delta^1 = eta0", sp.simplify(coeff(sen5, 1) - eta0))
check("Senichev (5) delta^2 = kappa", sp.simplify(coeff(sen5, 2) - kappa_sen))
check(
    "lost term in eta1 is 3 beta_s**2/(2 gamma_s**2)",
    sp.simplify(eta1 - kappa_sen - three_half * bs2 / gs**2),
)
# restoring the second-order piece of Delta v does not change (Delta v)**2 at delta**2
dv_full = a * delta + b * delta**2  # (v-v_s)/v_s through delta**2
sen_full = sp.expand(dC - dv_full - dC * dv_full + dv_full**2)
check(
    "full Delta v restores eta1 at delta**2",
    sp.simplify(coeff(sen_full, 2) - eta1),
)
check(
    "(Delta v)**2 differs from (delta/gamma**2)**2 only at delta**3",
    sp.simplify(coeff(dv_full**2 - dv**2, 2)),
)

# prefactor of printed (10)
prefactor = gs**2 / (gs**2 * a0 - 1)
check("gamma**2/(gamma**2 alpha0 - 1) = 1/eta0", sp.simplify(prefactor - 1 / eta0))

# --- stability, both signs of eta0 ---
def stability_case(sign: int) -> None:
    eta0_c = sign * eta_abs
    cos_phis = -sign  # -sgn(eta0)
    ws2 = K * wrf * eta_abs
    tag = "eta0>0" if sign > 0 else "eta0<0"
    check(
        f"{tag}: cos phi_s * eta0 = -|eta0|",
        sp.simplify(cos_phis * eta0_c + eta_abs),
    )
    check(
        f"{tag}: omega_s**2 = K omega_rf |eta0|",
        sp.simplify(ws2 - K * wrf * sp.Abs(eta0_c)),
    )
    check(
        f"{tag}: K omega_rf cos phi_s = -omega_s**2/eta0",
        sp.simplify(K * wrf * cos_phis + ws2 / eta0_c),
    )
    # linear coefficient of delta_ddot is -omega_s**2
    check(
        f"{tag}: linear delta_ddot = -omega_s**2 delta",
        sp.simplify(K * wrf * cos_phis * eta0_c + ws2),
    )


stability_case(+1)
stability_case(-1)

# Senichev at phi=0: delta_ddot = -K*omega*eta0*delta, lambda**2 = -K*omega*eta0
# unstable when eta0<0 (lambda**2>0)
eta0_neg = -eta_abs
lam2_sen0 = sp.simplify(-K * wrf * eta0_neg)
check("Senichev phi=0, eta0<0: lambda**2 = +K omega |eta0| > 0", sp.simplify(lam2_sen0 - K * wrf * eta_abs))

# Senichev at phi=pi, eta0<0:
# delta_ddot = K*omega*tau, omega_s**2 = -K*omega*eta0
# average: <delta> = -(kappa<delta**2>+Lam)/eta0
ws2_pi = sp.simplify(-K * wrf * eta0_neg)
check("Senichev phi=pi, eta0<0: omega_s**2 > 0 equals K omega |eta0|", sp.simplify(ws2_pi - K * wrf * eta_abs))
# delta_ddot + omega_s**2 delta = K*omega*(kappa*delta**2+Lam)
# <delta> = K*omega/omega_s**2 * ( ... ) = -1/eta0 * (...)
avg_pi = sp.simplify(K * wrf / ws2_pi)
check(
    "Senichev stable point average factor = -1/eta0",
    sp.simplify(avg_pi + 1 / eta0_neg),
)

# printed (8) about phi=0, even though unstable for eta0<0:
# delta_ddot + (K*omega*eta0) delta = -K*omega*(kappa*delta**2+Lam)
# <delta> = -(...)/eta0
check(
    "average of printed (8) factor = -1/eta0",
    sp.simplify((-K * wrf) / (K * wrf * eta0_neg) + 1 / eta0_neg),
)

# phase-sign independence: phi_dot = s * omega * tau, <phi_dot>=0 => <tau>=0
s = sp.symbols("s", nonzero=True)
check(
    "either phase sign: <phi_dot>=0 implies <tau>=0",
    sp.simplify(sp.solve(s * wrf * sp.symbols("tau_mean"), sp.symbols("tau_mean"))[0]),
)

# <delta**2> = dm**2/2  =>  delta_eq = -eta1*dm**2/(2 eta0)
eta0_s, eta1_s = sp.symbols("eta_0 eta_1", nonzero=True)
d2 = dm**2 / 2
deq = sp.simplify(-eta1_s / eta0_s * d2)
check(
    "delta_eq = -eta1 delta_m**2 / (2 eta0)",
    sp.simplify(deq + eta1_s * dm**2 / (2 * eta0_s)),
)

# --- harmonic balance of the oscillator (draft; not used in the section) ---
# delta = b + a cos + c cos2
# delta_ddot + ws**2 delta = -(ws**2/eta0)*eta1*delta**2
a_amp, ws, e0, e1 = sp.symbols("a omega_s eta_0 eta_1", real=True)
b_const = -e1 * a_amp**2 / (2 * e0)
c_cos2 = e1 * a_amp**2 / (6 * e0)
# constant balance: ws**2 * b = -(ws**2/e0)*e1*(a**2/2)
rhs_const = -(ws**2 / e0) * e1 * (a_amp**2 / 2)
check("balance constant -> b", sp.simplify(ws**2 * b_const - rhs_const))
# cos2: c*(-4 ws**2 + ws**2) = -(ws**2/e0)*e1*(a**2/2)
rhs_cos2 = -(ws**2 / e0) * e1 * (a_amp**2 / 2)
check(
    "balance cos2 -> c",
    sp.simplify((-3 * ws**2) * c_cos2 - rhs_cos2),
)
# cos1 at order a**3: (ws**2 - omega**2)*a = -(ws**2/e0)*e1*(2*b*a + a*c)
# sigma = omega**2 - ws**2 = (ws**2/e0)*e1*(2*b+c)
two_b_plus_c = sp.simplify(2 * b_const + c_cos2)
sigma = sp.simplify((ws**2 / e0) * e1 * two_b_plus_c)
check(
    "omega**2 - omega_s**2 = -(5/6) omega_s**2 (eta1 a / eta0)**2",
    sp.simplify(sigma + sp.Rational(5, 6) * ws**2 * (e1**2 * a_amp**2) / e0**2),
)
# omega/omega_s = 1 + sigma/(2 omega_s**2) + ...
dw_over_w = sp.simplify(sigma / (2 * ws**2))
check(
    "Delta omega/omega_s = -(5/12) (eta1 delta_m / eta0)**2",
    sp.simplify(dw_over_w + sp.Rational(5, 12) * (e1**2 * a_amp**2) / e0**2),
)
# split of 5/12: constant gives -1/2, second harmonic gives +1/12
part_b = sp.simplify((e1 / e0) * (2 * b_const) / a_amp**0 * a_amp**0)
# (eta1/eta0)*(2b) / a**2 * a**2, then half for the frequency
contrib_b = sp.simplify((e1 / e0) * (2 * b_const) / (2 * a_amp**2) * a_amp**2 / 1)
# clearer:
# sigma/ws**2 = (e1/e0)*(2b+c)
# Delta w/w = half of that
piece_b = sp.simplify(half * (e1 / e0) * (2 * b_const) / a_amp**2)
piece_c = sp.simplify(half * (e1 / e0) * c_cos2 / a_amp**2)
check("frequency piece from constant = -1/2 * (eta1/eta0)**2", sp.simplify(piece_b + half * e1**2 / e0**2))
check("frequency piece from cos2 = +1/12 * (eta1/eta0)**2", sp.simplify(piece_c - sp.Rational(1, 12) * e1**2 / e0**2))
check(
    "5/12 = 1/2 - 1/12",
    sp.simplify(sp.Rational(5, 12) - (half - sp.Rational(1, 12))),
)

# phi**2 correction at eta1=0: Delta omega/omega_s = -phi_m**2/16, no constant
psi = sp.symbols("psi", real=True)
phi_m, delta_m = sp.symbols("phi_m delta_m", real=True)
phi = phi_m * sp.sin(psi)
delta_lin = delta_m * sp.cos(psi)
rhs_phi = sp.trigsimp(sp.expand(phi**2 / 2 * delta_lin))
expected_phi = phi_m**2 * delta_m / 8 * (sp.cos(psi) - sp.cos(3 * psi))
check("phi**2/2 * delta fourier", sp.simplify(sp.expand(rhs_phi - expected_phi)))
# no constant term
const_phi = sp.integrate(rhs_phi, (psi, 0, 2 * sp.pi)) / (2 * sp.pi)
check("phi**2 correction has no constant", sp.simplify(const_phi))
# (-Omega**2 + ws**2) delta_m = ws**2 * phi_m**2 * delta_m / 8
# Omega = ws + dw, -2 ws dw = ws**2 phi_m**2 / 8
check(
    "Delta omega/omega_s = -phi_m**2/16",
    sp.simplify(-phi_m**2 / 16 - (-(phi_m**2 / 8) / 2)),
)

# numeric spot-check, gamma_s=2, not a ring
num = {gs: 2}
exact_num = vs_over_v.series(delta, 0, 3).removeO().subs(num)
expect_num = vsv_expected.subs(num)
check("numeric gamma_s=2 series", sp.simplify(sp.expand(exact_num - expect_num)))

print()
n_fail = sum(1 for _, ok, _ in results if not ok)
n_pass = sum(1 for _, ok, _ in results if ok)
print(f"{n_pass} passed, {n_fail} failed, {len(results)} checks")
sys.exit(1 if n_fail else 0)
