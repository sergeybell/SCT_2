# Проверки главы 6. Только рациональные числа и короткие ряды.
from fractions import Fraction
import sympy as sp

d, g, n = sp.symbols("delta gamma2 n", rational=True)
# g = gamma_s**2


def check(name, got, expected):
    diff = sp.simplify(sp.together(got - expected))
    ok = diff == 0
    print(f"{'OK' if ok else 'FAIL'}  {name}")
    if not ok:
        print("   diff =", diff)
    return ok


all_ok = True

# --- 1. Постоянное решение уравнения для D против (eq:D-deflector) ---
# D'' + (2-n+1/g)/rho**2 * D = (1+1/g)/rho, D=const
# D/rho = (1+1/g) / (2-n+1/g)
nu2 = 2 - n + 1 / g
d_over_rho = (1 + 1 / g) / nu2
all_ok &= check(
    "D/rho = (1+1/g)/nu_x^2",
    d_over_rho,
    (1 + 1 / g) / (2 - n + 1 / g),
)

# Таблица главы 4: линейная дисперсия кольца из одного дефлектора
rows = [
    (1, Fraction(25, 16), Fraction(1)),
    (1, Fraction(4), Fraction(1)),
    (2, Fraction(25, 16), Fraction(41, 16)),
    (2, Fraction(4), Fraction(5)),
    (Fraction(3, 2), Fraction(25, 16), Fraction(82, 57)),
]
for n_v, g_v, expect in rows:
    got = sp.simplify(d_over_rho.subs({n: n_v, g: g_v}))
    all_ok &= check(f"D/rho n={n_v} g={g_v}", got, expect)

# --- 2. Независимый линейный баланс pv = k (rho/r)^{n-1} ---
# gamma = gamma_out - beta**2 * gamma * psi, psi = x/rho на линейном порядке
# (gamma - 1/gamma) / (gamma_s - 1/gamma_s) = (rho/r)^{n-1}
# На линейном порядке это (1-d)*(1+1/g) = (1-n)*d, откуда d=(1+1/g)/(2-n+1/g)


def balance_linear(n_v, g_v):
    dd = sp.symbols("dd")
    # левая часть: (1 - D/rho) * (1 + 1/g), правая: (1-n)*(D/rho)
    eq = (1 - dd) * (1 + 1 / g_v) - (1 - n_v) * dd
    sol = sp.solve(sp.Eq(eq, 0), dd)[0]
    return sp.simplify(sol)


for n_v, g_v, expect in rows:
    all_ok &= check(
        f"баланс pv n={n_v} g={g_v}",
        balance_linear(n_v, g_v),
        expect,
    )

# --- 3. Тождество «магнитная форма + поправка» ---
a0, a1, d1, d2, p2 = sp.symbols("alpha0 alpha1 d1 d2 p2")
eta0_m_form = a0 - 1 / g
eta0_e = (1 + 1 / g) * a0 - 1 / g
all_ok &= check("eta0^E - magnetic form", eta0_e - eta0_m_form, a0 / g)

eta1_e = (
    p2 / 2
    + (1 + 1 / g) * d2
    + (5 - n) / (2 * g) * d1**2
    - 4 * a0 / g
    + (3 * g - 1) / (2 * g**2)
)
# геометрические alpha на тех же D, D2
# a0 и d1 не отождествляем символьно: в кольце a0 = <d1>, а <d1**2> не есть a0**2
eta1_m_form = a1 - a0 / g + (3 * g - 1) / (2 * g**2)
# при постоянных d1, d2 и D'=0: a0=d1, a1=d2
sub_flat = {a0: d1, a1: d2, p2: 0}
corr = sp.simplify((eta1_e - eta1_m_form).subs(sub_flat))
expect_corr = sp.simplify((1 / g) * (d2 + (5 - n) / 2 * d1**2 - 3 * d1))
all_ok &= check("eta1 correction, flat orbit", corr, expect_corr)

# общий случай: a1 = <d2> + p2/2, поправка не содержит p2
# здесь a0, d1, d2 --- уже средние, и член d1**2 означает <(D/rho)**2>
a0s, d2s, d1sq = sp.symbols("a0 d2avg d1sq")
eta1_e_avg = (
    p2 / 2
    + (1 + 1 / g) * d2s
    + (5 - n) / (2 * g) * d1sq
    - 4 * a0s / g
    + (3 * g - 1) / (2 * g**2)
)
eta1_m_avg = (d2s + p2 / 2) - a0s / g + (3 * g - 1) / (2 * g**2)
corr_avg = sp.factor(sp.together(eta1_e_avg - eta1_m_avg))
expect_avg = (1 / g) * (d2s + (5 - n) / 2 * d1sq - 3 * a0s)
all_ok &= check("eta1 correction, averages", corr_avg, expect_avg)

# --- 4. Круговая орбита n=2, gamma_s=2: ряд tau до delta**2 ---
# rho=1, g=4, beta2=1-1/g=3/4, k=beta2*gamma = (3/4)*2 = 3/2
# psi = 1 - 1/r  (n=2)
# gamma_out**2 = 1 + (g-1)*(1+delta)**2
# gamma_loc = gamma_out - k*psi
# pv = gamma_loc - 1/gamma_loc = k / r
gv = sp.Integer(4)
nv = sp.Integer(2)
beta2 = 1 - 1 / gv
gam_s = sp.Integer(2)
k = beta2 * gam_s  # 3/2
rr = sp.symbols("r", positive=True)
delta = sp.symbols("delta")
gam_out = sp.sqrt(1 + (gv - 1) * (1 + delta) ** 2)
psi = 1 - 1 / rr  # n=2, rho=1
gam_loc = gam_out - k * psi
pv = gam_loc - 1 / gam_loc
# сила: pv = k * r**(1-n) = k/r
force = sp.series(pv - k / rr, delta, 0, 3).removeO()
# r = 1 + a*delta + b*delta**2
a, b = sp.symbols("a b")
force_ab = force.subs(rr, 1 + a * delta + b * delta**2)
force_ab = sp.series(force_ab, delta, 0, 3).removeO()
sol = sp.solve(
    [force_ab.coeff(delta, 1), force_ab.coeff(delta, 2)],
    [a, b],
    dict=True,
)[0]
a_sol = sp.simplify(sol[a])
b_sol = sp.simplify(sol[b])
all_ok &= check("orbit a = D/rho", a_sol, Fraction(5))
all_ok &= check("orbit b = D2/rho", b_sol, Fraction(109, 8))

# скорость на окружности и tau = r * v_s/v - 1
# beta_loc**2 = 1 - 1/gam_loc**2, v_s/v = beta_s/beta_loc
beta_s = sp.sqrt(beta2)
r_of_d = 1 + a_sol * delta + b_sol * delta**2
gam_l = sp.series(gam_loc.subs(rr, r_of_d), delta, 0, 3).removeO()
v_ratio = sp.series(beta_s / sp.sqrt(1 - 1 / gam_l**2), delta, 0, 3).removeO()
tau = sp.series(r_of_d * v_ratio - 1, delta, 0, 3).removeO()
eta0 = sp.simplify(tau.coeff(delta, 1))
eta1 = sp.simplify(tau.coeff(delta, 2))
all_ok &= check("tau eta0", eta0, 6)
all_ok &= check("tau eta1", eta1, Fraction(87, 4))

# общая формула главы 4 на этих D, D2
a0_v = a_sol
d2_v = b_sol
d1_v = a_sol
eta0_f = (1 + 1 / gv) * a0_v - 1 / gv
eta1_f = (
    (1 + 1 / gv) * d2_v
    + (5 - nv) / (2 * gv) * d1_v**2
    - 4 * a0_v / gv
    + (3 * gv - 1) / (2 * gv**2)
)
all_ok &= check("formula eta0", sp.simplify(eta0_f), eta0)
all_ok &= check("formula eta1", sp.simplify(eta1_f), eta1)

# магнитная форма на тех же alpha
eta0_mf = a0_v - 1 / gv
eta1_mf = d2_v - a0_v / gv + (3 * gv - 1) / (2 * gv**2)
all_ok &= check("delta eta0", sp.simplify(eta0 - eta0_mf), a0_v / gv)
all_ok &= check(
    "delta eta1",
    sp.simplify(eta1 - eta1_mf),
    (1 / gv) * (d2_v + (5 - nv) / 2 * d1_v**2 - 3 * a0_v),
)

# --- 5. Гладкое приближение сводится к (eq:Lambda-ring) ---
# C=2*pi*rho, beta=rho/nu, gamma_twiss = nu/rho
Jx, Jy = sp.symbols("Jx Jy", positive=True)
nu_x, nu_y, rho, gs, xs = sp.symbols("nux nuy rho gamma2 xnl", positive=True)
nn = sp.symbols("n")
Cs = 2 * sp.pi * rho
beta_x = rho / nu_x
beta_y = rho / nu_y
# интеграл (eq:Lambda-es-J) при постоянных функциях
Lam = (1 / Cs) * (
    (Jx / 2) * (nu_x / rho + (5 - nn) * beta_x / (gs * rho**2)) * Cs
    + (Jy / 2) * (nu_y / rho + (nn - 1) * beta_y / (gs * rho**2)) * Cs
    + (1 + 1 / gs) * xs / rho * Cs
)
Lam = sp.simplify(Lam)
Lam_ring = (
    Jx / (2 * rho) * (nu_x + (5 - nn) / (gs * nu_x))
    + Jy / (2 * rho) * (nu_y + (nn - 1) / (gs * nu_y))
    + (1 + 1 / gs) * xs / rho
)
all_ok &= check("smooth = Lambda-ring", Lam, sp.simplify(Lam_ring))

# вертикальная скобка: (n-1)/nu_y = nu_y при nu_y**2 = n-1
Lam_y = sp.simplify(
    (Jy / (2 * rho) * (nu_y + (nn - 1) / (gs * nu_y))).subs(nn, nu_y**2 + 1)
)
all_ok &= check(
    "vertical bracket",
    Lam_y,
    sp.simplify(Jy * nu_y * (1 + 1 / gs) / (2 * rho)),
)

# --- 6. Отношение отброшенного члена gamma к члену с delta_eq ---
# |eta0| / (2 |eta1| g) для n=1: eta0=1, eta1=(g+1)/(2g)
eta0_n1 = 1
eta1_n1 = (g + 1) / (2 * g)
ratio_n1 = sp.simplify(eta0_n1 / (2 * eta1_n1 * g))
all_ok &= check("ratio n=1", ratio_n1, 1 / (g + 1))
print("ratio at gamma_s=2:", ratio_n1.subs(g, 4))

# --- 7. Локальный множитель (1 - D/rho) при n=1 ---
all_ok &= check("n=1 cancels local linear term", 1 - d_over_rho.subs(n, 1), 0)

print("ALL OK" if all_ok else "SOME FAILED")
