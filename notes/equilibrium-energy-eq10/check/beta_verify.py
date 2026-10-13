"""Проверка выкладок раздела о бетатронном удлинении орбиты.

Каждая проверка — тождество sympy (или численное следствие того же
разложения). Код возврата 0, только если все проверки прошли.
"""

from __future__ import annotations

import sys

import sympy as sp

PASS = 0
FAIL = 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global PASS, FAIL
    if ok:
        PASS += 1
        print(f"PASS  {name}")
    else:
        FAIL += 1
        print(f"FAIL  {name}")
        if detail:
            print(f"      {detail}")


def series_coeff(expr, var, n):
    return sp.simplify(expr.series(var, 0, n + 1).removeO().coeff(var, n))


def main() -> int:
    # ------------------------------------------------------------------
    # 1. Элемент пути: (x/rho)^2 сокращается, остаётся 1 + x/rho + угол/2
    # ------------------------------------------------------------------
    eps, h, xp, yp = sp.symbols("varepsilon h x_p y_p", real=True)
    dl = sp.sqrt((1 + eps * h) ** 2 + (eps * xp) ** 2 + (eps * yp) ** 2)
    dl_series = dl.series(eps, 0, 3).removeO()
    target = 1 + eps * h + eps**2 * (xp**2 + yp**2) / 2
    check(
        "dl to 2nd order: 1 + x/rho + (x'^2+y'^2)/2",
        sp.expand(dl_series - target) == 0,
        str(sp.expand(dl_series - target)),
    )

    # Явная отмена h^2: член порядка eps^2 не содержит h^2.
    h2 = sp.expand(dl_series).coeff(eps, 2)
    check(
        "coefficient of (x/rho)^2 in dl is 0",
        sp.simplify(sp.diff(h2, h, 2)) == 0 and sp.simplify(h2.coeff(h, 2)) == 0,
        str(sp.expand(h2)),
    )

    # Поэтому перекрёстное (D*delta)*(x_beta)/rho^2 в длину на этом порядке не входит.
    # Подстановка h = h_d + h_b и сборка коэффициента при h_d*h_b.
    hd, hb = sp.symbols("h_delta h_beta", real=True)
    h2_cross = sp.expand(h2.subs(h, hd + hb))
    check(
        "no h_delta*h_beta cross term in 2nd-order dl",
        sp.simplify(h2_cross.coeff(hd, 1).coeff(hb, 1)) == 0,
        str(h2_cross),
    )

    # ------------------------------------------------------------------
    # 2. <x'^2> = J * gamma при alpha != 0
    # ------------------------------------------------------------------
    s = sp.symbols("s", real=True)
    J = sp.symbols("J", positive=True)
    beta = sp.Function("beta", real=True, positive=True)
    psi = sp.Function("psi", real=True)
    alpha = sp.symbols("alpha", real=True)

    x = sp.sqrt(2 * J * beta(s)) * sp.cos(psi(s))
    xp_expr = sp.diff(x, s)
    xp_expr = xp_expr.subs(sp.Derivative(beta(s), s), -2 * alpha)
    xp_expr = xp_expr.subs(sp.Derivative(psi(s), s), 1 / beta(s))
    xp_expr = sp.simplify(xp_expr)

    theta = sp.symbols("theta", real=True)
    xp_th = sp.simplify(xp_expr.subs(psi(s), theta))
    expected_xp = -sp.sqrt(2 * J / beta(s)) * (alpha * sp.cos(theta) + sp.sin(theta))
    check(
        "x' = -sqrt(2J/beta)*(alpha*cos + sin)",
        sp.simplify(xp_th - expected_xp) == 0,
        str(sp.simplify(xp_th - expected_xp)),
    )

    avg_xp2 = sp.integrate(sp.simplify(xp_th**2), (theta, 0, 2 * sp.pi)) / (2 * sp.pi)
    avg_xp2 = sp.simplify(avg_xp2)
    gamma = (1 + alpha**2) / beta(s)
    check(
        "<x'^2> = J*gamma with alpha kept",
        sp.simplify(avg_xp2 - J * gamma) == 0,
        str(avg_xp2),
    )

    # Средние линейных членов по фазе — нули (обнуление x_beta/rho и перекрёстного).
    avg_x = sp.integrate(sp.cos(theta), (theta, 0, 2 * sp.pi)) / (2 * sp.pi)
    avg_xp = sp.integrate(xp_th, (theta, 0, 2 * sp.pi)) / (2 * sp.pi)
    check("<x_beta> phase average is 0", sp.simplify(avg_x) == 0)
    check("<x'_beta> phase average is 0", sp.simplify(avg_xp) == 0)

    # Числовая точка с alpha != 0, чтобы тождество не выродилось.
    num = sp.simplify(avg_xp2).subs({J: sp.Rational(3, 7), alpha: sp.Rational(3, 2), beta(s): sp.Rational(5, 2)})
    gamma_num = (1 + sp.Rational(3, 2) ** 2) / sp.Rational(5, 2) * sp.Rational(3, 7)
    check(
        "numeric <x'^2> at alpha=3/2, beta=5/2",
        sp.simplify(num - gamma_num) == 0,
        f"{num} vs {gamma_num}",
    )

    # ------------------------------------------------------------------
    # 3. Скорость и коэффициент 1 перед Delta L / C_s
    # ------------------------------------------------------------------
    gs, delta = sp.symbols("gamma_s delta", positive=True)
    # gamma^2 = 1 + (beta*gamma)^2, beta*gamma = beta_s*gamma_s*(1+delta),
    # beta_s^2 * gamma_s^2 = gamma_s^2 - 1
    g2 = 1 + (gs**2 - 1) * (1 + delta) ** 2
    vs_over_v = sp.sqrt(g2) / (gs * (1 + delta))
    a1 = series_coeff(vs_over_v, delta, 1)
    a2 = series_coeff(vs_over_v, delta, 2)
    check("v_s/v linear coeff = -1/gamma^2", sp.simplify(a1 + 1 / gs**2) == 0, str(a1))

    beta2 = 1 - 1 / gs**2
    a2_expected = 1 / gs**4 + 3 * beta2 / (2 * gs**2)
    check(
        "v_s/v quadratic coeff = 1/gamma^4 + 3 beta^2/(2 gamma^2)",
        sp.simplify(a2 - a2_expected) == 0,
        str(sp.simplify(a2 - a2_expected)),
    )

    alpha0, alpha1, lam = sp.symbols("alpha_0 alpha_1 Lambda", real=True)
    vs_poly = 1 + a1 * delta + a2 * delta**2
    # lam не зависит от delta в удерживаемом порядке
    tau = sp.expand(vs_poly * (1 + alpha0 * delta + alpha1 * delta**2 + lam) - 1)
    eta0 = alpha0 - 1 / gs**2
    eta1 = alpha1 - alpha0 / gs**2 + 1 / gs**4 + 3 * beta2 / (2 * gs**2)
    kept = eta0 * delta + eta1 * delta**2 + lam
    rem = sp.expand(tau - kept)
    # Произведение обрезанных рядов содержит delta**3 и delta**4. До delta**2 остаток
    # при Lambda=0 обязан быть нулём: это и есть eta0, eta1.
    rem_disp = sp.series(rem.subs(lam, 0), delta, 0, 3).removeO()
    check(
        "dispersion part reproduces eta0, eta1 through delta**2",
        sp.simplify(rem_disp) == 0,
        str(sp.simplify(rem_disp)),
    )
    # Коэффициент при Lambda^1, delta^0 равен 0 в остатке? kept уже содержит lam,
    # так что постоянный по delta кусок остатка должен быть 0.
    rem_at_0 = sp.simplify(sp.series(rem, delta, 0, 1).removeO())
    check(
        "constant term of tau is exactly Lambda (coeff 1)",
        rem_at_0 == 0,
        str(rem_at_0),
    )
    # Линейный по delta кусок остатка равен -Lambda/gamma^2 * delta, без alpha0*Lambda.
    rem_s = sp.series(rem, delta, 0, 3).removeO()
    lin = series_coeff(rem_s, delta, 1)
    check(
        "O(delta*Lambda) remainder is -Lambda*delta/gamma^2",
        sp.simplify(lin + lam / gs**2) == 0,
        str(sp.simplify(lin)),
    )
    # alpha0 не остаётся в этой поправке: eta0 - alpha0 = -1/gamma^2.
    # Перенормировка на C(delta): (v_s/v)*(C/C_s)*lam_d, lam_d = lam/(C/C_s).
    C_over = 1 + alpha0 * delta + alpha1 * delta**2
    lam_d = lam / C_over  # Delta L / C(delta), ряд по delta
    tau_d = sp.series(vs_poly * C_over * (1 + lam_d) - 1, delta, 0, 2).removeO()
    # до первого порядка по delta: eta0*delta + lam - lam*delta/gamma^2
    tau_d_expected = eta0 * delta + lam - lam * delta / gs**2
    # vs_poly * C_over * (1 + lam_d) = vs_poly * C_over + vs_poly * lam
    # потому что C_over * lam_d = lam. Совпадает с нормировкой на C_s.
    check(
        "normalising Delta L by C(delta) gives the same -delta*Lambda/gamma^2",
        sp.simplify(sp.series(tau_d - tau_d_expected, delta, 0, 2).removeO()) == 0,
        str(sp.series(tau_d - (eta0 * delta + lam * (1 - delta / gs**2)), delta, 0, 2)),
    )

    # ------------------------------------------------------------------
    # 4. Гладкий предел и множитель pi/(2L)
    # ------------------------------------------------------------------
    Jx, Jy, nux, nuy, C, L = sp.symbols("J_x J_y nu_x nu_y C_s L", positive=True)
    # alpha = 0 => int gamma ds = int ds/beta = 2*pi*nu
    Lam_ang = (Jx / 2 * 2 * sp.pi * nux + Jy / 2 * 2 * sp.pi * nuy) / C
    Lam_smooth = sp.pi * (Jx * nux + Jy * nuy) / C
    check(
        "smooth Lambda = pi (J_x nu_x + J_y nu_y) / C_s",
        sp.simplify(Lam_ang - Lam_smooth) == 0,
    )

    ex, ey = sp.symbols("varepsilon_x varepsilon_y", positive=True)
    eps = sp.symbols("varepsilon", positive=True)
    avg_alpha0 = sp.simplify(avg_xp2.subs(alpha, 0))
    check(
        "<x'^2> = J/beta when alpha=0",
        sp.simplify(avg_alpha0 - J / beta(s)) == 0,
        str(avg_alpha0),
    )
    # Senichev: <x'^2> = eps/(2 beta). Вместе с J/beta это eps = 2J.
    eps_sol = sp.solve(sp.Eq(J / beta(s), eps / (2 * beta(s))), eps)[0]
    check(
        "epsilon = 2J from <x'^2> = J/beta = epsilon/(2 beta)",
        sp.simplify(eps_sol - 2 * J) == 0,
        str(eps_sol),
    )
    Lam_eps = sp.pi * (ex * nux + ey * nuy) / (2 * C)
    check(
        "J=epsilon/2 turns smooth Lambda into Senichev (11)",
        sp.simplify(Lam_eps.subs({ex: 2 * Jx, ey: 2 * Jy}) - Lam_smooth) == 0,
    )

    # Промежуточная цепочка (11): (1/L) * int (1/2) * (eps/(2 beta)) ds,
    # int ds/beta = 2*pi*nu.
    eps, nu = sp.symbols("varepsilon nu", positive=True)
    chain = (1 / L) * (sp.Rational(1, 2)) * (eps / 2) * (2 * sp.pi * nu)
    printed = sp.pi * eps * nu / (2 * L)
    check(
        "two factors 1/2 give pi*epsilon*nu/(2L)",
        sp.simplify(chain - printed) == 0,
        str(sp.simplify(chain - printed)),
    )

    # Отброшенный при замене gamma -> 1/beta вклад: (J/2) int alpha^2/beta ds.
    # int gamma = int (1+alpha^2)/beta = 2*pi*nu + int alpha^2/beta, если int ds/beta = 2*pi*nu.
    al = sp.symbols("alpha")
    check(
        "gamma - 1/beta = alpha^2/beta",
        sp.simplify((1 + al**2) / beta(s) - 1 / beta(s) - al**2 / beta(s)) == 0,
    )

    # ------------------------------------------------------------------
    # 5. delta_eq из <tau> = 0 и перекрёстный сдвиг частоты
    # ------------------------------------------------------------------
    dm, eta0_s, eta1_s, Lam_s = sp.symbols("delta_m eta_0 eta_1 Lambda", real=True)
    # <delta^2> = delta_m**2 / 2 на удерживаемом порядке
    deq = -(eta1_s * dm**2 / 2 + Lam_s) / eta0_s
    avg_tau = sp.simplify(eta0_s * deq + eta1_s * dm**2 / 2 + Lam_s)
    check("<tau>=0 fixes delta_eq", avg_tau == 0, str(avg_tau))
    prefactor = gs**2 / (gs**2 * alpha0 - 1)
    check(
        "printed prefactor gamma^2/(gamma^2*alpha0-1) equals 1/eta0",
        sp.simplify(prefactor - 1 / (alpha0 - 1 / gs**2)) == 0,
        str(sp.simplify(prefactor - 1 / eta0)),
    )

    # Гармонический баланс:
    # ddot + ws^2 delta = -(ws^2/eta0)*(eta1*delta**2 + Lambda)
    # delta = d0 + a cos + b cos 2psi, psi = omega t
    ws, a = sp.symbols("omega_s a", positive=True)
    d0, b, omega = sp.symbols("d0 b omega", real=True)
    # Коэффициенты delta^2: const a^2/2, cos: a*(2*d0+b), cos2: a^2/2
    # (члены d0^2, b^2, d0*b — следующий порядок)
    c2 = -ws**2 * eta1_s / eta0_s
    c0 = -ws**2 * Lam_s / eta0_s
    eq_const = sp.Eq(ws**2 * d0, c0 + c2 * a**2 / 2)
    eq_cos2 = sp.Eq(b * (ws**2 - 4 * omega**2), c2 * a**2 / 2)
    eq_cos = sp.Eq(a * (ws**2 - omega**2), c2 * a * (2 * d0 + b))

    d0_sol = sp.solve(eq_const, d0)[0]
    # omega^2 = ws^2 + O(2), поэтому в уравнении на b заменяем omega -> ws
    b_sol = sp.solve(eq_cos2.subs(omega, ws), b)[0]
    # ws^2 - omega^2 = c2 * (2*d0+b)
    omega2 = sp.simplify(ws**2 - c2 * (2 * d0_sol + b_sol))
    omega2_pred = ws**2 * (
        1
        - 2 * eta1_s * Lam_s / eta0_s**2
        - sp.Rational(5, 6) * eta1_s**2 * a**2 / eta0_s**2
    )
    check(
        "harmonic balance: omega^2",
        sp.simplify(omega2 - omega2_pred) == 0,
        str(sp.simplify(omega2 - omega2_pred)),
    )
    # Линеаризация около центра (только сдвиг от Lambda, без второй гармоники):
    # omega^2/ws^2 = 1 + 2*(eta1/eta0)*d0, d0 = -Lam/eta0 при a=0
    # => omega/ws = 1 - eta1*Lam/eta0**2
    d0_lam = sp.simplify(d0_sol.subs(a, 0))
    check(
        "fixed point at a=0 is d0=-Lambda/eta0",
        sp.simplify(d0_lam + Lam_s / eta0_s) == 0,
        str(d0_lam),
    )
    slope = sp.simplify(1 + 2 * (eta1_s / eta0_s) * d0_lam)
    check(
        "linearisation cross shift omega^2/ws^2 = 1 - 2 eta1 Lambda/eta0^2",
        sp.simplify(slope - (1 - 2 * eta1_s * Lam_s / eta0_s**2)) == 0,
        str(slope),
    )
    # Вклад второй гармоники меняет 1/2 на 5/12 и не трогает член с Lambda.
    # Разность (omega2/ws**2 - 1) и линеаризации около полного d0:
    # линеаризация дала бы -eta1**2 a**2/eta0**2 - 2 eta1 Lam/eta0**2,
    # баланс даёт -5/6 и тот же -2 по Lambda.
    full_d0 = sp.simplify(d0_sol)
    lin_full = sp.simplify(1 + 2 * (eta1_s / eta0_s) * full_d0)
    bal = sp.simplify(omega2 / ws**2)
    delta_from_harmonic = sp.simplify(bal - lin_full)
    check(
        "second harmonic corrects only the eta1**2 piece, by +1/6 in omega^2",
        sp.simplify(delta_from_harmonic - sp.Rational(1, 6) * eta1_s**2 * a**2 / eta0_s**2) == 0,
        str(delta_from_harmonic),
    )
    # И сам b.
    check(
        "second harmonic b = eta1 a^2 / (6 eta0)",
        sp.simplify(b_sol - eta1_s * a**2 / (6 * eta0_s)) == 0,
        str(sp.simplify(b_sol)),
    )
    # omega/ws до второго порядка (Lambda считаем величиной порядка a**2:
    # подставляем Lambda = lam2 * a**2 и раскладываем).
    lam2 = sp.symbols("lam2", real=True)
    # sqrt(1+u) = 1 + u/2 + O(u**2). Символ eta0 вещественный и может быть
    # отрицательным, поэтому series(sqrt(...)) тащит Abs; раскладываем вручную.
    u = sp.simplify(omega2 / ws**2 - 1).subs(Lam_s, lam2 * a**2)
    ratio_small = sp.series(1 + u / 2, a, 0, 3).removeO()
    ratio_pred = (
        1
        - eta1_s * lam2 * a**2 / eta0_s**2
        - sp.Rational(5, 12) * eta1_s**2 * a**2 / eta0_s**2
    )
    check(
        "omega/ws = 1 - 5/12 eta1^2 a^2/eta0^2 - eta1 Lambda/eta0^2",
        sp.simplify(ratio_small - ratio_pred) == 0,
        str(sp.simplify(ratio_small - ratio_pred)),
    )

    # Численная частота малой амплитуды около точной неподвижной точки.
    # Уравнение: ddot = -ws^2 * delta - (ws^2/eta0)*(eta1*delta**2 + Lam)
    # Неподвижная точка: delta + (eta1/eta0)*delta**2 + Lam/eta0 = 0.
    def small_osc_ratio(eta0_v, eta1_v, lam_v):
        # квадратное уравнение (eta1/eta0) d^2 + d + lam/eta0 = 0
        aa = eta1_v / eta0_v
        bb = 1.0
        cc = lam_v / eta0_v
        disc = bb**2 - 4 * aa * cc
        r1 = (-bb + disc**0.5) / (2 * aa)
        r2 = (-bb - disc**0.5) / (2 * aa)
        target = -lam_v / eta0_v
        d_fix = r1 if abs(r1 - target) < abs(r2 - target) else r2
        w2 = 1.0 + 2.0 * (eta1_v / eta0_v) * d_fix
        return w2**0.5, 1.0 - eta1_v * lam_v / eta0_v**2

    # Следующий член ряда: omega/ws = 1 - eta1*Lam/eta0**2 - (3/2)*eta1**2*Lam**2/eta0**4 + ...
    eta0_v, eta1_v, lam_v = -0.5, 1.2, 1.0e-4
    ratio_num, ratio_lead = small_osc_ratio(eta0_v, eta1_v, lam_v)
    ratio_next = ratio_lead - 1.5 * eta1_v**2 * lam_v**2 / eta0_v**4
    check(
        "numeric small-oscillation cross shift matches through Lambda**2",
        abs(ratio_num - ratio_next) < 1e-8,
        f"exact {ratio_num}, through Lambda^2 {ratio_next}, diff {ratio_num - ratio_next}",
    )

    # Интегрирование ОДУ: коэффициент 5/12 и перекрёстный член.
    measure = _measure_frequency_shift()
    print(
        f"INFO  ODE shifts: eta1-only coeff {measure['coeff_5_12']:.6f}"
        f" (expect {-5/12:.6f}),"
        f" cross coeff {measure['coeff_cross']:.6f} (expect -1)"
    )
    check(
        "ODE: pure eta1 shift coefficient is -5/12",
        abs(measure["coeff_5_12"] + sp.Rational(5, 12)) < 2e-3,
        f"measured {measure['coeff_5_12']}, expected {-5/12}",
    )
    check(
        "ODE: cross coefficient of eta1*Lambda/eta0^2 is -1",
        abs(measure["coeff_cross"] + 1) < 2e-3,
        f"measured {measure['coeff_cross']}, expected -1",
    )

    print(f"\n{PASS} passed, {FAIL} failed")
    return 0 if FAIL == 0 else 1


def _measure_frequency_shift() -> dict:
    """RK4 по уравнению ddot + ws^2 delta = -(ws^2/eta0)*(eta1 delta^2 + Lambda)."""
    eta0 = -0.8
    eta1 = 0.5
    ws = 1.0
    a = 0.02

    def period(lam: float) -> float:
        d0 = -(eta1 * a * a / 2 + lam) / eta0
        b = eta1 * a * a / (6 * eta0)
        # начальная точка — приближённый максимум
        y0 = d0 + a + b
        v0 = 0.0
        dt = 2.0e-4

        def f(y, v):
            acc = -ws**2 * y - (ws**2 / eta0) * (eta1 * y * y + lam)
            return v, acc

        y, v = y0, v0
        t = 0.0
        # сразу уходим с максимума: скорость становится отрицательной
        seen_negative = False
        seen_positive = False
        t_max = None
        # запас по времени: два периода
        nsteps = int(4 * 3.141592653589793 / dt) + 10
        for _ in range(nsteps):
            k1y, k1v = f(y, v)
            k2y, k2v = f(y + 0.5 * dt * k1y, v + 0.5 * dt * k1v)
            k3y, k3v = f(y + 0.5 * dt * k2y, v + 0.5 * dt * k2v)
            k4y, k4v = f(y + dt * k3y, v + dt * k3v)
            y_new = y + (dt / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y)
            v_new = v + (dt / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
            t_new = t + dt
            if v < 0.0:
                seen_negative = True
            if seen_negative and v > 0.0:
                seen_positive = True
            # следующий максимум: скорость из + в -
            if seen_positive and v > 0.0 and v_new <= 0.0:
                frac = v / (v - v_new)
                t_max = t + frac * dt
                break
            y, v, t = y_new, v_new, t_new
        if t_max is None:
            raise RuntimeError("period not found")
        return t_max

    def shift(lam: float) -> float:
        t = period(lam)
        omega = 2.0 * 3.141592653589793 / t
        return omega / ws - 1.0

    s0 = shift(0.0)
    scale_a = eta1**2 * a**2 / eta0**2
    coeff_5_12 = s0 / scale_a

    lam = 1.0e-4
    s1 = shift(lam)
    # s1 ≈ -5/12 * scale_a - eta1*lam/eta0**2 + O(a^2 lam, lam^2)
    cross = s1 - s0
    coeff_cross = cross / (eta1 * lam / eta0**2)
    return {"coeff_5_12": coeff_5_12, "coeff_cross": coeff_cross, "shift0": s0, "shift1": s1}


if __name__ == "__main__":
    sys.exit(main())
