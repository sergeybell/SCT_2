# Fixed-number checks for chapter 5. One 360-degree deflector.
# gamma_s=2, rho=1. No symbolic series.

import math
import numpy as np
from scipy.integrate import quad, solve_ivp

RHO = 1.0
GS = 2.0
B2 = 1.0 - 1.0 / GS**2
B = math.sqrt(B2)
PS = GS * B
VS = B
PV = PS * VS


def psi_r(r, n):
    return (1.0 - (r / RHO) ** (1.0 - n)) / (n - 1.0)


def gkin(r, ginf, n):
    return ginf - B2 * GS * psi_r(r, n)


def ginf(delta):
    return math.sqrt(1.0 + (GS**2 - 1.0) * (1.0 + delta) ** 2)


def circ(delta, n):
    def mis(r):
        g = gkin(r, ginf(delta), n)
        return (g * g - 1.0) - g * PS * VS * (RHO / r) ** (n - 1.0)

    a, b = 0.3, 3.0
    fa, fb = mis(a), mis(b)
    for _ in range(70):
        m = 0.5 * (a + b)
        fm = mis(m)
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    r = 0.5 * (a + b)
    return r, gkin(r, ginf(delta), n)


def nux(g, n):
    return math.sqrt(2.0 - n + 1.0 / g**2)


def xi_x(n):
    g = GS
    nu = nux(g, n)
    return B2 * (n - 1.0) / (GS**2 * nu**3)


def lam_over_j(A, n):
    r_out = RHO + A
    g_out = gkin(r_out, GS, n)
    L = r_out * math.sqrt(g_out * g_out - 1.0)

    def f(r):
        g = gkin(r, GS, n)
        return g * g - 1.0 - (L / r) ** 2

    a, b = 0.15, RHO
    for _ in range(80):
        m = 0.5 * (a + b)
        if f(m) > 0:
            b = m
        else:
            a = m
    r_in = 0.5 * (a + b)
    amp = 0.5 * (r_out - r_in)
    mid = 0.5 * (r_out + r_in)
    eps = 1e-9

    def piece(phi, mode):
        r = mid + amp * math.cos(phi)
        fv = f(r)
        pr = math.sqrt(fv if fv > 0.0 else 1e-30)
        g = gkin(r, GS, n)
        jac = amp * abs(math.sin(phi))
        dth = L / (pr * r * r) * jac
        if mode == 0:
            return g / pr * jac
        if mode == 1:
            return dth
        if mode == 2:
            return dth * (r - RHO)
        return dth * (r - RHO) ** 2

    T = quad(lambda p: piece(p, 0), eps, math.pi - eps, epsrel=1e-11, limit=400)[0]
    Th = quad(lambda p: piece(p, 1), eps, math.pi - eps, epsrel=1e-11, limit=400)[0]
    Thx = quad(lambda p: piece(p, 2), eps, math.pi - eps, epsrel=1e-11, limit=400)[0]
    Thx2 = quad(lambda p: piece(p, 3), eps, math.pi - eps, epsrel=1e-11, limit=400)[0]
    xm = Thx / Th
    x2 = Thx2 / Th
    Trev = 2 * math.pi * T / Th
    Ts = 2 * math.pi * RHO / VS
    Lam = (Trev - Ts) / Ts
    nu = nux(GS, n)
    J = (x2 - xm * xm) / (RHO / nu)
    return Lam / J, 2 * math.pi / (2 * Th), xm / J


def prescribed(A, n=1.5):
    nu = nux(GS, n)
    phi = np.linspace(0.0, 2 * math.pi, 30000, endpoint=False)
    x = A * np.cos(phi)
    xp = -A * nu / RHO * np.sin(phi)
    r = RHO + x
    dl = np.sqrt((1.0 + x / RHO) ** 2 + xp**2)
    g = GS - B2 * GS * np.array([psi_r(ri, n) for ri in r])
    v = np.sqrt(1.0 - 1.0 / g**2)
    Lam = np.mean(dl * (VS / v)) - 1.0
    J = (A * A / 2.0) / (RHO / nu)
    series = (J / (2.0 * RHO)) * (nu + (5.0 - n) / (GS**2 * nu))
    return Lam / series


def vertical(A, n=1.5, nturns=400):
    nuy = math.sqrt(n - 1.0)
    nu = nux(GS, n)

    def psi2(r):
        return (n - 1.0) / RHO**2 * (RHO / r) ** (n + 1.0)

    def dpsi0(r):
        return (RHO / r) ** n / RHO

    def dpsi2(r):
        return -(n * n - 1.0) / RHO**3 * (RHO / r) ** (n + 2.0)

    psi = 0.5 * A * A * psi2(RHO)
    g0 = GS - B2 * GS * psi
    L = math.sqrt(g0 * g0 - 1.0) * RHO

    def f(th, y):
        r, pr, z, pz, _t = y
        pphi = L / r
        g = math.sqrt(1.0 + pr * pr + pz * pz + pphi * pphi)
        qEr = -PV * (dpsi0(r) + 0.5 * z * z * dpsi2(r))
        qEz = -PV * (z * psi2(r))
        fac = g * r * r / L
        return [
            pr * r * r / L,
            L / r + qEr * fac,
            pz * r * r / L,
            qEz * fac,
            fac,
        ]

    y = [RHO, 0.0, A, 0.0, 0.0]
    ts, xs = [], []
    for _ in range(nturns):
        sol = solve_ivp(f, (0.0, 2 * math.pi), y, rtol=1e-9, atol=1e-11, method="DOP853")
        y = sol.y[:, -1].copy()
        ts.append(y[4])
        xs.append(y[0])
        y[4] = 0.0
    ts = np.array(ts[nturns // 5 :])
    xs = np.array(xs[nturns // 5 :])
    Ts = 2 * math.pi * RHO / VS
    Lam = (ts.mean() - Ts) / Ts
    J = A * A / (2.0 * (RHO / nuy))
    pref = (J / (2.0 * RHO)) * (nuy + (n - 1.0) / (GS**2 * nuy))
    xterm = (1.0 + 1.0 / GS**2) * (xs.mean() - RHO) / RHO
    return Lam / J, pref / J, xterm / J, (xs.mean() - RHO) / J


if __name__ == "__main__":
    n = 1.5
    d = 1e-4
    rp, gp = circ(d, n)
    rm, gm = circ(-d, n)
    D = (rp - rm) / (2 * d)
    xi_fd = (nux(gp, n) - nux(gm, n)) / (2 * d)
    print("D/rho", D / RHO, "formula", (1 + 1 / GS**2) / nux(GS, n) ** 2)
    print("xi_fd", xi_fd, "xi_form", xi_x(n), "target", -xi_x(n) / RHO)
    print("prescribed ratios", prescribed(1e-4), prescribed(1e-3))
    for A in (0.004, 0.008):
        lj, nu_th, xj = lam_over_j(A, n)
        print(f"n=1.5 A={A} Lam/J={lj:.8f} nu={nu_th:.8f} xnl/J={xj:.5f}")
    lj, nu_th, xj = lam_over_j(0.004, 2.0)
    print(f"n=2 A=0.004 Lam/J={lj:.8f} target={-xi_x(2.0):.8f} xnl/J={xj:.5f}")
    print("vertical", vertical(0.01))
