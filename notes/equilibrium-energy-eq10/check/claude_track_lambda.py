"""Independent supervisor check of the betatron delay in smooth deflector rings.

Planar motion in a radial field E ~ r^-k (k=1 cylindrical, k=2 spherical deflector),
units m=c=rho=1, total energy equal to the reference one (delta=0).
A particle starts at r=1+a with zero radial velocity; tau is the relative excess of the
time per 2*pi of azimuth over T_s = 2*pi/beta_s. Expected: tau = c * <(x/rho)^2>.
"""
import numpy as np
from scipy.integrate import solve_ivp


def run(k, g, a, turns=40):
    bs = np.sqrt(1 - 1 / g**2)
    qE0 = g * bs**2  # p_s v_s / rho
    qPhi = (lambda r: qE0 * np.log(r)) if k == 1 else (lambda r: qE0 * (r**(1 - k) - 1) / (1 - k))
    E_tot = g
    r0 = 1 + a
    G0 = E_tot - qPhi(r0)
    p0 = np.sqrt(G0**2 - 1)

    def rhs(t, s):
        x, y, px, py = s
        r = np.hypot(x, y)
        G = np.sqrt(1 + px**2 + py**2)
        F = qE0 * r**(-k)
        return [px / G, py / G, F * x / r * -1, F * y / r * -1]

    def cross(t, s):
        return np.unwrap([0.0])[0] if False else s[1]

    state = [r0, 0.0, 0.0, p0]
    theta_target = 2 * np.pi * turns

    def ang(t, s):
        return s[1]

    sol = solve_ivp(rhs, [0, 2 * np.pi * turns / bs * 1.2], state, rtol=1e-12, atol=1e-14,
                    dense_output=True, max_step=0.05)
    t = sol.t
    th = np.unwrap(np.arctan2(sol.y[1], sol.y[0]))
    T = np.interp(theta_target, th, t)
    grid = np.linspace(0, theta_target, 200001)
    x2 = np.mean((np.interp(grid, th, np.hypot(sol.y[0], sol.y[1])) - 1)**2)
    return T / (2 * np.pi * turns / bs) - 1, x2


for k, name in [(1, "cylindrical"), (2, "spherical")]:
    for g in [1.25, 2.0]:
        bs2 = 1 - 1 / g**2
        a = 1e-3
        (t1, x1), (t2, x2) = run(k, g, a), run(k, g, 2 * a)
        # end-phase error is linear in a: tau = A a + B a^2, while <x^2> = q a^2
        quad = (t2 - 2 * t1) / 2
        print(f"{name:11s} gamma={g}: tau_quad/<x^2> = {quad / x1:+.5f};  -beta_s^2 = {-bs2:+.5f}")
