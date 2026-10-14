"""Control check of the longitudinal system in the variable delta.

System (t in units of 1/omega_s, small-phase form):
    phi' = r (eta0 d + eta1 d^2 + Lam),     d' = -(1/(r eta0)) phi (1 - d/g^2)
r = omega_rf/omega_s. The factor (1 - d/g^2) = v_s/v is the exact gap kick in delta.
Expected (supervisor):
    d = d_eq + dm cos(psi) + (eta1/(6 eta0) - 1/(3 g^2)) dm^2 cos(2 psi)
    d_eq = -(eta1 dm^2/2 + Lam)/eta0
    Omega/omega_s = 1 - et Lam/eta0^2 - (5/12)(et dm/eta0)^2,  et = eta1 - eta0/(2 g^2)
"""
import numpy as np
from scipy.integrate import solve_ivp

r, eta0, eta1, g = 3.0, -0.53, 1.02, 1.1439
et = eta1 - eta0 / (2 * g**2)


def run(dm_guess, lam):
    d_c = -lam / eta0                       # stable point at small amplitude
    def f(t, y):
        phi, d = y
        return [r * (eta0 * d + eta1 * d**2 + lam), -phi * (1 - d / g**2) / (r * eta0)]
    ev = lambda t, y: y[0]
    ev.direction = 1 if eta0 > 0 else -1
    sol = solve_ivp(f, [0, 400], [0.0, d_c + dm_guess], rtol=1e-12, atol=1e-14,
                    dense_output=True, events=ev, max_step=0.05)
    te = sol.t_events[0]
    t0, t1 = te[-2], te[-1]
    T = (t1 - t0)
    tt = np.linspace(t0, t1, 20001)[:-1]
    d = sol.sol(tt)[1]
    psi = 2 * np.pi * (tt - t0) / T
    dmean = d.mean()
    dm = (d.max() - d.min()) / 2
    c2 = 2 * np.mean((d - dmean) * np.cos(2 * psi)) * np.sign(np.mean(d * np.cos(psi)))
    return dmean, dm, c2, 2 * np.pi / T


for lam in (0.0, 2e-5):
    for a in (0.004, 0.008):
        dmean, dm, c2, Om = run(a, lam)
        print(f"Lam={lam:g} dm={dm:.5f}: "
              f"d_eq num={dmean:.4e} th={-(eta1*dm**2/2 + lam)/eta0:.4e}; "
              f"c2 num={c2:.4e} th={(eta1/(6*eta0) - 1/(3*g**2))*dm**2:.4e}; "
              f"Om-1 num={Om-1:.4e} th={-et*lam/eta0**2 - 5/12*(et*dm/eta0)**2:.4e}")
