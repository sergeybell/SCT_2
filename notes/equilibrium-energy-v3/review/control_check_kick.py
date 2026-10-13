"""Exact gap kick: delta-system with the factor v_s/v versus the w-system.

w-system (exact kick): phi' = w_rf (eta0 w + et1 w^2), w' = Kc sin(phi), et1 = eta1 - eta0/(2 g^2).
Check: delta = w - w^2/(2 g^2) has second harmonic eta1 A^2/(6 eta0) - A^2/(3 g^2),
and the delta-system with dE = beta^2 E d(delta) (no v_s/v factor) misses -A^2/(3 g^2).
Numerical integration of the exact w-system, harmonic analysis of delta(t).
"""
import numpy as np
from scipy.integrate import solve_ivp

g, eta0, eta1, mu, A0 = 1.2, -0.4, 0.7, 1.0, 1e-2     # mu = w_rf/w_s, time in 1/w_s
et1 = eta1 - eta0/(2*g**2)
Kc = -1/(mu*eta0)                                       # K cos(phi_s) w_rf = -w_s^2/eta0


def rhs(t, y):
    phi, w = y
    return [mu*(eta0*w + et1*w**2), Kc*np.sin(phi)]


sol = solve_ivp(rhs, [0, 40], [0, A0], rtol=1e-12, atol=1e-15, dense_output=True,
                events=lambda t, y: y[0], method='DOP853')
tz = sol.t_events[0]
T = tz[2] - tz[0]                                        # two zero crossings of phi = one period
t = np.linspace(tz[0], tz[0] + T, 40001)[:-1]
w = sol.sol(t)[1]
d = w - w**2/(2*g**2) + w**3/(2*g**2)*0                  # delta(w) to second order
psi = 2*np.pi*(t - tz[0])/T
c1 = 2*np.mean(d*np.cos(psi)); s1 = 2*np.mean(d*np.sin(psi))
A = np.hypot(c1, s1); ph = np.arctan2(s1, c1)
c2 = 2*np.mean(d*np.cos(2*(psi - ph)))
print('A =', A)
print('mean delta  :', d.mean(), ' -eta1 A^2/(2 eta0) =', -eta1*A**2/(2*eta0))
print('2nd harm    :', c2, ' eta1 A^2/(6eta0) - A^2/(3g^2) =', eta1*A**2/(6*eta0) - A**2/(3*g**2),
      ' (without kick factor:', eta1*A**2/(6*eta0), ')')
print('Omega/ws - 1:', 2*np.pi/T - 1, ' -5/12 (et1 A/eta0)^2 =', -5/12*(et1*A/eta0)**2,
      ' (with eta1:', -5/12*(eta1*A/eta0)**2, ')')
