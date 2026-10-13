"""Convergence in amplitude of the betatron delay in a single-deflector ring.

Units m = c = q = rho = 1. Orbit at total energy of the synchronous particle (w = 0),
started at the outer turning point r0 = 1 + a. Averages are over theta over one radial period.
"""
import numpy as np
from scipy.integrate import solve_ivp


def orbit(nv, gs, a):
    bs2 = 1 - 1/gs**2
    kk = bs2*gs
    Phi = (lambda R: kk*np.log(R)) if nv == 1 else (lambda R: kk*(R**(1 - nv) - 1)/(1 - nv))
    r0 = 1 + a
    L = r0*np.sqrt((gs - Phi(r0))**2 - 1)

    def f(th, y):
        R, pr = y[0], y[1]
        gam = np.sqrt(1 + pr**2 + L**2/R**2)
        dth_dt = L/(gam*R**2)
        rp = (pr/gam)/dth_dt
        return [rp, (L**2/(gam*R**3) - kk*R**(-nv))/dth_dt, 1/dth_dt, R - 1, rp**2]

    ev = lambda th, y: y[1]
    ev.direction = -1
    s = solve_ivp(f, [0, 40], [r0, 0, 0, 0, 0], events=ev, rtol=1e-13, atol=1e-15,
                  method='DOP853', dense_output=True)
    th_e = [te for te in s.t_events[0] if te > 1e-3][0]
    ys = s.sol(np.linspace(0, th_e, 20001))
    R, pr, t, I_x, I_rp2 = s.sol(th_e)
    rmin, rmax = ys[0].min(), ys[0].max()
    T = 2*np.pi*t/th_e
    return dict(Lam=T/(2*np.pi/np.sqrt(bs2)) - 1, xnl=I_x/th_e, xp2=I_rp2/th_e,
                X=((rmax - rmin)/2)**2/2, nu=2*np.pi/th_e)


gs = 1.2483
bs2 = 1 - 1/gs**2
for nv in (1.0, 1.5, 2.0):
    print(f'n={nv}')
    for a in (4e-3, 2e-3, 1e-3, 5e-4):
        o = orbit(nv, gs, a)
        gen = (1 + 1/gs**2)*o['xnl'] + o['xp2']/2 + (5 - nv)/(2*gs**2)*o['X']
        print(f'  a={a:.0e}: Lambda/X={o["Lam"]/o["X"]: .6f}  general/X={gen/o["X"]: .6f}'
              f'  xnl/X={o["xnl"]/o["X"]: .6f}  xp2/X={o["xp2"]/o["X"]: .6f}  nu^2={o["nu"]**2:.6f}')
    if nv == 1.0:
        print(f'  table: xnl/X = {-(gs**2+5)/(2*(gs**2+1)):.6f}')
    if nv == 2.0:
        print(f'  table: xnl/X = -1, Lambda/X = -beta^2 = {-bs2:.6f}')
