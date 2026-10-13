"""Independent checks of the electrostatic formulas (chapters 4-5).

Units m = c = q = 1, reference radius rho = 1.
Deflector potential psi = q*Phi/(p_s v_s) = xi - n/2 xi^2 + (n-1)/2 eta^2.
Exact local Lorentz factor: gamma = gamma_s (1 + beta_s^2 (w - psi)).
"""
import sympy as sp
import mpmath as mp

# ---------- A. v_s/v as a series in u = w - psi ----------
u, g = sp.symbols('u gamma', positive=True)
b2 = 1 - 1/g**2
gl = g*(1 + b2*u)
vs_v = sp.sqrt(b2)/sp.sqrt(1 - 1/gl**2)
ser = sp.series(vs_v, u, 0, 3).removeO()
a1 = sp.simplify(ser.coeff(u, 1))
a2 = sp.simplify(ser.coeff(u, 2))
print('v_s/v = 1 + a1 u + a2 u^2:  a1 =', a1, '  a2 =', sp.factor(a2))

# ---------- B. eta0, eta1 for a ring with electric bends (index n) ----------
d, n = sp.symbols('delta n')
D, D2, Dp, rho_inv = sp.symbols('D D2 Dp rhoinv')   # rhoinv = 1/rho (0 outside deflectors)
w = d + d**2/(2*g**2)
x = D*d + D2*d**2
xi = x*rho_inv
psi = xi - n/2*xi**2
uu = w - psi
F = (1 + a1*uu + a2*uu**2)*(1 + xi + (Dp*d)**2/2) - 1
F = sp.expand(sp.series(F, d, 0, 3).removeO())
e0 = sp.collect(sp.expand(F.coeff(d, 1)), [D, rho_inv])
e1 = sp.expand(F.coeff(d, 2))
print('integrand coeff of delta  :', e0)
print('integrand coeff of delta^2:', sp.collect(e1, [rho_inv], sp.factor))

# averaged symbols: <D/rho>=al0, <D2/rho>=m2, <(D/rho)^2>=q2, <Dp^2>=p2
al0, m2, q2, p2 = sp.symbols('alpha0 m2 q2 p2')
def avg(expr):
    expr = sp.expand(expr)
    expr = expr.subs({D**2*rho_inv**2: q2}).subs({D2*rho_inv: m2, D*rho_inv: al0, Dp**2: p2})
    return sp.expand(expr)
eta0 = sp.simplify(avg(e0) - 0)   # constant terms average to themselves
eta1 = sp.expand(avg(e1))
print('eta0 =', sp.factor(eta0))
print('eta1 =', sp.collect(eta1, [m2, q2, al0, p2], sp.factor))
alpha1 = m2 + p2/2
eta1_mag = alpha1 - al0/g**2 + (1 + b2/2)/g**2
old = eta1_mag + (m2 + (5 - n)/2*q2 - 3*al0)/g**2
print('eta1 - old formula (es-eta1) =', sp.simplify(eta1 - old))

# ---------- C. exact circular orbits in a single-deflector ring ----------
# force F(r) = k r^-n, k = p_s v_s = beta_s^2 gamma_s; Phi'(r) = k r^-n, Phi(1)=0
# Exact rational checks for Pythagorean gamma (beta rational) and several n.
# Parametrize circular orbits by the local Lorentz factor G: p v = G - 1/G = k r^(1-n).
Gs = sp.symbols('G', positive=True)
def single_deflector(nv, gval):
    bv2 = 1 - 1/gval**2
    kk = bv2*gval
    r_G = ((Gs - 1/Gs)/kk)**(1/(1 - sp.nsimplify(nv))) if nv != 1 else None
    # n = 1: p v constant on all circular orbits -> G = gamma_s, r free; handle separately
    if nv == 1:
        xr_ = sp.symbols('x')
        # w = ln r exactly (Phi = k ln r, gamma unchanged), tau = r - 1
        w_of_x = sp.log(1 + xr_)
        tau_of_x = xr_
    else:
        xr_ = sp.symbols('x')
        rr = 1 + xr_
        G_r = ((kk*rr**(1 - nv)) + sp.sqrt((kk*rr**(1 - nv))**2 + 4))/2
        Phi_r = kk*(rr**(1 - nv) - 1)/(1 - nv)
        w_of_x = (G_r + Phi_r - gval)/(bv2*gval)
        v_r = sp.sqrt(1 - 1/G_r**2)
        tau_of_x = rr/v_r*sp.sqrt(bv2) - 1
    tay = lambda e: [sp.nsimplify(sp.simplify(sp.diff(e, xr_, j).subs(xr_, 0)/sp.factorial(j)))
                     for j in range(3)]
    W, Tt = tay(w_of_x), tay(tau_of_x)
    cc1, cc2 = sp.symbols('cc1 cc2')
    xw = cc1*ww + cc2*ww**2
    eqw = sp.expand(W[1]*xw + W[2]*xw**2 - ww)
    sl = sp.solve([eqw.coeff(ww, 1), eqw.coeff(ww, 2)], [cc1, cc2], dict=True)[0]
    xw = xw.subs(sl)
    tau_w = sp.expand(Tt[1]*xw + Tt[2]*xw**2)
    gs_ = sp.nsimplify(gval)
    tau_d = sp.expand(tau_w.subs(ww, d + d**2/(2*gs_**2)))
    x_d = sp.expand(xw.subs(ww, d + d**2/(2*gs_**2)))
    e0x = tau_d.coeff(d, 1); e1x = tau_d.coeff(d, 2)
    Dv, D2v = x_d.coeff(d, 1), x_d.coeff(d, 2)
    sub = {al0: Dv, m2: D2v, q2: Dv**2, p2: 0, g: gs_, n: sp.nsimplify(nv)}
    return (sp.nsimplify(e0x), sp.nsimplify(e1x),
            sp.simplify(eta0.subs(sub) - e0x), sp.simplify(eta1.subs(sub) - e1x))

ww = sp.symbols('w')
print('\nsingle deflector ring, exact circular orbits vs general formula:')
for nv in (1, sp.Rational(3, 2), 2):
    for gval in (sp.Rational(5, 4), sp.Rational(5, 3)):
        e0x, e1x, r0_, r1_ = single_deflector(nv, gval)
        print(f'  n={nv}, gamma={gval}: eta0={e0x}, eta1={e1x}; general-exact: {r0_}, {r1_}')
# closed forms of the old table for n=1, n=2
for gval in (sp.Rational(5, 4), sp.Rational(5, 3)):
    print(f'  old table at gamma={gval}: n=1 eta1=(g^2+1)/(2g^2)={(gval**2+1)/(2*gval**2)},'
          f' n=2 eta0=g^2+2={gval**2+2}, eta1=(9g^4+7g^2+2)/(2g^2)={(9*gval**4+7*gval**2+2)/(2*gval**2)}')

# ---------- D. betatron delay: general formula (symbolic) ----------
xb, xnl, yb, xp, yp, X, Y, XNL, XP2 = sp.symbols('xb xnl yb xp yp X Y XNL XP2')
eps = sp.symbols('epsilon')
xi_b = (eps*xb + eps**2*xnl)*rho_inv
eta_b = eps*yb*rho_inv
psi_b = xi_b - n/2*xi_b**2 + (n - 1)/2*eta_b**2
Fb = (1 + a1*(-psi_b) + a2*psi_b**2)*(1 + xi_b + eps**2*(xp**2 + yp**2)/2) - 1
Fb = sp.expand(sp.series(Fb, eps, 0, 3).removeO())
Fb = Fb.subs(eps, 1)
# betatron-phase average: <xb>=0, <xb^2>=X/rhoinv^2 ...
Fb = sp.expand(Fb).subs({xb**2: X/rho_inv**2, yb**2: Y/rho_inv**2}).subs({xb: 0, yb: 0})
Fb = sp.expand(Fb.subs({xnl*rho_inv: XNL}))
print('\nLambda integrand:', sp.collect(Fb, [XNL, X, Y], sp.factor))
old_L = (1 + 1/g**2)*XNL + (xp**2 + yp**2)/2 + (5 - n)/(2*g**2)*X + (n - 1)/(2*g**2)*Y
print('Lambda integrand - old (es-lambda):', sp.simplify(Fb - old_L))

# ---------- E. single-deflector ring: Lambda from numerical orbits ----------
import numpy as np
from scipy.integrate import solve_ivp

def lam_numeric(nv, gs, a):
    bs2 = 1 - 1/gs**2
    kk = bs2*gs
    Phi_ = (lambda R: kk*np.log(R)) if nv == 1 else (lambda R: kk*(R**(1 - nv) - 1)/(1 - nv))
    r0 = 1 + a
    gam0 = gs - Phi_(r0)          # total energy = gamma_s (w = 0)
    L = r0*np.sqrt(gam0**2 - 1)
    def f(th, y):                 # independent variable theta
        R, pr, t, I_r, I_rp2 = y
        gam = np.sqrt(1 + pr**2 + L**2/R**2)
        dth_dt = L/(gam*R**2)
        drdt = pr/gam
        dprdt = L**2/(gam*R**3) - kk*R**(-nv)
        rp = drdt/dth_dt
        return [rp, dprdt/dth_dt, 1/dth_dt, R, rp**2]
    ev = lambda th, y: y[1]
    ev.direction = -1             # p_r goes + -> - : outer turning point
    s = solve_ivp(f, [0, 40], [r0, 0, 0, 0, 0], events=ev, rtol=1e-12, atol=1e-14,
                  method='DOP853')
    th_e = [te for te in s.t_events[0] if te > 1e-6][0]
    ye = s.sol(th_e) if s.sol is not None else None
    s2 = solve_ivp(f, [0, th_e], [r0, 0, 0, 0, 0], rtol=1e-12, atol=1e-14, method='DOP853')
    R, pr, t, I_r, I_rp2 = s2.y[:, -1]
    T = 2*np.pi*t/th_e
    Ts = 2*np.pi/np.sqrt(bs2)
    rmin = s2.y[0].min()
    return dict(Lam=T/Ts - 1, xnl=I_r/th_e - 1, xp2=I_rp2/th_e, nu=2*np.pi/th_e,
                X=((r0 - rmin)/2)**2/2)

print('\nsingle deflector, betatron delay from orbits (gamma_s = 1.2483):')
gs = 1.2483
bs2 = 1 - 1/gs**2
for nv in (1.0, 1.5, 2.0):
    r1, r2 = lam_numeric(nv, gs, 1e-3), lam_numeric(nv, gs, 2e-3)
    # remove O(a^3) by Richardson on a^2 scaling
    cL = (4*r1['Lam'] - r2['Lam'])/3/r1['X'] if False else r1['Lam']/r1['X']
    gen = ((1 + 1/gs**2)*r1['xnl'] + r1['xp2']/2 + (5 - nv)/(2*gs**2)*r1['X'])/r1['X']
    # chromaticity form: Lambda = -X nu_x dnu_x/dw, nu_x^2 = 3 - n - beta_loc^2(w)
    wv = sp.symbols('w')
    print(f'  n={nv}: Lambda/X orbits = {cL:.6f}, general formula with orbit x_nl = {gen:.6f},'
          f' nu_x = {r1["nu"]:.8f} vs sqrt(3-n-beta^2) = {np.sqrt(3-nv-bs2):.8f}')

# chromaticity form: Lambda = -X nu_x xi_x, xi_x = d nu_x / dw, nu_x^2 = 3 - n - beta_loc^2(w)
# numerically: circular orbit with total energy w -> local gamma -> beta^2
from scipy.optimize import brentq
def beta2_loc(nv, gs, wv):
    bs2 = 1 - 1/gs**2; kk = bs2*gs
    Phi_ = (lambda R: kk*np.log(R)) if nv == 1 else (lambda R: kk*(R**(1 - nv) - 1)/(1 - nv))
    Gr = lambda R: (kk*R**(1 - nv) + np.sqrt((kk*R**(1 - nv))**2 + 4))/2
    Etot = gs + bs2*gs*wv
    R = brentq(lambda R: Gr(R) + Phi_(R) - Etot, 0.5, 2.0) if nv != 1 else np.exp(wv)
    return 1 - 1/Gr(R)**2
for nv in (1.0, 1.5, 2.0):
    h = 1e-5
    nu = lambda wv: np.sqrt(3 - nv - beta2_loc(nv, gs, wv))
    xi_x = (nu(h) - nu(-h))/(2*h)
    print(f'  n={nv}: chromaticity form Lambda/X = -nu_x xi_x = {-nu(0)*xi_x:.6f}')
