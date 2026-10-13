"""Independent symbolic checks of the control formulas (chapters 1-3)."""
import sympy as sp
import mpmath as mp

d, g, a0, a1 = sp.symbols('delta gamma alpha0 alpha1', positive=True)
b2 = 1 - 1/g**2

# 1. tau(delta) for a magnetic ring: T/Ts = (C/Cs)(vs/v)
bg_s = sp.sqrt(g**2 - 1)
gam = sp.sqrt(1 + bg_s**2*(1 + d)**2)
vs_over_v = gam/(g*(1 + d))
tau = sp.expand(sp.series((1 + a0*d + a1*d**2)*vs_over_v - 1, d, 0, 3).removeO())
eta0 = sp.simplify(tau.coeff(d, 1))
eta1 = sp.simplify(tau.coeff(d, 2))
print('eta0 =', eta0, ' check:', sp.simplify(eta0 - (a0 - 1/g**2)))
print('eta1 check:', sp.simplify(eta1 - (a1 - a0/g**2 + (1 + b2/2)/g**2)))

# 2. w = dE/(beta^2 E) as series in delta
E = gam  # E/(m c^2)
w = sp.series((E - g)/(b2*g), d, 0, 3).removeO()
print('w - (d + d^2/(2g^2)) =', sp.simplify(sp.expand(w) - (d + d**2/(2*g**2))))

# 3. Lindstedt solution of  x'' + w0^2 x = c2 x^2  (x = delta), to second order
eps, A, th, w0, c2, W2, b0, b2c, a3, b1c = sp.symbols('epsilon A theta omega0 c2 W2 b0 b2c a3 b1c')
Om = w0*(1 + eps**2*W2)
x = eps*A*sp.cos(th) + eps**2*(b0 + b2c*sp.cos(2*th)) + eps**3*(a3*sp.cos(3*th) + b1c*sp.cos(th))
res = sp.expand(Om**2*sp.diff(x, th, 2) + w0**2*x - c2*x**2)
res2 = sp.expand(sp.fu(res.coeff(eps, 2), measure=lambda e: e.count_ops()))
sol2 = sp.solve([sp.integrate(res.coeff(eps, 2), (th, 0, 2*sp.pi)),
                 sp.integrate(res.coeff(eps, 2)*sp.cos(2*th), (th, 0, 2*sp.pi))], [b0, b2c])
r3 = res.coeff(eps, 3).subs(sol2)
sol3 = sp.solve(sp.integrate(r3*sp.cos(th), (th, 0, 2*sp.pi)), W2)
print('b0 =', sol2[b0], ' b2 =', sol2[b2c], ' W2 =', sp.simplify(sol3[0]))
e0, e1 = sp.symbols('eta0 eta1')
sub = {c2: -w0**2*e1/e0}
print('  with c2=-w0^2 eta1/eta0: delta_eq =', sp.simplify(sol2[b0].subs(sub)),
      ' 2nd harm =', sp.simplify(sol2[b2c].subs(sub)), ' Omega2 =', sp.simplify(sol3[0].subs(sub)))
# half-span: x(0) - x(pi) over 2 at eps^2 order
xs = (eps*A*sp.cos(th) + eps**2*(sol2[b0] + sol2[b2c]*sp.cos(2*th)))
print('  half-span - A =', sp.simplify((xs.subs(th, 0) - xs.subs(th, sp.pi))/2 - eps*A))

# 4. phi(t): phi' = w_rf tau, tau = eta0 x + eta1 x^2, time-average zero
wrf = sp.symbols('omega_rf')
xd = (A*sp.cos(th) + sol2[b0].subs(sub) + sol2[b2c].subs(sub)*sp.cos(2*th))
taud = sp.expand(e0*xd + e1*xd**2)
taud2 = sp.fourier_series  # placeholder to keep namespace tidy
c_const = sp.simplify(sp.integrate(taud, (th, 0, 2*sp.pi))/(2*sp.pi))
c_cos2 = sp.simplify(sp.integrate(taud*sp.cos(2*th), (th, 0, 2*sp.pi))/sp.pi)
print('  <tau> through A^2 (should be O(A^3)) =', sp.series(c_const, A, 0, 3))
print('  tau cos2 coeff =', sp.series(c_cos2, A, 0, 3), ' -> phi sin2 coeff = (w_rf/Om)*that/2')

# 5. w_eq
print('w_eq = delta_eq + A^2/(4 g^2) = -(eta1 - eta0/(2g^2)) A^2/(2 eta0):',
      sp.simplify(-e1*A**2/(2*e0) + A**2/(4*g**2) + (e1 - e0/(2*g**2))*A**2/(2*e0)))

# 6. Pendulum: <delta^2> in units of the half-span A, vs phase half-amplitude
mp.mp.dps = 30
def ratio(phihat):
    k2 = mp.sin(phihat/2)**2
    return (mp.ellipe(k2)/mp.ellipk(k2) - 1 + k2)/k2
def ratio_quad(phihat):
    # direct time average over a pendulum orbit: phi'' = -sin phi, delta ~ phi'
    T = 4*mp.ellipk(mp.sin(phihat/2)**2)
    f = lambda p: mp.sqrt(2*(mp.cos(p) - mp.cos(phihat)))
    num = 4*mp.quad(f, [0, phihat])            # \oint phi' dphi
    A2 = 2*(1 - mp.cos(phihat))                 # max phi'^2
    return num/T/A2
for ph in [0.05, 0.2, 0.5]:
    r = ratio(ph)
    print(f'phihat={ph}: elliptic {mp.nstr(r,12)}, quadrature {mp.nstr(ratio_quad(ph),12)},'
          f' 1/2(1-phi^2/32) {mp.nstr(0.5*(1-ph**2/32),12)}')
