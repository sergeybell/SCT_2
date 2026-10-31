"""Electrostatic deflector: checks of chapters 4 and 5 from one Hamiltonian.

H = -(1+x/r) sqrt(Pi^2 - px^2 - py^2), momenta normalised by p_s,
Pi^2 = kinetic momentum^2 from exact total-energy conservation, delta = outer momentum.
Time: v_s dt/ds = -(v_s/v_out) dH/d delta.
Checks:
  1. potential (5.12) satisfies Laplace in sector coordinates;
  2. closed-orbit integrand (4.16) == v_s dt/ds on x = D1 d + D2 d^2;
  3. dispersion equation (4.10);
  4. second-order force (5.19)-(5.20);
  5. chromaticity integrands (5.22)-(5.23);
  6. Lambda integrand + J*I/2 == -d/ds(2 Jx D1 ax / r)   (5.25);
  7. orbit kink at the face: x'_in - x'_out = 2 x x'/r (also for D2 at order d^2).
"""
import sympy as sp

x, px, y, py, d, r, e, n = sp.symbols('x p_x y p_y delta rho epsilon n', real=True)
g = sp.symbols('gamma_s', positive=True)
X, Px, Y, Py = sp.symbols('X P_X Y P_Y', real=True)
D, Dp, D2, D2p = sp.symbols('D1 Dp1 D2 Dp2', real=True)
Jx, Jy, bx, by, ax, ay = sp.symbols('J_x J_y beta_x beta_y alpha_x alpha_y', positive=True)
a, b = sp.symbols('psi_x psi_y', real=True)
b2 = 1 - 1/g**2


def taylor(expr, order):
    out, cur = 0, expr
    for k in range(order):
        out += sp.simplify(cur.subs(e, 0))*e**k/sp.factorial(k)
        cur = sp.diff(cur, e)
    return out


chi = x/r - n*x**2/(2*r**2) + (n - 1)*y**2/(2*r**2) + n*(n + 1)*x**3/(6*r**3) + (1 - n**2)*x*y**2/(2*r**3)

# 1. Laplace in sector coordinates, through second order
lap = sp.diff(chi, x, 2) + sp.diff(chi, y, 2) + sp.diff(chi, x)/(r + x)
lap = taylor(lap.subs({x: e*x, y: e*y}), 2)
assert sp.simplify(sp.expand(lap)) == 0

gout = g*sp.sqrt(1 + b2*(2*d + d**2))
gam = gout - b2*g*chi
Pi2 = (gam**2 - 1)/(g**2*b2)
H = -(1 + x/r)*sp.sqrt(Pi2 - px**2 - py**2)
vs_over_vout = (gout/(1 + d))/g          # v_s/v_out = gamma_out/(gamma_s (1+d))


def ser(expr, subs, order):
    return sp.expand(taylor(expr.subs(subs, simultaneous=True), order).subs(e, 1))

# 2. closed orbit time of flight
xs, xps = sp.symbols("x_c x'_c", real=True)
Pi2_xy0 = Pi2.subs(y, 0)
px_of = xps*sp.sqrt(Pi2_xy0)/sp.sqrt((1 + xs/r)**2 + xps**2)
dtds = (-vs_over_vout*sp.diff(H, d)).subs(py, 0).subs(y, 0)
dtds = dtds.subs(px, px_of.subs({}, simultaneous=True)).subs({x: xs})
orb = {xs: e*D*d + e**2*D2*d**2, xps: e*Dp*d + e**2*D2p*d**2, d: e*d}
lhs2 = ser(dtds.subs(x, xs), orb, 3)
rhs2 = (1 + ((1 + 1/g**2)*D/r - 1/g**2)*d
        + (Dp**2/2 + (1 + 1/g**2)*D2/r + (5 - n)/(2*g**2)*(D/r)**2 - 4/g**2*D/r
           + (3*g**2 - 1)/(2*g**4))*d**2)
assert sp.simplify(lhs2 - sp.expand(rhs2)) == 0
print('(4.16) ok')

# 3. linear equations: x' = dH/dpx, px' = -dH/dx
lin = {x: e*x, px: e*px, y: e*y, py: e*py, d: e*d}
xp1 = ser(sp.diff(H, px), lin, 2)
pxp1 = ser(-sp.diff(H, x), lin, 2)
Kx = (2 - n + 1/g**2)/r**2
print("x' =", sp.simplify(xp1), "  px' =", sp.simplify(pxp1))
assert sp.simplify(pxp1 - (-Kx*x + (1 + 1/g**2)*d/r)) == 0 and sp.simplify(xp1 - px) == 0
print('(4.10) ok')

# 4. second-order force at d = 0 in geometric slopes
H0 = H.subs(d, 0)
q3 = {x: e*x, px: e*px, y: e*y, py: e*py}
f = {v: ser(expr, q3, 3) for v, expr in
     [(x, sp.diff(H0, px)), (px, -sp.diff(H0, x)), (y, sp.diff(H0, py)), (py, -sp.diff(H0, y))]}
xpp = ser(sum(sp.diff(f[x], v)*f[v] for v in (x, px, y, py)), q3, 3)
ypp = ser(sum(sp.diff(f[y], v)*f[v] for v in (x, px, y, py)), q3, 3)
u, w = sp.symbols("x' y'", real=True)
inv = {px: u - 2*x*u/r, py: w - 2*x*w/r}     # px = x'(1 - 2x/r) at second order
geo = {x: e*x, u: e*u, y: e*y, w: e*w}
Fx = ser(xpp.subs(inv, simultaneous=True) + Kx.subs(g, g)*x, geo, 3)
Fy = ser(ypp.subs(inv, simultaneous=True) + (n - 1)/r**2*y, geo, 3)
Axx = -(g**2*(n**2 - 6*n + 8) - 3*n + 10)/(2*g**2*r**3)
Ayy = (n - 1)*(n*g**2 - 1)/(2*g**2*r**3)
assert sp.simplify(Fx - (Axx*x**2 + Ayy*y**2 + (u**2 - w**2)/r)) == 0
Axy = (n - 1)*(n - 2)/r**3 + (1 - n)/(g**2*r**3)
assert sp.simplify(Fy - (2*u*w/r + Axy*x*y)) == 0
print('(5.19)-(5.20) ok; F^y also has an xy term, phase average 0')

# 5. chromaticity integrands
H2 = sp.diff(H.subs({x: D*d + e*X, px: Dp*d + e*Px, y: e*Y, py: e*Py}, simultaneous=True), e, 2).subs(e, 0)/2
H21 = sp.expand(sp.simplify(sp.diff(H2, d).subs(d, 0)))
tw = {X: sp.sqrt(2*Jx*bx)*sp.cos(a), Px: -sp.sqrt(2*Jx/bx)*(sp.sin(a) + ax*sp.cos(a)),
      Y: sp.sqrt(2*Jy*by)*sp.cos(b), Py: -sp.sqrt(2*Jy/by)*(sp.sin(b) + ay*sp.cos(b))}
def avg(expr, sub):
    h = sp.expand(sp.expand_trig(expr.subs(sub, simultaneous=True)))
    h = sp.integrate(h, (a, 0, 2*sp.pi))/(2*sp.pi)
    return sp.simplify(sp.integrate(h, (b, 0, 2*sp.pi))/(2*sp.pi))
A21 = avg(H21, tw)
Ix = sp.simplify(2*sp.diff(A21, Jx)); Iy = sp.simplify(2*sp.diff(A21, Jy))
gx = (1 + ax**2)/bx; gy = (1 + ay**2)/by
Ix_t = -gx - (5 - n)*bx/(g**2*r**2) + D/r*(2*gx + (n - 2)*(n*g**2 - 3)*bx/(g**2*r**2)) - 4*ax*Dp/r
Iy_t = -gy - (n - 1)*by/(g**2*r**2) + D/r*(2*gy - (n - 1)*(n*g**2 - 1)*by/(g**2*r**2))
assert sp.simplify(Ix - Ix_t) == 0 and sp.simplify(Iy - Iy_t) == 0
print('(5.22)-(5.23) ok')

# 6. Lambda integrand vs chromaticity inside the deflector
tg = {x: sp.sqrt(2*Jx*bx)*sp.cos(a), u: -sp.sqrt(2*Jx/bx)*(sp.sin(a) + ax*sp.cos(a)),
      y: sp.sqrt(2*Jy*by)*sp.cos(b), w: -sp.sqrt(2*Jy/by)*(sp.sin(b) + ay*sp.cos(b))}
Fav = avg(Fx, tg)
Lint = Jx*gx/2 + Jy*gy/2 + (5 - n)*Jx*bx/(2*g**2*r**2) + (n - 1)*Jy*by/(2*g**2*r**2) + Fav*D
dax = Kx*bx - gx
assert sp.simplify(Lint + (Jx*Ix + Jy*Iy)/2 + 2*Jx/r*(Dp*ax + D*dax)) == 0
print('(5.25) ok: Lambda-integrand + J I/2 = -d/ds(2 Jx D1 ax/rho)')

# 7. kink at the face from continuity of px
xo, xi_ = sp.symbols("x'_out x'_in", real=True)
px_in = (xi_*sp.sqrt(Pi2)/sp.sqrt((1 + x/r)**2 + xi_**2)).subs(y, 0)
px_out = (xo*sp.sqrt(Pi2.subs({x: 0, y: 0}))/sp.sqrt(1 + xo**2))
k = sp.symbols('k', real=True)
sol = ser((px_in - px_out).subs(xi_, xo + k), {x: e*x, xo: e*xo, k: e**2*k, d: e*d}, 3)
kk = sp.solve(sp.expand(sol), k)[0]
print("x'_in - x'_out =", sp.simplify(kk))
assert sp.simplify(kk - 2*x*xo/r) == 0
print('jump ok; on the closed orbit at order d^2: D2\'_in - D2\'_out = 2 D1 D1\'/rho')
print('average at entrance (after - before):', avg(2*x*u/r, tg))
