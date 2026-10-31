"""Sector dipole: Lambda vs -2*pi*J*xi, term by term, with hard-edge boundary terms.

Hamiltonian (s independent, momenta normalised by p_s, uniform field, perpendicular faces):
    H = -(1+x/r) sqrt((1+d)^2 - px^2 - py^2) + x/r + x^2/(2 r^2).
Checks:
  1. linear equation and dK_x/dd;
  2. second-order force F^x at d = 0 in geometric slopes;
  3. dipole chromaticity integrands I_x, I_y from the canonical quadratic Hamiltonian;
  4. Lambda-integrand + J*I/2 + (Jx ax' + Jy ay')/2 == -d/ds(Jx D1 ax / r) inside the dipole;
  5. phase-averaged slope jump at the faces reproduces the boundary values of Jx D1 ax / r.
"""
import sympy as sp

x, px, y, py, d, r, e = sp.symbols('x p_x y p_y delta rho epsilon', real=True)
X, Px, Y, Py = sp.symbols('X P_X Y P_Y', real=True)
D, Dp = sp.symbols('D1 Dp1', real=True)
Jx, Jy, bx, by, ax, ay = sp.symbols('J_x J_y beta_x beta_y alpha_x alpha_y', positive=True)
a, b = sp.symbols('psi_x psi_y', real=True)

H = -(1 + x/r)*sp.sqrt((1 + d)**2 - px**2 - py**2) + x/r + x**2/(2*r**2)

def series2(expr, order):
    return sp.expand(sp.series(expr.subs({x: e*x, px: e*px, y: e*y, py: e*py, d: e*d}), e, 0, order).removeO().subs(e, 1))

# 1. equations of motion
xp = sp.diff(H, px); pxp = -sp.diff(H, x)
xp1 = series2(xp, 2)
pxp1 = series2(pxp, 2)
print('x\'  (lin):', xp1)
print('px\' (lin):', pxp1)
# x' = px/(1+d) at linear order -> x'' = px'/(1+d) = d/(r(1+d)) - x/(r^2 (1+d)); dK/dd = -1/r^2

# 2. second-order force at d = 0, in geometric variables
xps, yps = sp.symbols("x' y'", real=True)
H0 = H.subs(d, 0)
xp_2 = series2(sp.diff(H0, px), 3)
pxp_2 = series2(-sp.diff(H0, x), 3)
yp_2 = series2(sp.diff(H0, py), 3)
pyp_2 = series2(-sp.diff(H0, y), 3)
# x'' = d/ds x'(x,px,y,py) along the flow
xpp = sum(sp.diff(xp_2, v)*f for v, f in [(x, xp_2), (px, pxp_2), (y, yp_2), (py, pyp_2)])
xpp = series2(xpp, 3)
ypp = sum(sp.diff(yp_2, v)*f for v, f in [(x, xp_2), (px, pxp_2), (y, yp_2), (py, pyp_2)])
ypp = series2(ypp, 3)
# invert px(x,x') to 2nd order: px = x'(1 - x/r), py = y'(1 - x/r)
sub_geo = {px: xps*(1 - x/r), py: yps*(1 - x/r)}
Fx = series2(xpp.subs(sub_geo) + x/r**2, 3).subs({px: 0})
Fx = sp.expand(series2(sp.expand(xpp.subs(sub_geo)), 3) + x/r**2)
Fx = sp.expand(sp.series(Fx.subs({x: e*x, xps: e*xps, yps: e*yps, y: e*y}), e, 0, 3).removeO().subs(e, 1))
Fy = sp.expand(sp.series(sp.expand(ypp.subs(sub_geo)).subs({x: e*x, xps: e*xps, yps: e*yps, y: e*y}), e, 0, 3).removeO().subs(e, 1))
print('F^x =', Fx)
print('y\'\' =', Fy)
assert sp.simplify(Fx - (-x**2/r**3 + (xps**2 - yps**2)/(2*r))) == 0
assert sp.simplify(Fy - xps*yps/r) == 0

# 3. chromaticity: expand around closed orbit x = D1 d + X, px = D1' d + Px
Hs = H.subs({x: D*d + e*X, px: Dp*d + e*Px, y: e*Y, py: e*Py})
Hs = sp.series(Hs, e, 0, 3).removeO()
H2 = sp.expand(Hs.coeff(e, 2))
H2_0 = sp.expand(H2.subs(d, 0))
H2_1 = sp.expand(sp.diff(H2, d).subs(d, 0))
print('H2(d=0) =', H2_0)
print('dH2/dd  =', sp.factor_terms(H2_1))
assert sp.simplify(H2_0 - (Px**2/2 + Py**2/2 + X**2/(2*r**2))) == 0

tw = {X: sp.sqrt(2*Jx*bx)*sp.cos(a), Px: -sp.sqrt(2*Jx/bx)*(sp.sin(a) + ax*sp.cos(a)),
      Y: sp.sqrt(2*Jy*by)*sp.cos(b), Py: -sp.sqrt(2*Jy/by)*(sp.sin(b) + ay*sp.cos(b))}
def avg(expr):
    f = sp.expand(sp.expand_trig(expr.subs(tw)))
    f = sp.integrate(f, (a, 0, 2*sp.pi))/(2*sp.pi)
    f = sp.integrate(f, (b, 0, 2*sp.pi))/(2*sp.pi)
    return sp.simplify(f)
gx = (1 + ax**2)/bx; gy = (1 + ay**2)/by
H21 = avg(H2_1)
Ix = sp.simplify(2*sp.diff(H21, Jx))
Iy = sp.simplify(2*sp.diff(H21, Jy))
print('I_x =', sp.simplify(Ix.subs(ax**2, gx*bx - 1)))
print('I_y =', sp.simplify(Iy.subs(ay**2, gy*by - 1)))

# 4. Lambda integrand inside dipole: <x'^2+y'^2>/2 + <F^x> D1
geo = {x: sp.sqrt(2*Jx*bx)*sp.cos(a), xps: -sp.sqrt(2*Jx/bx)*(sp.sin(a) + ax*sp.cos(a)),
       y: sp.sqrt(2*Jy*by)*sp.cos(b), yps: -sp.sqrt(2*Jy/by)*(sp.sin(b) + ay*sp.cos(b))}
def avg_geo(expr):
    f = sp.expand(sp.expand_trig(expr.subs(geo)))
    f = sp.integrate(f, (a, 0, 2*sp.pi))/(2*sp.pi)
    f = sp.integrate(f, (b, 0, 2*sp.pi))/(2*sp.pi)
    return sp.simplify(f)
Fav = avg_geo(Fx)
print('<F^x> =', sp.simplify(Fav))
Lint = Jx*gx/2 + Jy*gy/2 + Fav*D

# Twiss inside the dipole: K_x = 1/r^2, alpha_x' = K_x beta_x - gamma_x, beta_x' = -2 alpha_x
dax = bx/r**2 - gx
lhs = Lint + (Jx*Ix + Jy*Iy)/2
rhs = -(Jx/r)*(Dp*ax + D*dax)   # -d/ds (Jx D1 alpha_x / rho), rho constant inside
assert sp.simplify(lhs - rhs) == 0
print('inside dipole: Lambda-integrand + J*I/2 = -d/ds(Jx D1 ax/rho)')

# quadrupole + sextupole (no bend), same canonical recipe: Lambda-integrand + J*I/2 = 0
K1, K2 = sp.symbols('K1 K2', real=True)
Hq = -sp.sqrt((1 + d)**2 - px**2 - py**2) + K1*(x**2 - y**2)/2 + K2*(x**3 - 3*x*y**2)/6
Hqs = sp.series(Hq.subs({x: D*d + e*X, px: Dp*d + e*Px, y: e*Y, py: e*Py}), e, 0, 3).removeO()
Hq21 = avg(sp.expand(sp.diff(sp.expand(Hqs.coeff(e, 2)), d).subs(d, 0)))
Iqx = sp.simplify(2*sp.diff(Hq21, Jx)); Iqy = sp.simplify(2*sp.diff(Hq21, Jy))
Fq_av = -K2*(Jx*bx - Jy*by)/2
Lq = Jx*gx/2 + Jy*gy/2 + Fq_av*D
assert sp.simplify(Lq + (Jx*Iqx + Jy*Iqy)/2) == 0
print('quad/sextupole: I_x =', sp.simplify(Iqx.subs(ax**2, gx*bx - 1)), ' I_y =', sp.simplify(Iqy.subs(ay**2, gy*by - 1)))

# 5. slope jump at faces: entrance x'_after = x'_before (1 + x/r); exit x'_after = x'_before (1 - x/r)
jump_entr = avg_geo(x*xps/r)
print('<dx\'> at entrance =', jump_entr, ' at exit =', -jump_entr)
# boundary term in Green identity: sum D1 <dx'> = (Jx/r)(D1 ax|exit - D1 ax|entr)
assert sp.simplify(jump_entr + Jx*ax/r) == 0
print('all checks passed')
