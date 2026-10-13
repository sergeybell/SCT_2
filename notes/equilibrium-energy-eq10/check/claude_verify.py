"""Independent supervisor checks (Claude). Run: python claude_verify.py"""
import sympy as sp

d, g, a0, a1 = sp.symbols('delta gamma alpha0 alpha1', positive=True)
b2 = 1 - 1/g**2

print("=== 1. Magnetic ring: tau(delta) from exact velocity ===")
# p = p_s(1+d); beta*gamma = bg_s (1+d); v_s/v = beta_s/beta
bg_s = sp.sqrt(g**2 - 1)
bg = bg_s*(1 + d)
beta = bg/sp.sqrt(1 + bg**2)
beta_s = bg_s/g
ratio = sp.series((beta_s/beta)*(1 + a0*d + a1*d**2), d, 0, 3).removeO()
eta0 = sp.simplify(ratio.coeff(d, 1))
eta1 = sp.simplify(ratio.coeff(d, 2))
print("eta0 =", eta0)
eta1_expected = a1 - a0/g**2 + 1/g**4 + sp.Rational(3, 2)*b2/g**2
print("eta1 - expected =", sp.simplify(eta1 - eta1_expected))
vs_v = sp.series(beta_s/beta, d, 0, 3).removeO()
print("v_s/v coeffs:", sp.simplify(vs_v.coeff(d, 1)), sp.simplify(vs_v.coeff(d, 2) - 1/g**4 - sp.Rational(3, 2)*b2/g**2))
senichev_bracket = a1 - a0/g**2 + 1/g**4
print("eta1 - Senichev bracket =", sp.simplify(eta1 - senichev_bracket))

print("\n=== 2. Ng (18.20), (18.21) ===")
e0, e1 = sp.symbols('eta0 eta1')
w = sp.series(1/(1 + e0*d + e1*d**2), d, 0, 3).removeO()  # omega/omega0
eta_1820 = sp.expand(-(w - 1)/d)
print("(18.20) eta1' =", sp.expand(eta_1820.coeff(d, 1)))
eta_1821 = sp.expand(-sp.diff(w, d))
print("(18.21) eta1'' =", sp.expand(eta_1821.coeff(d, 1)))
# Ng (18.18) with alpha1 = alpha0*alpha1Ng
a1ng = sp.symbols('alpha1Ng')
ng1818 = a0*a1ng + sp.Rational(3, 2)*b2/g**2 - (a0 - 1/g**2)/g**2
print("Ng(18.18) - ours =", sp.simplify(ng1818 - eta1_expected.subs(a1, a0*a1ng)))

print("\n=== 3. Lindstedt for x'' + w^2 x = c0 + c2 x^2 (2nd order) ===")
eps, t, A, W = sp.symbols('epsilon t A W', positive=True)
c0, c2, om = sp.symbols('c0 c2 omega', positive=True)
th = sp.symbols('theta')
# x = A cos th + eps x1 + eps^2 x2, Omega = om + eps w1 + eps^2 w2; c0->eps*eps*c0? treat c0 ~ eps*A^2 order
w1, w2 = sp.symbols('w1 w2')
X1 = sp.Function('X1')(th)
X2 = sp.Function('X2')(th)
x = A*sp.cos(th) + eps*X1 + eps**2*X2
Om = om + eps*w1 + eps**2*w2
expr = Om**2*sp.diff(x, th, 2) + om**2*x - eps*(c0 + c2*x**2)
ser = sp.expand(sp.series(expr, eps, 0, 3).removeO())
o1 = ser.coeff(eps, 1)
# order 1: om^2 (X1''+X1) = c0 + c2 A^2 cos^2 + 2 om w1 A cos
print("order-1 secular requires w1 = 0")
X1s = (c0 + c2*A**2/2)/om**2 - c2*A**2/(6*om**2)*sp.cos(2*th)
chk = sp.simplify(om**2*(sp.diff(X1s, th, 2) + X1s) - (c0 + c2*A**2*sp.cos(th)**2))
print("X1 check:", chk)
o2 = ser.coeff(eps, 2).subs(w1, 0).subs(X1, X1s).doit()
# o2 = om^2 (X2''+X2) + rest ; secular cos term of rest must vanish
rest = sp.expand(o2 - om**2*(sp.diff(X2, th, 2) + X2))
cos_coef = sp.integrate(sp.expand(rest*sp.cos(th)), (th, 0, 2*sp.pi))/sp.pi
w2sol = sp.solve(sp.simplify(cos_coef), w2)[0]
print("w2 =", sp.simplify(w2sol))
print("w2 - expected (-c0 c2/om^3 - 5/12 c2^2 A^2/om^3) =",
      sp.simplify(w2sol - (-c0*c2/om**3 - sp.Rational(5, 12)*c2**2*A**2/om**3)))
print("mean x (order eps) =", sp.simplify((c0 + c2*A**2/2)/om**2))

print("\n=== 4. Smooth electrostatic rings, radial field E ~ r^-k ===")
# units m=c=rho=1; reference orbit r=1 with gamma=g; qE0 = p_s v_s
r1, r2, G1, G2 = sp.symbols('r1 r2 G1 G2')
rr = 1 + r1*d + r2*d**2
GG = g + G1*d + G2*d**2
qE0 = bg_s*beta_s
E_out = sp.series(sp.sqrt(1 + bg**2), d, 0, 3).removeO()  # total energy, defined at Phi=0
for kk in [1, 2, sp.Rational(3, 2)]:
    qPhi = qE0*sp.log(rr) if kk == 1 else qE0*(rr**(1 - kk) - 1)/(1 - kk)
    eq_orbit = sp.series(GG**2 - 1 - GG*qE0*rr**(1 - kk), d, 0, 3).removeO()  # p v = qE r
    eq_energy = sp.series(GG + qPhi - E_out, d, 0, 3).removeO()
    sol = sp.solve([eq_orbit.coeff(d, n) for n in (1, 2)] + [eq_energy.coeff(d, n) for n in (1, 2)],
                   [r1, r2, G1, G2], dict=True)[0]
    v_loc = sp.sqrt(1 - 1/GG**2)
    tau = sp.series((rr/v_loc).subs(sol)*beta_s - 1, d, 0, 3).removeO()
    al0, al1 = sp.simplify(sol[r1]), sp.simplify(sol[r2])
    e0s, e1s = sp.simplify(tau.coeff(d, 1)), sp.simplify(tau.coeff(d, 2))
    print(f"k={kk}: alpha0={al0}, alpha1={al1}")
    print(f"   eta0={e0s};  eta0 - [alpha0(1+1/g^2)-1/g^2] =", sp.simplify(e0s - (al0*(1 + 1/g**2) - 1/g**2)))
    print(f"   eta1={sp.factor(e1s)}")
    print("   eta1 - eta1_mag(alpha0,alpha1) =", sp.factor(e1s - eta1_expected.subs({a0: al0, a1: al1})))
