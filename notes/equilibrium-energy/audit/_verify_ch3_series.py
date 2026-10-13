"""Superseded by _verify_ch3_eta.py.

The expansion below treats beta2 and gamma as independent symbols and breaks
the cos/sin substitution, so its printed series is not used in the chapter.

Analytic second-order time and path for a smooth cylindrical deflector.

Independent variable is the laboratory angle. The radial deviation is expanded
through second order in the amplitude and the revolution time is averaged.
"""

import sympy as sp

rho, eps, psi = sp.symbols("rho eps psi", positive=True)
beta2, gam = sp.symbols("beta2 gamma", positive=True)
nu = sp.sqrt(2 - beta2)
c0, c2 = sp.symbols("c0 c2")

# xi = x/rho. psi = nu * theta, so d/dtheta = nu d/dpsi.
xi = c0 * eps**2 + eps * sp.cos(psi) + c2 * eps**2 * sp.cos(2 * psi)
dxi_dpsi = sp.diff(xi, psi)
dxi_dth = nu * dxi_dpsi

E = gam * (1 - beta2 * sp.log(1 + xi))
E_tp = gam * (1 - beta2 * sp.log(1 + eps))  # outer turning point, xi=eps
L = rho * (1 + eps) * sp.sqrt(E_tp**2 - 1)
r = rho * (1 + xi)

lhs = (rho * dxi_dth) ** 2
rhs = r**4 / L**2 * (E**2 - 1) - r**2
resid = sp.series(sp.expand(lhs - rhs), eps, 0, 3).removeO()
resid = sp.expand(resid)

# Express everything as a polynomial in z = cos(psi).
z = sp.symbols("z")
resid_z = resid.subs(sp.cos(2 * psi), 2 * z**2 - 1).subs(sp.sin(2 * psi), 0)
# sin(psi)**2 and cos(psi) derivatives already became cos/sin.
# dxi has sin(psi). lhs has sin**2. Replace sin**2 = 1-z**2, and any remaining sin.
resid_z = sp.expand(resid_z.rewrite(sp.sin))
resid_z = resid_z.subs(sp.sin(psi) ** 2, 1 - z**2).subs(sp.sin(psi), 0)
resid_z = resid_z.subs(sp.cos(psi), z)
resid_z = sp.expand(resid_z)

# Kill odd powers? Collect z**0, z**1, z**2 at order eps**2, and z**0? at eps**1.
poly = sp.Poly(sp.collect(sp.expand(resid_z), eps), eps)
# Order eps**1 must vanish (linear tune already built in) and eps**2 coefficients of z^0 and z^2.
e1 = sp.expand(resid_z.coeff(eps, 1))
e2 = sp.expand(resid_z.coeff(eps, 2))
print("order eps^1:", sp.factor(e1))
# e2 is quadratic in z. Coefficients of 1, z, z**2.
e2p = sp.Poly(sp.expand(e2), z)
print("e2 degrees", e2p.degree())
eqs = [sp.simplify(e2p.nth(k)) for k in range(e2p.degree() + 1)]
for k, eq in enumerate(eqs):
    print(f"z^{k}:", sp.factor(eq))

sols = sp.solve(eqs, [c0, c2], dict=True)
print("solutions", len(sols))
sol = sols[0]
print("c0 =", sp.factor(sol[c0]))
print("c2 =", sp.factor(sol[c2]))

# Time: dt/dtheta = E * r**2 / L
dt = sp.series(sp.expand(E * r**2 / L), eps, 0, 3).removeO()
dt = sp.expand(dt)
# reference dt/dtheta = rho / beta, beta=sqrt(beta2), E=gam, r=rho, L=rho*beta*gam
# T/Tref - 1 averaged. dt_ref = rho/sp.sqrt(beta2)
dt_ref = rho / sp.sqrt(beta2)
rel = sp.series(sp.expand(dt / dt_ref - 1), eps, 0, 3).removeO()
rel = sp.expand(rel.subs(sol))
rel_z = rel.subs(sp.cos(2 * psi), 2 * z**2 - 1).subs(sp.cos(psi), z)
rel_z = sp.expand(rel_z.rewrite(sp.sin).subs(sp.sin(psi) ** 2, 1 - z**2).subs(sp.sin(psi), 0).subs(sp.cos(psi), z))
rel_avg = sp.expand(rel_z).subs(z**2, sp.Rational(1, 2)).subs(z, 0)
# The substitution z^2 -> 1/2 then z->0 is the average only if we already reduced powers.
# Do it properly: poly in z, average z^2=1/2, z^4=3/8, odd=0.
rp = sp.Poly(sp.expand(rel_z), z)
avg = 0
for k in range(rp.degree() + 1):
    if k % 2 == 0:
        # average of z^k = (k-1)!!/(k!!) 
        avg += rp.nth(k) * sp.binomial(k, k // 2) / (2**k)
avg = sp.simplify(sp.together(avg))
print("mean (dt-dtref)/dtref =", avg)

# geometric action: x = rho*xi, x' = dx/ds = dxi/dtheta since s=rho*theta, dx/ds=dxi/dth
# J = A^2 / (2 beta_twiss), beta_twiss = rho/nu, A = eps*rho
# so J = eps**2 * rho**2 / (2 * rho/nu) = eps**2 * rho * nu / 2
# (dT/T)/J = avg / J
J = eps**2 * rho * nu / 2
kT = sp.simplify(sp.expand(avg / J))
print(" (dT/T)/J =", kT)

# Path: dl/dtheta = r * sqrt(1 + (dr/(r dtheta))**2) = r * sqrt(1 + (rho dxi/dth / r)**2)
dl = r * sp.sqrt(1 + (rho * dxi_dth / r) ** 2)
dl = sp.series(sp.expand(dl), eps, 0, 3).removeO()
path_rel = sp.series(sp.expand(dl / rho - 1), eps, 0, 3).removeO()
path_rel = sp.expand(path_rel.subs(sol))
path_z = path_rel.subs(sp.cos(2 * psi), 2 * z**2 - 1)
path_z = sp.expand(path_z.rewrite(sp.sin).subs(sp.sin(psi) ** 2, 1 - z**2).subs(sp.sin(psi), 0).subs(sp.cos(psi), z))
pp = sp.Poly(sp.expand(path_z), z)
pavg = 0
for k in range(pp.degree() + 1):
    if k % 2 == 0:
        pavg += pp.nth(k) * sp.binomial(k, k // 2) / (2**k)
pavg = sp.simplify(pavg)
kC = sp.simplify(sp.expand(pavg / J))
print(" (dC/C)/J =", kC)
print("c0 numeric-ready", sol[c0])
