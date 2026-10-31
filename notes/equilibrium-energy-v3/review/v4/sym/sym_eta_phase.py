"""Symbolic identities for eta1 conventions, ch.1 kinematics, ch.3 phase series.

No numerics: series, expand, simplify, cancel only.
"""
import sympy as sp

g, a0, a1, d = sp.symbols("gamma_s alpha_0 alpha_1 delta", positive=True)
b2 = 1 - 1 / g**2  # beta_s**2

# --- kinematic series of the note -----------------------------------------
# gamma**2 / gamma_s**2 exact
g2_ratio = (1 + (g**2 - 1) * (1 + d) ** 2) / g**2
g2_ratio_s = sp.expand(g2_ratio)
# claimed exact polynomial
claimed_g2 = 1 + 2 * b2 * d + b2 * d**2
print("gamma2 ratio minus claimed:", sp.simplify(g2_ratio_s - claimed_g2))

x = 2 * b2 * d + b2 * d**2
g_over = sp.series(sp.sqrt(1 + x), d, 0, 3).removeO()
claimed_g_over = 1 + b2 * d + (g**2 - 1) / (2 * g**4) * d**2
print("gamma/gamma_s coeff check:", sp.simplify(sp.expand(g_over - claimed_g_over)))

vsv = sp.series(g_over / (1 + d), d, 0, 3).removeO()
claimed_vsv = 1 - d / g**2 + (3 * g**2 - 1) / (2 * g**4) * d**2
print("v_s/v check:", sp.simplify(sp.expand(vsv - claimed_vsv)))

# Senichev bracket vs note eta1
eta1 = a1 - a0 / g**2 + (3 * g**2 - 1) / (2 * g**4)
sen = a1 + 3 * b2 / (2 * g**2) - a0 / g**2 + 1 / g**4
print("eta1 - Senichev:", sp.simplify(sp.expand(eta1 - sen)))
print("kinematic piece as 3 beta^2/(2 gamma^2)+1/gamma^4:",
      sp.simplify(sp.expand((3 * g**2 - 1) / (2 * g**4) - (3 * b2 / (2 * g**2) + 1 / g**4))))

# incomplete Senichev 2013 bracket (no 3 beta^2 term)
sen_incomplete = a1 - a0 / g**2 + 1 / g**4
print("eta1 - incomplete:", sp.factor(sp.simplify(eta1 - sen_incomplete)))

# period -> frequency
eta0 = a0 - 1 / g**2
tau = eta0 * d + eta1 * d**2
omega_ratio = sp.series(1 / (1 + tau), d, 0, 3).removeO()
# omega/omega_s = 1 - eta0*d - eta1_F * d**2
# collect
c1 = sp.simplify(omega_ratio.coeff(d))
c2 = sp.simplify(omega_ratio.coeff(d**2))
eta1_F = sp.simplify(-c2)  # from 1 - eta0 d - eta1_F d**2
print("freq linear coeff (expect -eta0):", sp.simplify(c1 + eta0))
print("eta1_F - (eta1 - eta0**2):", sp.simplify(sp.expand(eta1_F - (eta1 - eta0**2))))

# derivative definition eta = -(1/omega_s) d omega / d delta = eta0 + eta1_D * delta + ...
# omega/omega_s series, derivative
domega = sp.diff(omega_ratio, d)
eta_deriv = sp.series(-domega, d, 0, 2).removeO()
eta1_D = sp.simplify(eta_deriv.coeff(d))
print("eta1_D - 2(eta1-eta0**2):", sp.simplify(sp.expand(eta1_D - 2 * (eta1 - eta0**2))))

# (1/omega) d omega / d delta, omega = omega_s * omega_ratio
# d ln omega / d delta = (1/omega_ratio) * d(omega_ratio)/d delta
dln = sp.series(sp.diff(omega_ratio, d) / omega_ratio, d, 0, 2).removeO()
print("-(d ln omega / d delta) series:", sp.series(-dln, d, 0, 2).removeO())
eta1_ln = sp.simplify(sp.series(-dln, d, 0, 2).removeO().coeff(d))
print("eta1_ln - 2(eta1 - eta0**2/2):", sp.simplify(sp.expand(eta1_ln - 2 * (eta1 - eta0**2 / 2))))

# gamma expansion
g_part = sp.series(g * sp.sqrt(1 + x), d, 0, 3).removeO()
claimed_g = g + b2 * g * d + b2 / (2 * g) * d**2
print("gamma particle check:", sp.simplify(sp.expand(g_part - claimed_g)))

# ratio of energy pieces
# |beta^2 dm^2/(4 gamma)| / |beta^2 gamma delta_eq|
# delta_eq = -eta1 dm**2 / (2 eta0), absolute ratio |eta0|/(2|eta1| gamma**2)
dm, eta0s, eta1s = sp.symbols("delta_m eta_0 eta_1", real=True)
deq = -eta1s * dm**2 / (2 * eta0s)
kin = b2 * dm**2 / (4 * g)
orb = b2 * g * deq
ratio = sp.simplify(kin / (-orb))  # both positive if eta1/eta0 > 0? keep symbolic
print("ratio kin/|orbit| structure:", sp.simplify(sp.together(kin / sp.Abs(orb))))
# algebraic without abs: kin / (beta^2 gamma |deq|) 
print("ratio simplified:", sp.factor(sp.simplify(kin / (b2 * g * sp.Abs(deq)))))

# --- ch.3 cubic BM and half-amplitude --------------------------------------
phi, ws, a, psi = sp.symbols("phihat omega_s a psi", positive=True)
# a = phi + phi**3/192
a_of = phi + phi**3 / 192
# Omega/ws = 1 - a**2/16 -> 1 - phi**2/16 + O(phi**4)
Om = sp.series(1 - a_of**2 / 16, phi, 0, 4).removeO()
print("Omega/ws:", Om)

# delta_m factor (1 - phi**2/16)*(1 + phi**2/48)
dm_fac = sp.series((1 - phi**2 / 16) * (1 + phi**2 / 48), phi, 0, 4).removeO()
print("dm factor:", dm_fac)

cos1 = sp.series((1 + phi**2 / 192) * (1 - phi**2 / 48), phi, 0, 4).removeO()
print("cos1 coeff:", cos1)

# elliptic
k = sp.symbols("k", positive=True)
# use known series of E, K through k**4 and divide
Kser = sp.pi / 2 * (1 + k**2 / 4 + sp.Rational(9, 64) * k**4)
Eser = sp.pi / 2 * (1 - k**2 / 4 - sp.Rational(3, 64) * k**4)
ratio_EK = sp.series(Eser / Kser, k, 0, 6).removeO()
print("E/K:", ratio_EK)
frac = sp.series((ratio_EK - 1 + k**2) / k**2, k, 0, 4).removeO()
print("elliptic fraction:", frac)
k2 = sp.series(sp.sin(phi / 2) ** 2, phi, 0, 6).removeO()
frac_phi = sp.series(frac.subs(k**2, k2).subs(k**4, k2**2), phi, 0, 6).removeO()
# safer substitution
frac_phi = sp.series(frac.subs(k, sp.sin(phi / 2)), phi, 0, 6).removeO()
print("<d^2>/dm^2 from elliptic:", frac_phi)
target = sp.Rational(1, 2) - phi**2 / 64
print("matches d2-ph:", sp.series(frac_phi - target, phi, 0, 6).removeO())

# exact pendulum frequency vs cubic: K(k)=pi/2(1+(1/2)^2 k^2+...)
# Omega/ws = pi/(2K) = 1/(1+k^2/4+...) 
Kfull = sp.elliptic_k(k**2)
Om_ex = sp.series(sp.pi / (2 * Kfull), k, 0, 4).removeO()
print("exact pendulum Omega/ws:", Om_ex)
Om_ex_phi = sp.series(Om_ex.subs(k, sp.sin(phi / 2)), phi, 0, 4).removeO()
print("exact pendulum in phihat:", Om_ex_phi)

# phase-shift identity cos(3*(psi-pi/2))
psi = sp.symbols("psi", real=True)
lhs = sp.cos(3 * (psi - sp.pi / 2))
rhs = -sp.sin(3 * psi)
print("phase shift 3rd harmonic:", sp.simplify(lhs - rhs))
