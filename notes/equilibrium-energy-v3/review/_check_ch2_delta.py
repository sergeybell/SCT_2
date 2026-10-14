# Fixed-number checks for chapter 2. No fits, no scans.
import math

# ---------- 1. path element: (x/rho)^2 cancels through second order ----------
x_rho = 0.02
xp = 0.03
yp = -0.04
exact = math.sqrt((1 + x_rho) ** 2 + xp ** 2 + yp ** 2)
approx = 1 + x_rho + 0.5 * (xp ** 2 + yp ** 2)
# next geometric piece is third order: (x/rho)*(xp**2+yp**2)/2
third = 0.5 * x_rho * (xp ** 2 + yp ** 2)
print("path exact", exact)
print("path approx", approx)
print("path diff", exact - approx)
print("third-order piece", third)

# ---------- 2. <x'^2> = J * gamma ----------
beta = 2.5
alpha = 0.4
J = 1.7e-6
gamma = (1 + alpha ** 2) / beta
# exact phase average, 4-point Gauss on [0, 2pi] is unnecessary:
# <(sin + alpha cos)^2> = (1+alpha^2)/2
avg_factor = (1 + alpha ** 2) / 2
xp2 = (2 * J / beta) * avg_factor
print("xp2", xp2, "J*gamma", J * gamma, "diff", xp2 - J * gamma)

# ---------- 3. uniform focusing channel ----------
K = 0.25
C = 40.0
nu = C * math.sqrt(K) / (2 * math.pi)
d = 1.0e-4
nu_p = C * math.sqrt(K / (1 + d)) / (2 * math.pi)
nu_m = C * math.sqrt(K / (1 - d)) / (2 * math.pi)
xi_fd = (nu_p - nu_m) / (2 * d)
xi_exact = -nu / 2
Jch = 2.5e-6
Lam_g = (Jch / 2) * (C / (1 / math.sqrt(K))) / C   # (J/2)*oint gamma ds / C, gamma=sqrt(K), beta=1/sqrt(K)
# gamma = 1/beta = sqrt(K), oint gamma ds = C*sqrt(K) = 2*pi*nu
Lam_g = (Jch / 2) * (2 * math.pi * nu) / C
Lam_xi = -(2 * math.pi / C) * xi_exact * Jch
print("nu", nu, "xi_fd", xi_fd, "xi_exact", xi_exact)
print("Lambda_g", Lam_g, "Lambda_xi", Lam_xi, "diff", Lam_g - Lam_xi)

# ---------- 4. FODO: oint gamma = oint beta K = -4 pi xi, not 2 pi nu ----------
L = 2.0
q = 0.4  # 1/m, thin quad
# sin(mu/2) = L*q/2
s2 = L * q / 2
mu = 2 * math.asin(s2)
nu_cell = mu / (2 * math.pi)
# beta at thin quads (mid-lens formula for thin FODO)
bF = 2 * L * (1 + s2) / math.sin(mu)
bD = 2 * L * (1 - s2) / math.sin(mu)
xi = -q * (bF - bD) / (4 * math.pi)
# integrate gamma through the two drifts by propagating Twiss
# Start just after QF: beta=bF, alpha = alpha_after
# Thin lens QF: alpha_+ = alpha_- + beta * q
# Periodic: the drift-QD-drift-QF closes.
# Propagate with matrices from the periodic solution of the cell matrix.

def mat_drift(L):
    return ((1.0, L), (0.0, 1.0))

def mat_quad(k):
    # x' := x' - k x   (k>0 focusing)
    return ((1.0, 0.0), (-k, 1.0))

def mm(A, B):
    return (
        (A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
        (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]),
    )

# cell: QF/2, drift L, QD, drift L, QF/2  (symmetric) OR QF, drift, QD, drift
# Use QF then drift L then QD then drift L, Twiss at entrance of QF.
M = mm(mat_drift(L), mm(mat_quad(-q), mm(mat_drift(L), mat_quad(q))))
# one-turn (cell) matrix M = [[m11,m12],[m21,m22]]
m11, m12 = M[0]
m21, m22 = M[1]
cosmu = (m11 + m22) / 2
# beta, alpha at start: beta = m12/sin(mu), alpha = (m11-m22)/(2 sin)
sinmu = math.sqrt(max(0.0, 1 - cosmu ** 2))
# sign of sin: mu in (0, pi) for this FODO, sin>0
beta0 = m12 / sinmu
alpha0 = (m11 - m22) / (2 * sinmu)
print("FODO nu", math.acos(cosmu) / (2 * math.pi), "analytic", nu_cell)
print("beta0", beta0, "bF", bF, "alpha0", alpha0)

def propagate(beta, alpha, M):
    # Twiss transport: (beta, alpha, gamma) by matrix
    a, b = M[0]
    c, d = M[1]
    g = (1 + alpha ** 2) / beta
    b2 = a * a * beta - 2 * a * b * alpha + b * b * g
    a2 = -a * c * beta + (a * d + b * c) * alpha - b * d * g
    return b2, a2

# slice each drift into N steps and sum gamma*ds; quads are thin so oint beta*K = beta*q
N = 2000
ds = L / N
Md = mat_drift(ds)
# start at entrance, apply QF first (matches M construction: QF, drift, QD, drift)
beta, alpha = beta0, alpha0
beta, alpha = propagate(beta, alpha, mat_quad(q))
sum_g = 0.0
sum_bK = beta * q  # QF, K integrated = +q, but beta is AFTER the thin lens; beta continuous
# beta unchanged by thin lens, so beta*q uses beta at the quad. We already propagated alpha only.
# Recompute beta at quad: thin quad doesn't change beta. beta after = beta before = beta0 before quad?
# Order in M: first mat_quad(q) is applied to the vector on the right, so the cell STARTS with QF.
# I set beta,alpha = periodic at the entrance BEFORE QF, then propagated through QF.
# Beta is continuous, so beta at QF is beta0. But I added beta*q AFTER propagation.
# propagate through quad changes alpha, not beta. So beta at QF is still beta0.
# I incorrectly used post-quad beta. Fix: use beta0.
sum_bK = beta0 * q

def accumulate_drift(beta, alpha, L, N):
    ds = L / N
    Md = mat_drift(ds)
    sg = 0.0
    for _ in range(N):
        g = (1 + alpha ** 2) / beta
        sg += g * ds
        beta, alpha = propagate(beta, alpha, Md)
    return beta, alpha, sg

beta, alpha, sg1 = accumulate_drift(beta, alpha, L, N)
# now at QD
g_here = (1 + alpha ** 2) / beta
sum_bK += beta * (-q)
beta, alpha = propagate(beta, alpha, mat_quad(-q))
beta, alpha, sg2 = accumulate_drift(beta, alpha, L, N)
sum_g = sg1 + sg2
# close check: beta, alpha should return to pre-QF periodic values
print("close beta", beta, beta0, "alpha", alpha, alpha0)
print("oint gamma", sum_g, "oint beta K", sum_bK, "diff", sum_g - sum_bK)
print("-4 pi xi", -4 * math.pi * xi, "xi", xi)
print("2 pi nu", 2 * math.pi * nu_cell)
# (1/2) oint gamma vs -2 pi xi vs pi nu
print("half oint gamma", 0.5 * sum_g, "-2pi xi", -2 * math.pi * xi, "pi nu", math.pi * nu_cell)

# ---------- 5. closed orbit: oint x/rho ds = D * theta ----------
# Uniform channel discretized, plus 1/rho constant, plus a thin kick.
N = 400
C = 40.0
K = 0.25
rho = 10.0
ds = C / N
# state vector (x, xp), map of one step: drift-like with focusing
# x'' = -K x + 1/rho   for dispersion, homogeneous for betatron
# step matrix for x'' + K x = 0
om = math.sqrt(K)
# exact drift-focus step of length ds
cd = math.cos(om * ds)
sd = math.sin(om * ds)
# particular: constant D = 1/(K*rho), D'=0 is exact. Check numerically anyway.
# For a kick theta at index 0: jump xp by theta, then find periodic orbit.
# Solve (M_turn - I) (x,xp) = kick inhomogeneity.
# Build one-turn of the homogeneous map and the kick as xp -> xp + theta at s=0,
# periodicity: transport around the ring brings (x, xp+theta) back to (x, xp).
# Equivalent: (M - I) X = - M_from_after_kick_to_end applied to (0, theta) ... 
# Simpler shooting: transport full turn starting AFTER adding theta to xp,
# require end == start_before_kick.
# Let X = (x, xp) just before the kick.
# After kick: (x, xp+theta)
# After one turn of focusing (no more kicks): M X_after = X
# M (x, xp+theta) = (x, xp)

def turn_matrix(nsteps, ds, om):
    cd = math.cos(om * ds)
    sd = math.sin(om * ds)
    # one step [[cd, sd/om], [-om*sd, cd]]
    A = ((1.0, 0.0), (0.0, 1.0))
    step = ((cd, sd / om), (-om * sd, cd))
    for _ in range(nsteps):
        A = mm(step, A)
    return A

Mturn = turn_matrix(N, ds, om)
# M (x, xp+theta) = (x, xp)
# [[m11, m12],[m21, m22]] (x, xp+theta) = (x, xp)
theta = 1.0e-4
m11, m12 = Mturn[0]
m21, m22 = Mturn[1]
# m11 x + m12 (xp+theta) = x
# m21 x + m22 (xp+theta) = xp
# (m11-1) x + m12 xp = -m12 theta
# m21 x + (m22-1) xp = -m22 theta
A11, A12 = m11 - 1, m12
A21, A22 = m21, m22 - 1
b1 = -m12 * theta
b2 = -m22 * theta
det = A11 * A22 - A12 * A21
x0 = (b1 * A22 - A12 * b2) / det
xp0 = (A11 * b2 - b1 * A21) / det
# integrate x(s) and compare oint x/rho with D*theta
D = 1 / (K * rho)
# sample x along the ring starting after the kick
x, xp = x0, xp0 + theta
acc = 0.0
step = ((cd, sd / om), (-om * sd, cd))
for i in range(N):
    acc += x * ds / rho
    nx = step[0][0] * x + step[0][1] * xp
    nxp = step[1][0] * x + step[1][1] * xp
    x, xp = nx, nxp
print("orbit close", x, x0, xp, xp0)
print("oint x/rho", acc, "D*theta", D * theta, "diff", acc - D * theta)

# ---------- 6. sextupole feed-down: Delta L = -2 pi xi J ----------
# theta_avg = -(m/2) * J * beta, xi = beta * m * D / (4 pi)
# Delta L = D * theta
m_sext = 0.3  # 1/m^2
J = 1.7e-6
beta_u = 1 / math.sqrt(K)
theta_avg = -(m_sext / 2) * J * beta_u
dL = D * theta_avg
xi_s = beta_u * m_sext * D / (4 * math.pi)
dL_xi = -2 * math.pi * xi_s * J
print("sext dL", dL, "from xi", dL_xi, "diff", dL - dL_xi)
print("beta", beta_u, "D", D, "xi_s", xi_s)

# thin-lens tune shift check: extra focusing q at one point
# nu from trace of (uniform turn * thin lens)
q_extra = 0.01  # 1/m
# M_new = Mturn * quad(q)   (kick at the start)
Mq = mm(Mturn, mat_quad(q_extra))
cos1 = (Mq[0][0] + Mq[1][1]) / 2
# unperturbed
cos0 = (Mturn[0][0] + Mturn[1][1]) / 2
# nu = acos(cos)/(2pi), derivative: careful near the branch
nu0 = math.acos(cos0) / (2 * math.pi)
nu1 = math.acos(max(-1, min(1, cos1))) / (2 * math.pi)
# beta at the lens from Mturn: beta = m12/sin
sin0 = math.sin(2 * math.pi * nu0)
beta_at = Mturn[0][1] / sin0
dnu = nu1 - nu0
dnu_formula = beta_at * q_extra / (4 * math.pi)
print("tune shift", dnu, "formula", dnu_formula, "beta_at", beta_at, "beta_u", beta_u)

# ---------- 7. tau product: coefficient of Lambda is 1; Lambda*delta is the cross ----------
# symbolic-free rational check with fractions
from fractions import Fraction
gs2 = Fraction(9)  # gamma_s^2 = 9, gamma_s=3
a0 = Fraction(1, 5)
a1 = Fraction(1, 7)
# v_s/v coefficients from chapter 1
# eta0 = a0 - 1/gs2
# eta1 = a1 - a0/gs2 + (3*gs2 - 1)/(2*gs2**2)
eta0 = a0 - 1 / gs2
eta1 = a1 - a0 / gs2 + (3 * gs2 - 1) / (2 * gs2 ** 2)
# multiply (1 + a0 d + a1 d^2 + lam) * (1 - d/gs2 + c2 d^2) and read coef of lam and lam*d
c2 = (3 * gs2 - 1) / (2 * gs2 ** 2)
# coef of lam: 1
# coef of lam*d: -1/gs2
print("eta0", eta0, "eta1", eta1)
print("cross Lambda*delta coef", -1 / gs2)

# ---------- 8. ODE residual of the BM solution through order a^2 and Lambda ----------
# Substitute delta = d0 + a cos + b cos2 into the ODE and compare Fourier modes
# kept orders: const, cos, cos2. Products a^3, a*d0 if d0 is O(a^2) are O(a^3).
import sympy as sp

a, w, eta0s, eta1s, g2, Lam, t = sp.symbols("a omega eta0 eta1 g2 Lambda t", positive=True)
b = a ** 2 / 6 * (eta1s / eta0s - 2 / g2)
d0 = -Lam / eta0s - eta1s * a ** 2 / (2 * eta0s)
psi = w * t
delta = d0 + a * sp.cos(psi) + b * sp.cos(2 * psi)
dd = sp.diff(delta, t)
ddd = sp.diff(dd, t)
rhs = -(w ** 2 / eta0s) * (eta1s - eta0s / g2) * delta ** 2 - (1 / g2) * dd ** 2 - (w ** 2 / eta0s) * Lam
resid = sp.expand(ddd + w ** 2 * delta - rhs)
# drop O(a^3), O(Lam*a), O(Lam**2), O(a**4)
resid2 = resid.subs({a ** 3: 0})
# better: series in a and Lam
resid_s = sp.series(resid.subs(Lam, sp.symbols("Lsmall")), a, 0, 3).removeO()
# treat Lam as order a**2: set Lam = lam2 * a**2, expand to a**2
lam2 = sp.symbols("lam2")
resid_a = sp.expand(resid.subs(Lam, lam2 * a ** 2))
resid_a = sum(term for term in sp.Add.make_args(sp.expand(resid_a)) if sp.degree(term, a) <= 2)
# collect cos/sin
resid_a = sp.expand(resid_a)
# replace powers of cos, sin by multiple angles and drop nothing within a^2
resid_trig = sp.trigsimp(sp.expand_trig(resid_a))
# Fourier: integrate against 1, cos, sin, cos2, sin2
def mode(expr, harm):
    per = 2 * sp.pi / w
    if harm == 0:
        return sp.simplify(sp.integrate(expr, (t, 0, per)) / per)
    c = sp.simplify(sp.integrate(expr * sp.cos(harm * w * t), (t, 0, per)) * 2 / per)
    s = sp.simplify(sp.integrate(expr * sp.sin(harm * w * t), (t, 0, per)) * 2 / per)
    return c, s

# faster numerical Fourier with fixed numbers
subs = {a: sp.Rational(1, 500), w: 1, eta0s: sp.Rational(1, 50), eta1s: sp.Rational(1, 20), g2: 9, Lam: sp.Rational(1, 10**6)}
# At these numbers a^2 = 4e-6, Lam=1e-6, eta0=0.02, so Lam/eta0=5e-5, a=0.002
# Residual should be O(a^3)~8e-9 plus Lam*a ~ 2e-9
rn = sp.lambdify(t, sp.simplify(resid.subs(subs)), "numpy")
import numpy as np
ts = np.linspace(0, 2 * np.pi, 4001, endpoint=False)
rv = rn(ts)
c0 = np.mean(rv)
c1 = 2 * np.mean(rv * np.cos(ts))
s1 = 2 * np.mean(rv * np.sin(ts))
c2 = 2 * np.mean(rv * np.cos(2 * ts))
s2 = 2 * np.mean(rv * np.sin(2 * ts))
c3 = 2 * np.mean(rv * np.cos(3 * ts))
print("resid modes", c0, c1, s1, c2, s2, c3, "rms", np.sqrt(np.mean(rv ** 2)))

# analytic constant and cos2 of residual should vanish at kept order.
# Compute residual with symbols truncated.
delta_tr = d0 + a * sp.cos(psi) + b * sp.cos(2 * psi)
# Expand rhs keeping a^2, Lam, and NOT a^2*cos * etc that generate a^3 or Lam*a
# Do polynomial truncation via series in a formal eps
eps = sp.symbols("eps")
delta_e = delta_tr.subs({a: eps * a, Lam: eps ** 2 * Lam})
dd_e = sp.diff(delta_e, t)
ddd_e = sp.diff(dd_e, t)
rhs_e = -(w ** 2 / eta0s) * (eta1s - eta0s / g2) * delta_e ** 2 - (1 / g2) * dd_e ** 2 - (w ** 2 / eta0s) * (eps ** 2 * Lam)
resid_e = sp.expand(ddd_e + w ** 2 * delta_e - rhs_e)
# coefficient of eps^2 (order a^2 and Lam). eps^1 should be the linear oscillator, already zero.
def coef(expr, n):
    return sp.expand(expr).coeff(eps, n)

r1 = sp.simplify(coef(resid_e, 1))
r2 = sp.simplify(coef(resid_e, 2))
print("order eps^1", r1)
print("order eps^2", r2)
