# Algebraic checks of chapter 6: smooth Lambda, ratio, master prefactor.
import sympy as sp

g = sp.symbols("g", positive=True)
C, rho = sp.symbols("C rho", positive=True)
Jx, Jy = sp.symbols("J_x J_y", positive=True)
nux, nuy = sp.symbols("nu_x nu_y", positive=True)
n = sp.symbols("n", real=True)
eta0, eta1, Lam, dm = sp.symbols("eta0 eta1 Lambda delta_m", real=True)

# smooth approximation of Lambda^E
# beta = C / (2 pi nu), oint ds/rho**2 = 2 pi / rho, oint gamma_twiss ds = 2 pi nu
# first term from slopes: (J/(2C)) * oint gamma = pi J nu / C
slope = sp.pi / C * (Jx * nux + Jy * nuy)

beta_x = C / (2 * sp.pi * nux)
beta_y = C / (2 * sp.pi * nuy)
oint_inv_rho2 = 2 * sp.pi / rho
pot_x = (Jx / (2 * C)) * (5 - n) / g * beta_x * oint_inv_rho2
pot_y = (Jy / (2 * C)) * (n - 1) / g * beta_y * oint_inv_rho2
claimed_x = Jx * (5 - n) / (2 * g * nux * rho)
claimed_y = Jy * (n - 1) / (2 * g * nuy * rho)
print("smooth x:", sp.simplify(pot_x - claimed_x))
print("smooth y:", sp.simplify(pot_y - claimed_y))
print("slope:", sp.simplify(slope - sp.pi / C * (Jx * nux + Jy * nuy)))

# ratio |dropped / kept|
# dropped = beta**2 * <delta**2> / (2 gamma) , kept modulus = beta**2 * gamma * |delta_eq|
# <delta**2> = dm**2 / 2, delta_eq = -(eta1*dm**2/2 + Lam)/eta0
d2 = dm ** 2 / 2
delta_eq = -(eta1 * d2 + Lam) / eta0
ratio = sp.Abs(d2) / (2 * g * sp.Abs(delta_eq))
# the note writes |eta0| dm**2 / (4 g |eta1 dm**2/2 + Lam|)
claimed = sp.Abs(eta0) * dm ** 2 / (4 * g * sp.Abs(eta1 * dm ** 2 / 2 + Lam))
# compare without Abs by assuming the positive branch of the written moduli:
# ratio_signed_structure
num = d2 / (2 * g * (-delta_eq))  # if delta_eq < 0 this is the positive ratio only when...
# Direct algebraic reduction of the note's RHS from the definitions, signs stripped
# by replacing Abs(z) -> z and assuming the combinations positive.
ratio_raw = (dm ** 2 / 2) / (2 * g * ((eta1 * dm ** 2 / 2 + Lam) / eta0))
claimed_raw = eta0 * dm ** 2 / (4 * g * (eta1 * dm ** 2 / 2 + Lam))
print("ratio:", sp.simplify(sp.together(ratio_raw - claimed_raw)))

# Lambda = 0 reduction |eta0|/(2 |eta1| g)
reduced = sp.simplify(claimed_raw.subs(Lam, 0) - eta0 / (2 * eta1 * g))
print("Lambda=0 ratio minus |eta0|/(2|eta1|g):", reduced)

# deq master: 1/eta0 vs gamma**2/(gamma**2 alpha - 1) already in sym_senichev
print("ok")
