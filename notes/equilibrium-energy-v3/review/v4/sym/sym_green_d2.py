# Green identity for D2 in the electrostatic ring and the sextupole piece of eta1.
# Integrals are inert symbols: this is an algebraic identity, not a lattice sum.
import sympy as sp

g, C, K2 = sp.symbols("g C K2", real=True)
# I(f) means oint f ds. Linearity is encoded by treating independent integrals
# as basis symbols and imposing the Green relation as an equation.
I_S_D1, I_D2_over_rho = sp.symbols("I_S_D1 I_D2_over_rho")
factor = 1 + 1 / g

# Sextupole source in the D2 equation, outside the deflector:
# D2'' + K D2 = S2,   S2 = -K2/2 * D1**2,
# so oint S2 D1 ds = -K2/2 * oint D1**3 ds.
I_D1_cubed = sp.symbols("I_D1_cubed")
I_S_sext = -K2 / 2 * I_D1_cubed

# Correct Green identity: factor * oint D2/rho ds = oint S2 D1 ds
green = sp.Eq(factor * I_D2_over_rho, I_S_sext)
I_D2 = sp.solve(green, I_D2_over_rho)[0]

# eta1 depends on D2 only through factor * <D2/rho> = factor * I_D2 / C
delta_eta = sp.simplify(factor * I_D2 / C)
magnetic = sp.simplify(-K2 * I_D1_cubed / (2 * C))
print("ES sextupole shift:", delta_eta)
print("magnetic shift:", magnetic)
print("difference:", sp.simplify(delta_eta - magnetic))

# What goes wrong if the magnetic Green identity is copied without the factor:
wrong_I_D2 = I_S_sext  # oint D2/rho = oint S D1, factor forgotten
wrong_eta = sp.simplify(factor * wrong_I_D2 / C)
print("forgotten-factor shift minus magnetic:", sp.simplify(wrong_eta - magnetic))
print("forgotten-factor equals (1+1/g) times magnetic:", sp.simplify(wrong_eta - factor * magnetic))

# Cancellation is the same for ANY source that enters S2 the way a magnetic element does.
S_any = sp.symbols("S_any")
I_D2_any = sp.solve(sp.Eq(factor * I_D2_over_rho, S_any), I_D2_over_rho)[0]
print("generic shift - S/C:", sp.simplify(factor * I_D2_any / C - S_any / C))
