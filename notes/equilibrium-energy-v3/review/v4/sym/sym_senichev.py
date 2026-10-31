# Senichev WEPEA036 bracket and prefactor versus eq:deq-L / eq:cmp-master.
# The Senichev expression is the one quoted in the audit task. No numeric values.
import sympy as sp

g = sp.symbols("g", positive=True)
alpha0, alpha1 = sp.symbols("alpha0 alpha1", real=True)
Lam, dm = sp.symbols("Lambda delta_m", real=True, positive=True)

beta2 = 1 - 1 / g
eta0 = alpha0 - 1 / g
eta1 = alpha1 - alpha0 / g + (3 * g - 1) / (2 * g ** 2)

# quoted parenthesis of (delta_m**2 / 2)
sen_quad = alpha1 + 3 * beta2 / (2 * g) - alpha0 / g + 1 / g ** 2
pref = g / (g * alpha0 - 1)  # gamma**2 / (gamma**2 alpha0 - 1)

delta_sen = pref * (dm ** 2 / 2 * sen_quad + Lam)
delta_note = -(1 / eta0) * (eta1 * dm ** 2 / 2 + Lam)


def main():
    print("parenthesis - eta1:", sp.simplify(sp.together(sen_quad - eta1)))
    print("prefactor - 1/eta0:", sp.simplify(sp.together(pref - 1 / eta0)))
    print("delta_sen + delta_note:", sp.simplify(sp.together(delta_sen + delta_note)))
    print("delta_sen / delta_note:", sp.simplify(sp.together(delta_sen / delta_note)))
    # kinematic piece alone: 3 beta^2/(2 gamma^2) + 1/gamma^4 vs (3 gamma^2 - 1)/(2 gamma^4)
    kin = sp.simplify(sp.together(3 * beta2 / (2 * g) + 1 / g ** 2 - (3 * g - 1) / (2 * g ** 2)))
    print("kinematic difference:", kin)


if __name__ == "__main__":
    main()
