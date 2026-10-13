"""Numbers for chapter 3 after every quadrupole is kicked once per turn."""

import importlib.util
import math

import numpy as np

spec = importlib.util.spec_from_file_location(
    "v", r"notes/equilibrium-energy/audit/_verify_ch3_ring.py"
)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)

K = 0.02


def summary(magnetic, S=0.0, J=1e-8):
    el, leng = v.build_ring(K, sextupole=S, magnetic_quad=magnetic)
    info = v.dispersion_alpha(el, magnetic, leng)
    x0, a0, Ts, Cs, _ = v.closed_orbit(0.0, el, magnetic)
    eps = 2e-5

    def tr(x, a, d):
        return v.track_xt(x, a, d, el, magnetic)

    def mat(d, x, a):
        dx = (np.array(tr(x + eps, a, d)[:2]) - np.array(tr(x - eps, a, d)[:2])) / (2 * eps)
        da = (np.array(tr(x, a + eps, d)[:2]) - np.array(tr(x, a - eps, d)[:2])) / (2 * eps)
        return np.column_stack([dx, da])

    nu, beta, alpha = v.tune_twiss(mat(0.0, x0, a0))
    ds = np.array([-3e-4, -1.5e-4, 1.5e-4, 3e-4])
    nus = []
    for d in ds:
        x, a, _, _, _ = v.closed_orbit(float(d), el, magnetic)
        nus.append(v.tune_twiss(mat(float(d), x, a))[0])
    design = np.column_stack([np.ones(4), ds, ds**2])
    coef = np.linalg.lstsq(design, np.array(nus), rcond=None)[0]
    xi = float(coef[1])
    relt, relc = [], []
    for ph in np.linspace(0.0, 2.0 * math.pi, 12, endpoint=False):
        xx = x0 + math.sqrt(2 * J * beta) * math.cos(ph)
        aa = a0 - math.sqrt(2 * J / beta) * (math.sin(ph) + alpha * math.cos(ph))
        _, _, T, path = tr(xx, aa, 0.0)
        relt.append((T - Ts) / Ts)
        relc.append((path - Cs) / Cs)
    st = float(np.mean(relt) / J)
    sc = float(np.mean(relc) / J)
    pred = -2.0 * math.pi * xi / Cs
    g2 = v.GAMMA**2
    g4 = g2 * g2
    b2 = v.BETA2
    a0, a1 = info["alpha0"], info["alpha1"]
    tag = "magnetic" if magnetic else "electrostatic"
    print(f"--- {tag} S={S} ---")
    print(f"nu={nu:.5f} beta={beta:.3f} alpha_tw={alpha:.3f} xi={xi:.5f}")
    print(f"alpha0={a0:.8f} alpha1={a1:.8f}")
    print(f"eta0={info['eta0']:.8f} eta1={info['eta1']:.8f} resid={info['eta_resid']:.2e}")
    print(f"eta0_mag={a0 - 1/g2:.8f} eta1_mag={a1 - a0/g2 + 1/g4 + 1.5*b2/g2:.8f}")
    print(f"eta0_es={(2 - b2)*a0 - 1/g2:.8f}")
    ratio_t = st / pred if abs(pred) > 1e-8 else float("nan")
    ratio_c = sc / pred if abs(pred) > 1e-8 else float("nan")
    print(
        f"slopeT={st:.6e} slopeC={sc:.6e} pred={pred:.6e} "
        f"ratioT={ratio_t:.4f} ratioC={ratio_c:.4f}"
    )


if __name__ == "__main__":
    summary(True)
    summary(False)
    print("electrostatic sextupoles")
    for S in (0.002, 0.004, 0.008, 0.012):
        try:
            summary(False, S)
        except Exception as exc:
            print("S", S, type(exc).__name__, exc)
    print("magnetic sextupoles")
    for S in (0.001, 0.002, 0.004):
        try:
            summary(True, S)
        except Exception as exc:
            print("S", S, type(exc).__name__, exc)
