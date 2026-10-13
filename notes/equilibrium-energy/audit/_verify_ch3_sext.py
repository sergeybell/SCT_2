"""Sextupole scan: does zero chromaticity zero the revolution-time shift?"""

import importlib.util
import math

import numpy as np

spec = importlib.util.spec_from_file_location(
    "v", r"notes/equilibrium-energy/audit/_verify_ch3_ring.py"
)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def measure(S, J=2e-8):
    el, _ = v.build_ring(0.12, sextupole=S, magnetic_quad=False)
    x0, a0, Ts, Cs, err = v.closed_orbit(0.0, el, False)
    if err > 1e-8:
        raise RuntimeError(f"orbit {err}")

    def mat(delta, x, a):
        eps = 2e-5

        def tr(xx, aa):
            return v.track_xt(xx, aa, delta, el, False)

        dx = (np.array(tr(x + eps, a)[:2]) - np.array(tr(x - eps, a)[:2])) / (2 * eps)
        da = (np.array(tr(x, a + eps)[:2]) - np.array(tr(x, a - eps)[:2])) / (2 * eps)
        return np.column_stack([dx, da])

    nu, beta, alpha = v.tune_twiss(mat(0.0, x0, a0))
    nus = []
    for d in (-2e-4, 2e-4):
        x, a, _, _, _ = v.closed_orbit(d, el, False)
        nus.append(v.tune_twiss(mat(d, x, a))[0])
    xi = (nus[1] - nus[0]) / 4e-4
    relt, relc = [], []
    for ph in np.linspace(0.0, 2.0 * math.pi, 12, endpoint=False):
        xx = x0 + math.sqrt(2 * J * beta) * math.cos(ph)
        aa = a0 - math.sqrt(2 * J / beta) * (math.sin(ph) + alpha * math.cos(ph))
        _, _, T, path = v.track_xt(xx, aa, 0.0, el, False)
        relt.append((T - Ts) / Ts)
        relc.append((path - Cs) / Cs)
    st = float(np.mean(relt) / J)
    sc = float(np.mean(relc) / J)
    pred = -2.0 * math.pi * xi / Cs
    ratio = st / pred if abs(pred) > 1e-5 else float("nan")
    print(
        f"S={S:+.4f} nu={nu:.4f} xi={xi:+.4f} slopeT={st:.4e} "
        f"slopeC={sc:.4e} pred={pred:.4e} ratio={ratio:.3f} Lphi={st - sc:.4e}"
    )


for S in np.linspace(0.02, 0.14, 7):
    try:
        measure(float(S))
    except Exception as exc:
        print(f"S={S:+.4f} FAIL {type(exc).__name__}: {exc}")
