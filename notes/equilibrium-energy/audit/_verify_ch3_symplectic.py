"""Symplectic test: mean revolution time versus Twiss action.

Phase average of one turn at a resolved amplitude, and a multi-turn average
of one particle. The prediction is -2*pi*xi/C for the time in both rings
and for the path only in the magnet. J is the linear Twiss action on the
straight, at delta = 0.
"""

from __future__ import annotations

import math

import numpy as np

import importlib.util

spec = importlib.util.spec_from_file_location(
    "v", r"notes/equilibrium-energy/audit/_verify_ch3_ring.py"
)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)

K = 0.02


def twiss_and_xi(elements, magnetic):
    x0, a0, ts, cs, err = v.closed_orbit(0.0, elements, magnetic)
    if err > 1e-8:
        raise RuntimeError(f"orbit {err}")

    def matrix(delta, x, a):
        eps = 2e-5

        def tr(xx, aa):
            return v.track_xt(xx, aa, delta, elements, magnetic)

        dx = (np.array(tr(x + eps, a)[:2]) - np.array(tr(x - eps, a)[:2])) / (2 * eps)
        da = (np.array(tr(x, a + eps)[:2]) - np.array(tr(x, a - eps)[:2])) / (2 * eps)
        return np.column_stack([dx, da])

    nu, beta, alpha = v.tune_twiss(matrix(0.0, x0, a0))
    det = float(np.linalg.det(matrix(0.0, x0, a0)))
    deltas = np.array([-3e-4, -1.5e-4, 1.5e-4, 3e-4])
    nus = []
    for d in deltas:
        x, a, _, _, e = v.closed_orbit(float(d), elements, magnetic)
        if e > 1e-7:
            raise RuntimeError(f"off-momentum orbit {e}")
        nus.append(v.tune_twiss(matrix(float(d), x, a))[0])
    design = np.column_stack([np.ones(4), deltas, deltas**2])
    coef = np.linalg.lstsq(design, np.array(nus), rcond=None)[0]
    return {
        "x0": x0,
        "a0": a0,
        "T": ts,
        "C": cs,
        "nu": nu,
        "beta": beta,
        "alpha": alpha,
        "xi": float(coef[1]),
        "det": det,
        "matrix": matrix,
    }


def ellipse(info, action, phase):
    amp_x = math.sqrt(2.0 * action * info["beta"])
    amp_a = math.sqrt(2.0 * action / info["beta"])
    xx = info["x0"] + amp_x * math.cos(phase)
    aa = info["a0"] - amp_a * (math.sin(phase) + info["alpha"] * math.cos(phase))
    return xx, aa


def phase_average(elements, magnetic, info, action, nphase=24):
    rel_t = []
    rel_c = []
    actions = []
    gamma = (1.0 + info["alpha"] ** 2) / info["beta"]
    for phase in np.linspace(0.0, 2.0 * math.pi, nphase, endpoint=False):
        xx, aa = ellipse(info, action, phase)
        xf, af, time, path = v.track_xt(xx, aa, 0.0, elements, magnetic)
        rel_t.append((time - info["T"]) / info["T"])
        rel_c.append((path - info["C"]) / info["C"])
        inv = 0.5 * (
            gamma * (xx - info["x0"]) ** 2
            + 2 * info["alpha"] * (xx - info["x0"]) * (aa - info["a0"])
            + info["beta"] * (aa - info["a0"]) ** 2
        )
        actions.append(inv)
    return float(np.mean(rel_t)), float(np.mean(rel_c)), float(np.mean(actions))


def multi_turn(elements, magnetic, info, action, nturn=80):
    xx, aa = ellipse(info, action, 0.3)
    gamma = (1.0 + info["alpha"] ** 2) / info["beta"]
    rel_t = []
    rel_c = []
    actions = []
    for _ in range(nturn):
        inv = 0.5 * (
            gamma * (xx - info["x0"]) ** 2
            + 2 * info["alpha"] * (xx - info["x0"]) * (aa - info["a0"])
            + info["beta"] * (aa - info["a0"]) ** 2
        )
        xf, af, time, path = v.track_xt(xx, aa, 0.0, elements, magnetic)
        rel_t.append((time - info["T"]) / info["T"])
        rel_c.append((path - info["C"]) / info["C"])
        actions.append(inv)
        xx, aa = xf, af
    return float(np.mean(rel_t)), float(np.mean(rel_c)), float(np.mean(actions))


def report(magnetic, sextupole, action, nturn=60):
    elements, _ = v.build_ring(K, sextupole=sextupole, magnetic_quad=magnetic)
    info = twiss_and_xi(elements, magnetic)
    pred = -2.0 * math.pi * info["xi"] / info["C"]
    mt, mc, mj = phase_average(elements, magnetic, info, action)
    nt, nc, nj = multi_turn(elements, magnetic, info, action, nturn=nturn)
    tag = "magnetic" if magnetic else "electrostatic"
    print(f"--- {tag} S={sextupole} J={action:.3e} ---")
    print(
        f"nu={info['nu']:.5f} beta={info['beta']:.3f} xi={info['xi']:.5f} "
        f"det={info['det']:.8f} C={info['C']:.4f}"
    )
    print(f"pred={pred:.6e}")
    print(
        f"phase  slopeT={mt/mj:.6e} slopeC={mc/mj:.6e} "
        f"ratioT={mt/mj/pred:.4f} ratioC={mc/mj/pred:.4f} "
        f"Jmeas={mj:.6e}"
    )
    print(
        f"turns  slopeT={nt/nj:.6e} slopeC={nc/nj:.6e} "
        f"ratioT={nt/nj/pred:.4f} ratioC={nc/nj/pred:.4f} "
        f"Jmeas={nj:.6e} N={nturn}"
    )


def main():
    report(True, 0.0, 1.0e-5, nturn=80)
    report(True, 0.0, 4.0e-5, nturn=40)
    report(True, 0.002, 1.0e-5, nturn=40)
    report(False, 0.0, 1.0e-5, nturn=80)
    report(False, 0.0, 4.0e-5, nturn=40)
    report(False, 0.0025, 1.0e-5, nturn=40)


if __name__ == "__main__":
    main()
