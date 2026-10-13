"""Chapter-3 hard ring: linear-ellipse time versus the invariant-curve correction.

The one-turn average on the linear Twiss ellipse is the quadratic form of the
revolution time. The nonlinear centre (x_c, a_c) = (I-R)^{-1} <M(z)-z> adds
(dT/dx) x_c + (dT/da) a_c. For a symplectic map the sum equals -2*pi*xi/C.
"""

from __future__ import annotations

import math

import numpy as np

from importlib.machinery import SourceFileLoader

v = SourceFileLoader(
    "v3", r"notes/equilibrium-energy/audit/_verify_ch3_ring.py"
).load_module()

K = 0.02
J = 1.0e-6
N_PHASE = 16


def matrix_and_orbit(elements, magnetic, delta):
    x, a, _, _, err = v.closed_orbit(delta, elements, magnetic)
    if err > 1e-8:
        raise RuntimeError(f"orbit {err}")
    mat = v.one_turn_matrix(delta, elements, magnetic, x, a)
    nu, beta, alpha = v.tune_twiss(mat)
    gamma = (1.0 + alpha * alpha) / beta
    return {
        "x": x,
        "a": a,
        "nu": nu,
        "beta": beta,
        "alpha": alpha,
        "gamma": gamma,
        "mat": mat,
        "T": None,
    }


def chromaticity(elements, magnetic):
    deltas = (-3e-4, -1.5e-4, 1.5e-4, 3e-4)
    nus = []
    for d in deltas:
        nus.append(matrix_and_orbit(elements, magnetic, d)["nu"])
    design = np.column_stack([np.ones(4), deltas, np.array(deltas) ** 2])
    return float(np.linalg.lstsq(design, np.array(nus), rcond=None)[0][1])


def averages(elements, magnetic, info, action):
    amp_x = math.sqrt(2.0 * action * info["beta"])
    amp_a = math.sqrt(2.0 * action / info["beta"])
    rel_t = []
    dxs = []
    das = []
    for phase in np.linspace(0.0, 2.0 * math.pi, N_PHASE, endpoint=False):
        xx = info["x"] + amp_x * math.cos(phase)
        aa = info["a"] - amp_a * (math.sin(phase) + info["alpha"] * math.cos(phase))
        xf, af, time, path = v.track_xt(xx, aa, 0.0, elements, magnetic)
        rel_t.append(time)
        dxs.append(xf - xx)
        das.append(af - aa)
    return float(np.mean(rel_t)), np.array([float(np.mean(dxs)), float(np.mean(das))])


def time_derivatives(elements, magnetic, info):
    eps = 1e-5
    _, _, t0, _ = v.track_xt(info["x"], info["a"], 0.0, elements, magnetic)
    _, _, tx, _ = v.track_xt(info["x"] + eps, info["a"], 0.0, elements, magnetic)
    _, _, txm, _ = v.track_xt(info["x"] - eps, info["a"], 0.0, elements, magnetic)
    _, _, ta, _ = v.track_xt(info["x"], info["a"] + eps, 0.0, elements, magnetic)
    _, _, tam, _ = v.track_xt(info["x"], info["a"] - eps, 0.0, elements, magnetic)
    dtdx = (tx - txm) / (2.0 * eps)
    dtda = (ta - tam) / (2.0 * eps)
    return t0, dtdx, dtda


def one(sextupole, magnetic):
    elements, _ = v.build_ring(K, sextupole=sextupole, magnetic_quad=magnetic)
    info = matrix_and_orbit(elements, magnetic, 0.0)
    t_ref, dtdx, dtda = time_derivatives(elements, magnetic, info)
    info["T"] = t_ref
    t_ell, kick = averages(elements, magnetic, info, J)
    # Kick was computed on the ellipse around the fixed point, so it is the
    # displacement relative to a point that already closes. Centre in that frame:
    centre = np.linalg.solve(np.eye(2) - info["mat"], kick)
    t_centre = dtdx * centre[0] + dtda * centre[1]
    xi = chromaticity(elements, magnetic)
    # Length for the prediction: reference period times reference speed is not
    # stored. Use the reference path from a on-axis turn.
    _, _, _, path = v.track_xt(info["x"], info["a"], 0.0, elements, magnetic)
    pred = -2.0 * math.pi * xi / path * t_ref  # absolute time, since Lambda*T = d<T>
    # Identity in chapter 3 is (1/T) d<T>/dJ = -2*pi*xi/C, and C = v T on the
    # reference orbit of this magnetic model, so d<T>/dJ = -2*pi*xi/C * T
    # = -2*pi*xi * T / path if path is C.
    slope_ell = (t_ell - t_ref) / J
    slope_sum = slope_ell + t_centre / J
    return {
        "xi": xi,
        "pred": pred / t_ref,  # in units of (T-Ts)/Ts per J, if path=C
        "ell": slope_ell / t_ref,
        "sum": slope_sum / t_ref,
        "C": path,
        "T": t_ref,
    }


def from_chapter3():
    """Use the chapter-3 chromaticity routine and add the centre term."""
    strengths = [0.0, 0.002, 0.004, 0.008, 0.012]
    J_loc = 1.0e-5
    for magnetic, tag in ((True, "magnetic"), (False, "cylinder")):
        print(f"--- {tag} via chapter-3 chromaticity ---")
        preds, ells, sums = [], [], []
        for s in strengths:
            elements, length = v.build_ring(K, sextupole=float(s), magnetic_quad=magnetic)
            info = v.dispersion_alpha(elements, magnetic, length)
            ch = v.chromaticity_and_action(elements, magnetic, info)
            x0, a0, t_ref, _c, _err = v.closed_orbit(0.0, elements, magnetic)
            mat = v.one_turn_matrix(0.0, elements, magnetic, x0, a0)
            _nu, beta, alpha = v.tune_twiss(mat)
            amp_x = math.sqrt(2.0 * J_loc * beta)
            amp_a = math.sqrt(2.0 * J_loc / beta)
            rels, dxs, das = [], [], []
            for phase in np.linspace(0.0, 2.0 * math.pi, 24, endpoint=False):
                xx = x0 + amp_x * math.cos(phase)
                aa = a0 - amp_a * (math.sin(phase) + alpha * math.cos(phase))
                xf, af, time, _path = v.track_xt(xx, aa, 0.0, elements, magnetic)
                rels.append((time - t_ref) / t_ref)
                dxs.append(xf - xx)
                das.append(af - aa)
            ell = float(np.mean(rels)) / J_loc
            centre = np.linalg.solve(
                np.eye(2) - mat, np.array([float(np.mean(dxs)), float(np.mean(das))])
            )
            eps = 1.0e-5

            def rel(xx, aa, elements=elements, magnetic=magnetic, t_ref=t_ref):
                _xf, _af, time, _path = v.track_xt(xx, aa, 0.0, elements, magnetic)
                return (time - t_ref) / t_ref

            dtdx = (rel(x0 + eps, a0) - rel(x0 - eps, a0)) / (2.0 * eps)
            dtda = (rel(x0, a0 + eps) - rel(x0, a0 - eps)) / (2.0 * eps)
            extra = (dtdx * centre[0] + dtda * centre[1]) / J_loc
            pred = -2.0 * math.pi * ch["xi"] / info["C_s"]
            preds.append(pred)
            ells.append(ell)
            sums.append(ell + extra)
            print(
                f"S={s:.3f} xi={ch['xi']:.4f} pred={pred:.4e} ell={ell:.4e} "
                f"sum={ell+extra:.4e} re={ell/pred:.3f} rs={(ell+extra)/pred:.3f}"
            )
        pred_a = np.array(preds)
        design = np.column_stack([pred_a, np.ones(len(pred_a))])
        for name, y in (("ellipse", np.array(ells)), ("sum", np.array(sums))):
            kappa, b = np.linalg.lstsq(design, y, rcond=None)[0]
            print(f"{tag} fit {name}: kappa={kappa:.4f} b={b:.4e}")


def main():
    from_chapter3()
    return
    strengths = [0.0, 0.002, 0.004, 0.008]
    for magnetic, tag in ((True, "magnetic"), (False, "cylinder")):
        rows = []
        print(f"--- {tag} ---")
        for s in strengths:
            row = one(s, magnetic)
            row["s"] = s
            rows.append(row)
            print(
                f"S={s:.3f} xi={row['xi']:.5f} pred={row['pred']:.6e} "
                f"ellipse={row['ell']:.6e} sum={row['sum']:.6e} "
                f"ell/pred={row['ell']/row['pred']:.4f} sum/pred={row['sum']/row['pred']:.4f}"
            )
        pred = np.array([r["pred"] for r in rows])
        for key in ("ell", "sum"):
            y = np.array([r[key] for r in rows])
            kappa, b = np.linalg.lstsq(
                np.column_stack([pred, np.ones(len(rows))]), y, rcond=None
            )[0]
            print(f"{tag} {key}: kappa={kappa:.4f} b={b:.6e}")


if __name__ == "__main__":
    main()
