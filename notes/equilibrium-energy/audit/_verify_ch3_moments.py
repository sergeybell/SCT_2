"""Time one revolution and sample the dispersion moments of the model ring.

eta1^es from the series in _verify_ch3_eta.py, after reduction to ring averages:

    eta1 = eta1_mag + ( <x2/rho> + 2<(D/rho)^2> - 3 alpha0 ) / gamma**2

where x = D*delta + x2*delta**2 on the design arc and the averages are
(1/C) * integral ds. Straights do not contribute to x/rho.
"""

from __future__ import annotations

import math
import time

import numpy as np

import importlib.util

spec = importlib.util.spec_from_file_location(
    "v", r"notes/equilibrium-energy/audit/_verify_ch3_ring.py"
)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)

K = 0.02
DELTA = 8.0e-4


def eta1_mag(a0, a1, g):
    b2 = 1.0 - 1.0 / g**2
    return a1 - a0 / g**2 + 1.0 / g**4 + 1.5 * b2 / g**2


def eta1_es(a0, a1, d2_avg, q, g):
    return eta1_mag(a0, a1, g) + (d2_avg + 2.0 * q - 3.0 * a0) / g**2


def sample_orbit(elements, magnetic, delta):
    x, a, _, _, err = v.closed_orbit(delta, elements, magnetic)
    if err > 1e-8:
        raise RuntimeError(f"orbit residual {err} at delta={delta}")
    pos, mom = v.state_from_xt(x, a, delta, elements[0])
    samples = []
    v.one_turn(pos, mom, elements, magnetic, samples=samples)
    return samples


def bend_profile(sol, center, alpha0, alpha1, ngrid=81):
    rel_x = sol.y[0] - center[0]
    rel_z = sol.y[1] - center[1]
    ang = np.unwrap(np.arctan2(rel_z, rel_x))
    # Design arc runs from alpha0 to alpha1, increasing.
    if ang[-1] < ang[0]:
        ang = ang + 2.0 * math.pi * (ang < ang[0])
    order = np.argsort(ang)
    ang = ang[order]
    rad = np.hypot(rel_x, rel_z)[order]
    # Drop duplicate angles from the integrator.
    keep = np.concatenate([[True], np.diff(ang) > 1e-10])
    ang = ang[keep]
    rad = rad[keep]
    grid = np.linspace(alpha0, alpha1, ngrid)
    # The stored angle can sit 2*pi away from the element angles.
    shift = round((ang[0] - alpha0) / (2.0 * math.pi)) * 2.0 * math.pi
    grid_s = grid + shift
    x = np.interp(grid_s, ang, rad - v.RHO)
    return grid, x


def moments(magnetic, sextupole, delta=DELTA, ngrid=81):
    elements, length = v.build_ring(K, sextupole=sextupole, magnetic_quad=magnetic)
    info = v.dispersion_alpha(elements, magnetic, length)
    t0 = time.perf_counter()
    plus = sample_orbit(elements, magnetic, delta)
    minus = sample_orbit(elements, magnetic, -delta)
    # Reference orbit, used only as a sanity check that x vanishes.
    zero = sample_orbit(elements, magnetic, 0.0)
    bends_p = [s for s in plus if True]
    acc_d = 0.0
    acc_q = 0.0
    acc_d2 = 0.0
    acc_angle = 0.0
    x_zero_max = 0.0
    for (el_p, sol_p), (el_m, sol_m), (el_0, sol_0) in zip(plus, minus, zero):
        grid, xp = bend_profile(
            sol_p, el_p["center"], el_p["alpha0"], el_p["alpha1"], ngrid=ngrid
        )
        _, xm = bend_profile(
            sol_m, el_m["center"], el_m["alpha0"], el_m["alpha1"], ngrid=ngrid
        )
        _, xz = bend_profile(
            sol_0, el_0["center"], el_0["alpha0"], el_0["alpha1"], ngrid=ngrid
        )
        d_over_rho = (xp - xm) / (2.0 * delta) / v.RHO
        x2_over_rho = (xp + xm) / (2.0 * delta**2) / v.RHO
        # ∫ f ds = ∫ f * rho dθ. trapz, not a rectangular sum of the nodes.
        acc_d += float(np.trapezoid(d_over_rho, grid) * v.RHO)
        acc_q += float(np.trapezoid(d_over_rho**2, grid) * v.RHO)
        acc_d2 += float(np.trapezoid(x2_over_rho, grid) * v.RHO)
        acc_angle += float(grid[-1] - grid[0])
        x_zero_max = max(x_zero_max, float(np.max(np.abs(xz))))
    c_ref = info["C_s"]
    a0 = acc_d / c_ref
    q = acc_q / c_ref
    d2 = acc_d2 / c_ref
    g = v.GAMMA
    # Path check: ∫ x dθ on the +delta orbit versus the tracked length.
    bend_extra = 0.0
    for el_p, sol_p in plus:
        grid_b, x_b = bend_profile(
            sol_p, el_p["center"], el_p["alpha0"], el_p["alpha1"], ngrid=ngrid
        )
        bend_extra += float(np.trapezoid(x_b, grid_b))
    _, _, _, path_plus, _ = v.closed_orbit(delta, elements, magnetic)
    path_lin = (path_plus - c_ref) / delta
    pred = eta1_es(info["alpha0"], info["alpha1"], d2, q, g)
    tag = "magnetic" if magnetic else "electrostatic"
    print(f"--- {tag} S={sextupole}  ({time.perf_counter()-t0:.1f}s) ---")
    print(
        f"C={c_ref:.4f} angle_sum={acc_angle:.6f} x0_max={x_zero_max:.3e}"
    )
    print(
        f"alpha0_track={info['alpha0']:.8f} alpha0_moments={a0:.8f} "
        f"alpha1={info['alpha1']:.8f}"
    )
    print(
        f"path_excess/delta={path_lin:.8f} bend_int_x_dth/delta={bend_extra/delta:.8f}"
    )
    print(f"<(D/rho)^2>={q:.8f} <x2/rho>={d2:.8f}")
    print(
        f"eta1_track={info['eta1']:.8f} eta1_mag={eta1_mag(info['alpha0'], info['alpha1'], g):.8f}"
    )
    print(f"eta1_es_formula={pred:.8f} residual={pred - info['eta1']:.3e}")
    if not magnetic:
        # Same moments must NOT be fed to the magnetic ring; printed for the cylinder only.
        pass
    return info


def main():
    # One turn timing.
    elements, _ = v.build_ring(K, sextupole=0.0, magnetic_quad=True)
    t0 = time.perf_counter()
    v.closed_orbit(0.0, elements, True)
    print(f"one closed orbit (magnetic, includes Newton): {time.perf_counter()-t0:.2f}s")
    for delta in (2.0e-4, 5.0e-4, 8.0e-4):
        moments(False, 0.0, delta=delta, ngrid=161)


if __name__ == "__main__":
    main()
