"""Smooth cylindrical ring: time delay versus betatron amplitude.

Every circular orbit has the reference kinetic energy, and the linear tune
sqrt(2-beta^2) does not depend on the radius. The symplectic relation then
requires a vanishing time slope. The path slope stays finite.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp

import importlib.util

spec = importlib.util.spec_from_file_location(
    "v", r"notes/equilibrium-energy/audit/_verify_ch3_ring.py"
)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


NU = math.sqrt(2.0 - v.BETA2)
BETA_TW = v.RHO / NU
T_REF = 2.0 * math.pi * v.RHO / v.V_S
C_REF = 2.0 * math.pi * v.RHO


def track(amp, nturn=12):
    r0 = v.RHO + amp
    phi = v.potential(r0)
    gamma_k = v.GAMMA - phi
    p = math.sqrt(gamma_k * gamma_k - 1.0)
    # X, Z, PX, PZ, polar angle, path
    y0 = [r0, 0.0, 0.0, p, 0.0, 0.0]

    def rhs(_t, y):
        energy = math.sqrt(y[2] * y[2] + y[3] * y[3] + 1.0)
        r2 = y[0] * y[0] + y[1] * y[1]
        scale = -v.FIELD / r2
        vx, vz = y[2] / energy, y[3] / energy
        ang = (y[0] * vz - y[1] * vx) / r2
        speed = math.hypot(vx, vz)
        return [vx, vz, scale * y[0], scale * y[1], ang, speed]

    def hit(_t, y):
        return y[4] - 2.0 * math.pi * nturn

    hit.terminal = True
    hit.direction = 1.0
    sol = solve_ivp(
        rhs,
        (0.0, nturn * T_REF * 1.5),
        y0,
        events=hit,
        rtol=1e-11,
        atol=1e-12,
        method="DOP853",
        dense_output=False,
    )
    if sol.status != 1:
        raise RuntimeError(sol.message)
    t_end = float(sol.t_events[0][0])
    path = float(sol.y_events[0][0][5])
    # Mean radial shift from the samples.
    rad = np.hypot(sol.y[0], sol.y[1])
    return t_end / nturn, path / nturn, float(np.mean(rad) - v.RHO)


def main():
    print(f"nu={NU:.6f} beta_twiss={BETA_TW:.4f} gamma={v.GAMMA:.6f}")
    print(f"Tref={T_REF:.8f} Cref={C_REF:.8f}")
    print("amp_m  J  dT/T  dC/C  slopeT  slopeC  <x>/J  Lphi")
    for amp in (1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3):
        # Turning point of the linear ellipse: amp = sqrt(2 J beta).
        action = amp**2 / (2.0 * BETA_TW)
        t_mean, c_mean, x_mean = track(amp)
        dt = (t_mean - T_REF) / T_REF
        dc = (c_mean - C_REF) / C_REF
        print(
            f"{amp:.4f} {action:.6e} {dt:.6e} {dc:.6e} "
            f"{dt/action:.6e} {dc/action:.6e} {x_mean/action:.6e} "
            f"{(dt-dc)/action:.6e}"
        )
    # Angle part of the path: nu*J/(2*rho) is (dC/C)/J from x'^2 only.
    print(f"angle-only slopeC = {NU / (2.0 * v.RHO):.6e}")
    # xi from two nearby circular orbits is identically zero for the linear tune.
    print("linear tune is sqrt(2-beta^2) on every circular orbit, xi=0")
    print(f"symplectic slopeT = 0")


if __name__ == "__main__":
    main()
