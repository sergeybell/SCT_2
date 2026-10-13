"""Linear optics of the model FODO: angle integral versus chromaticity.

Explains the magnetic ratio near 1.2 without another tracking run.
Quadrupoles are thin. The dipole is a sector of fixed angle; its body
matrix is the design-orbit sector matrix, and the quadrupole scales as 1/(1+delta).
"""

from __future__ import annotations

import math

import numpy as np

RHO = 20.0
THETA = math.pi / 4.0
L_HALF = 4.0
K = 0.02
N = 8


def sector(delta):
    # Body of a hard-edge sector at fixed geometry. Off momentum the
    # closed-orbit matrix used here is the standard one with rho*(1+delta):
    # length along the design arc stays rho*theta, focusing 1/rho_d**2 weakens.
    rho = RHO * (1.0 + delta)
    # Angle of the design arc is fixed, so the phase advance of the sector
    # is theta_design * rho_design / rho_particle = THETA / (1+delta).
    ang = THETA / (1.0 + delta)
    c, s = math.cos(ang), math.sin(ang)
    return np.array([[c, rho * s], [-s / rho, c]])


def drift():
    return np.array([[1.0, L_HALF], [0.0, 1.0]])


def quad(k, delta):
    return np.array([[1.0, 0.0], [-k / (1.0 + delta), 1.0]])


def oneturn(delta):
    m = np.eye(2)
    d = drift()
    for i in range(N):
        k = K if i % 2 == 0 else -K
        m = quad(k, delta) @ m
        m = d @ m
        m = sector(delta) @ m
        m = d @ m
    return m


def twiss(matrix):
    cos_mu = 0.5 * (matrix[0, 0] + matrix[1, 1])
    mu = math.acos(float(np.clip(cos_mu, -1, 1)))
    beta = matrix[0, 1] / math.sin(mu)
    alpha = (matrix[0, 0] - matrix[1, 1]) / (2.0 * math.sin(mu))
    if beta < 0:
        mu = 2 * math.pi - mu
        beta = -beta
        alpha = -alpha
    return mu / (2 * math.pi), beta, alpha


def propagate_gamma(delta=0.0):
    """Integrate gamma_twiss along the design, starting from the quad."""
    m = oneturn(delta)
    nu, beta, alpha = twiss(m)
    gamma = (1.0 + alpha**2) / beta
    # state (beta, alpha, gamma) propagation through each element, accumulate ∫gamma ds
    acc = 0.0
    b, a = beta, alpha
    for i in range(N):
        # thin quad does not add length; it jumps alpha
        k = (K if i % 2 == 0 else -K) / (1.0 + delta)
        a = a + k * b  # alpha_out = alpha_in + k*beta? 
        # Transfer of Twiss through a thin kick x' -> x'-k x:
        # beta unchanged, alpha_new = alpha_old + k*beta? Let's use the matrix.
        # Standard: M = [[1,0],[-k,1]], 
        # beta2 = m11**2 beta - 2 m11 m12 alpha + m12**2 gamma
        # For thin quad m12=0, m11=1, m21=-k, m22=1
        # beta2 = beta
        # alpha2 = -m11 m21 beta + (m11 m22+m12 m21) alpha - m12 m22 gamma
        #        = k beta + alpha
        # Yes alpha := alpha + k*beta. gamma updates.
        a = a  # already... wait I added before defining. Recalculate cleanly below.
    return nu, beta, alpha, acc


def integrate():
    nu, beta0, alpha0 = twiss(oneturn(0.0))
    b, a = beta0, alpha0
    acc = 0.0
    length = 0.0
    for i in range(N):
        k = K if i % 2 == 0 else -K
        # thin quad
        a = a + k * b
        g = (1.0 + a * a) / b
        # drift L/2: beta(s) = beta - 2 alpha s + gamma s^2
        # ∫_0^L gamma(s) ds, gamma(s) = gamma0 (constant? no)
        # gamma(s) = gamma - 2 * (something). Easier stepwise matrix.
        nstep = 40
        ds = L_HALF / nstep
        for _ in range(nstep):
            g = (1.0 + a * a) / b
            acc += g * ds
            length += ds
            # drift step
            b = b - 2 * a * ds + g * ds * ds
            a = a - g * ds
        # sector, phase ψ = sqrt(k)*s with k=1/rho**2, like a quadrupole of strength 1/rho**2
        # propagate with the sector matrix in small steps of the design arc
        nstep = 80
        ds = (RHO * THETA) / nstep
        kk = 1.0 / RHO**2
        for _ in range(nstep):
            g = (1.0 + a * a) / b
            acc += g * ds
            length += ds
            # exact drift-like step of a focusing sector over ds: matrix
            ang = ds / RHO
            c, s = math.cos(ang), math.sin(ang)
            # M = [[c, rho s], [-s/rho, c]]
            # Twiss transport
            m11, m12, m21, m22 = c, RHO * s, -s / RHO, c
            b2 = m11**2 * b - 2 * m11 * m12 * a + m12**2 * g
            a2 = -m11 * m21 * b + (m11 * m22 + m12 * m21) * a - m12 * m22 * g
            b, a = b2, a2
        # drift L/2
        for _ in range(nstep_d):
            pass
    return nu, acc, length


# The loop above is messy. Rewrite cleanly in main.
def main():
    def transport_twiss(b, a, m11, m12, m21, m22):
        g = (1.0 + a * a) / b
        b2 = m11**2 * b - 2 * m11 * m12 * a + m12**2 * g
        a2 = -m11 * m21 * b + (m11 * m22 + m12 * m21) * a - m12 * m22 * g
        return b2, a2

    nu, b, a = twiss(oneturn(0.0))
    # Check the matrix convention reproduces the one-turn Twiss by going around.
    b0, a0 = b, a
    acc = 0.0
    length = 0.0
    for i in range(N):
        k = K if (i % 2 == 0) else -K
        b, a = transport_twiss(b, a, 1.0, 0.0, -k, 1.0)
        for length_el, kind in ((L_HALF, "d"), (RHO * THETA, "s"), (L_HALF, "d")):
            nstep = 20 if kind == "d" else 40
            ds = length_el / nstep
            for _ in range(nstep):
                g = (1.0 + a * a) / b
                acc += g * ds
                length += ds
                if kind == "d":
                    b, a = transport_twiss(b, a, 1.0, ds, 0.0, 1.0)
                else:
                    ang = ds / RHO
                    c, s = math.cos(ang), math.sin(ang)
                    b, a = transport_twiss(b, a, c, RHO * s, -s / RHO, c)
    c_ring = N * (RHO * THETA + 2 * L_HALF)
    slope_angle = acc / (2.0 * c_ring)
    # chromaticity of this analytic model
    deltas = np.array([-2e-4, -1e-4, 1e-4, 2e-4])
    nus = []
    for d in deltas:
        nus.append(twiss(oneturn(float(d)))[0])
    design = np.column_stack([np.ones(4), deltas, deltas**2])
    xi = float(np.linalg.lstsq(design, np.array(nus), rcond=None)[0][1])
    pred = -2.0 * math.pi * xi / c_ring
    xi_angle = -acc / (4.0 * math.pi)
    print(f"nu={nu:.5f} beta0={b0:.3f} alpha0={a0:.3f}")
    print(f"closed Twiss error beta {b-b0:.3e} alpha {a-a0:.3e}")
    print(f"C={c_ring:.4f} int_gamma={acc:.4f}")
    print(f"slope from angles (int gamma)/(2C) = {slope_angle:.6e}")
    print(f"xi_angle from int gamma = {xi_angle:.5f}")
    print(f"xi_model={xi:.5f} pred={pred:.6e}")
    print(f"tracked magnetic slopeC={8.3616e-3:.6e} tracked xi={-0.20683:.5f}")
    print(f"angle/tracked_xi_pred = {slope_angle / ( -2*math.pi*(-0.20683)/c_ring ):.4f}")
    print(f"angle/model_pred = {slope_angle / pred:.4f}")


if __name__ == "__main__":
    main()
