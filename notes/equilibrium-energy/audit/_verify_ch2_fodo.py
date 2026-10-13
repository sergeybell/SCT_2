"""FODO ring, thin lenses: mean path lengthening vs (J/2)*int(gamma ds).

Horizontal and vertical planes are tracked together. Inside a drift the
slope is constant, and the Twiss gamma is constant, so both integrals are exact.
"""
import math


def mat_mul(a, b):
    return (
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    )


def mat_drift(length):
    return (1.0, length, 0.0, 1.0)


def mat_quad(k):
    # Delta x' = -k * x. Positive k focuses this plane.
    return (1.0, 0.0, -k, 1.0)


def twiss_from_oneturn(m):
    cos_mu = 0.5 * (m[0] + m[3])
    if abs(cos_mu) >= 1.0:
        raise RuntimeError("unstable cell, cos mu = %s" % cos_mu)
    sin_mu = math.sqrt(1.0 - cos_mu * cos_mu)
    # M12 = beta * sin(mu); take sin(mu) with the sign of M12.
    if m[1] < 0.0:
        sin_mu = -sin_mu
    beta = m[1] / sin_mu
    alpha = (m[0] - m[3]) / (2.0 * sin_mu)
    gamma = (1.0 + alpha * alpha) / beta
    mu = math.atan2(sin_mu, cos_mu)
    if mu < 0.0:
        mu += 2.0 * math.pi
    return beta, alpha, gamma, mu


def propagate_twiss(beta, alpha, gamma, m):
    # (beta, alpha, gamma) through M = [[a, b], [c, d]].
    a, b, c, d = m
    beta2 = (a * a) * beta - 2.0 * a * b * alpha + (b * b) * gamma
    alpha2 = -a * c * beta + (a * d + b * c) * alpha - b * d * gamma
    gamma2 = (c * c) * beta - 2.0 * c * d * alpha + (d * d) * gamma
    return beta2, alpha2, gamma2


def apply_vec(m, x, xp):
    return m[0] * x + m[1] * xp, m[2] * x + m[3] * xp


def main():
    n_cells = 16
    length = 0.8  # drift between thin quadrupoles, m
    k_focus = 0.40  # 1/m, horizontal focusing strength of the F quad
    j_x = 1.0e-6
    j_y = 0.7e-6
    n_turns = 8000

    # One cell: drift, F, drift, D. F focuses horizontally.
    # Elements: (matrix_x, matrix_y, drift_length or 0)
    cell = [
        (mat_drift(length), mat_drift(length), length),
        (mat_quad(k_focus), mat_quad(-k_focus), 0.0),
        (mat_drift(length), mat_drift(length), length),
        (mat_quad(-k_focus), mat_quad(k_focus), 0.0),
    ]
    ring = cell * n_cells

    # One-turn map is M_last ... M_first, so each element is multiplied on the left.
    one_x = (1.0, 0.0, 0.0, 1.0)
    one_y = (1.0, 0.0, 0.0, 1.0)
    for mx, my, _ in ring:
        one_x = mat_mul(mx, one_x)
        one_y = mat_mul(my, one_y)

    bx, ax, gx, mu_x = twiss_from_oneturn(one_x)
    by, ay, gy, mu_y = twiss_from_oneturn(one_y)
    nu_x = mu_x / (2.0 * math.pi)
    nu_y = mu_y / (2.0 * math.pi)

    # Particle with the given actions, initial phase zero.
    x = math.sqrt(2.0 * j_x * bx)
    xp = -math.sqrt(2.0 * j_x / bx) * ax
    y = math.sqrt(2.0 * j_y * by)
    yp = -math.sqrt(2.0 * j_y / by) * ay

    path_sum = 0.0
    # Also collect x'^2 at the ring entrance to test <x'^2> = J*gamma.
    xp2_sum = 0.0
    for _ in range(n_turns):
        xp2_sum += xp * xp
        for mx, my, drift in ring:
            if drift > 0.0:
                path_sum += 0.5 * (xp * xp + yp * yp) * drift
            x, xp = apply_vec(mx, x, xp)
            y, yp = apply_vec(my, y, yp)
    mean_path = path_sum / n_turns
    # recompute plane-wise gamma integrals for the expected value
    beta, alpha, gamma = bx, ax, gx
    betay, alphay, gammay = by, ay, gy
    int_gx = 0.0
    int_gy = 0.0
    for mx, my, drift in ring:
        if drift > 0.0:
            int_gx += gamma * drift
            int_gy += gammay * drift
        beta, alpha, gamma = propagate_twiss(beta, alpha, gamma, mx)
        betay, alphay, gammay = propagate_twiss(betay, alphay, gammay, my)
    expected = 0.5 * (j_x * int_gx + j_y * int_gy)
    # alpha=0 surrogate: (J/2) * int ds/beta, plane by plane
    beta, alpha, gamma = bx, ax, gx
    betay, alphay, gammay = by, ay, gy
    int_bx = 0.0
    int_by = 0.0
    n_sample = 80
    for mx, my, drift in ring:
        if drift > 0.0:
            acc = 0.0
            accy = 0.0
            for i in range(n_sample):
                s = (i + 0.5) * drift / n_sample
                md = mat_drift(s)
                b, _, _ = propagate_twiss(beta, alpha, gamma, md)
                bys, _, _ = propagate_twiss(betay, alphay, gammay, md)
                acc += 1.0 / b
                accy += 1.0 / bys
            int_bx += acc * drift / n_sample
            int_by += accy * drift / n_sample
        beta, alpha, gamma = propagate_twiss(beta, alpha, gamma, mx)
        betay, alphay, gammay = propagate_twiss(betay, alphay, gammay, my)
    expected_alpha0 = 0.5 * (j_x * int_bx + j_y * int_by)
    mean_xp2 = xp2_sum / n_turns

    print("nu_x", nu_x)
    print("nu_y", nu_y)
    print("beta_x0", bx, "alpha_x0", ax, "gamma_x0", gx)
    print("int_gamma_x", int_gx, "2*pi*nu_x", 2.0 * math.pi * nu_x)
    print("int_gamma_y", int_gy, "2*pi*nu_y", 2.0 * math.pi * nu_y)
    print("int_ds_over_beta_x", int_bx)
    print("mean_path", mean_path)
    print("expected_J_gamma", expected)
    print("ratio_to_gamma", mean_path / expected)
    print("expected_J_over_beta", expected_alpha0)
    print("ratio_to_1_over_beta", mean_path / expected_alpha0)
    print("mean_xp2_entrance", mean_xp2, "J*gamma", j_x * gx)
    print("ratio_xp2", mean_xp2 / (j_x * gx))
    # closure of Twiss
    print("twiss_closure_beta", beta, bx)


if __name__ == "__main__":
    main()
