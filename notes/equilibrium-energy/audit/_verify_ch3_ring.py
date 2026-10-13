"""Chapter 3: revolution period in a cylindrical electrostatic ring.

Exact relativistic tracking in the plane, hard-edge cylindrical deflectors
(radial field E proportional to 1/r plus the longitudinal fringe kick required
by the potential step), thin quadrupoles and an optional thin sextupole.
A magnetic ring with the same geometry is tracked the same way.

Prints eta0, eta1, the betatron time coefficient, chromaticity, and the
symplectic relation dT/dJ = -2*pi*xi/v.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

# Proton, 233 MeV. Rest energy is the unit, c = 1, length in metres.
M_REST_MEV = 938.2720813
EKIN_MEV = 233.0
GAMMA = 1.0 + EKIN_MEV / M_REST_MEV
BETA = math.sqrt(1.0 - 1.0 / GAMMA**2)
BETA2 = BETA * BETA
P_S = BETA * GAMMA
V_S = BETA
# q = 1. Electrode scale: q*Phi = BETA2*GAMMA*ln(r/rho), E_r = -BETA2*GAMMA/r.
FIELD = BETA2 * GAMMA

RHO = 20.0
N_CELL = 8
THETA = 2.0 * math.pi / N_CELL
# Straight length per cell, split equally before and after the bend.
L_STRAIGHT = 8.0
# Horizontal focusing strength of the thin quad, 1/m. Alternates in sign.
K_QUAD = 0.12

RTOL = 1e-11
ATOL = 1e-12


def pack_report(lines):
    text = "\n".join(lines)
    print(text)
    return text


def gamma_of_p(p_mag):
    return math.sqrt(p_mag * p_mag + 1.0)


def p_of_delta(delta):
    return P_S * (1.0 + delta)


def total_energy_of_delta(delta):
    """Relativistic energy in the field-free gap, in units of the rest energy."""
    return gamma_of_p(p_of_delta(delta))


def build_ring(k_quad, sextupole=0.0, magnetic_quad=False):
    """Element list for one turn, starting on the first quadrupole plane.

    sextupole is the geometric thin-lens strength S in the kick
    Delta x' = (S/2) x^2 / (1+delta), applied on the focusing quad planes
    (even cells) as a magnetic sextupole.
    """
    pos = np.array([0.0, 0.0])
    phi = 0.0
    elements = []
    for i in range(N_CELL):
        normal = np.array([math.cos(phi), math.sin(phi)])
        trans = np.array([-math.sin(phi), math.cos(phi)])
        sign = 1.0 if (i % 2 == 0) else -1.0
        s_here = sextupole if (i % 2 == 0) else 0.0
        elements.append(
            {
                "kind": "quad",
                "Q": pos.copy(),
                "n": normal.copy(),
                "t": trans.copy(),
                "k": sign * k_quad,
                "S": s_here,
                "magnetic_quad": magnetic_quad,
            }
        )
        # Straight L/2 from the quad to the entrance of the arc.
        pos = pos + (L_STRAIGHT / 2.0) * normal
        center = pos + RHO * np.array(
            [math.cos(phi + math.pi / 2.0), math.sin(phi + math.pi / 2.0)]
        )
        alpha0 = phi - math.pi / 2.0
        elements.append(
            {
                "kind": "bend",
                "center": center.copy(),
                "alpha0": alpha0,
                "alpha1": alpha0 + THETA,
                "Q_in": pos.copy(),
                "n_in": normal.copy(),
            }
        )
        alpha1 = alpha0 + THETA
        pos = center + RHO * np.array([math.cos(alpha1), math.sin(alpha1)])
        phi = alpha1 + math.pi / 2.0
        normal_out = np.array([math.cos(phi), math.sin(phi)])
        pos = pos + (L_STRAIGHT / 2.0) * normal_out
    if abs(phi - 2.0 * math.pi) > 1e-8 or np.linalg.norm(pos) > 1e-8:
        raise RuntimeError(f"reference did not close: pos={pos}, phi={phi}")
    length = N_CELL * (RHO * THETA + L_STRAIGHT)
    return elements, length


def potential(r):
    return FIELD * math.log(r / RHO)


def fringe_kick(pos, mom, center, entering):
    """Longitudinal momentum jump from the potential step. p_radial stays."""
    rel = pos - center
    r = np.hypot(rel[0], rel[1])
    r_hat = rel / r
    theta_hat = np.array([-r_hat[1], r_hat[0]])
    p_r = float(np.dot(mom, r_hat))
    p_th = float(np.dot(mom, theta_hat))
    e_kin = math.sqrt(p_r * p_r + p_th * p_th + 1.0)
    phi = potential(r)
    e_new = e_kin - phi if entering else e_kin + phi
    p2 = e_new * e_new - 1.0 - p_r * p_r
    if p2 <= 0.0:
        raise RuntimeError("fringe kick turned the particle around")
    p_th_new = math.copysign(math.sqrt(p2), p_th if p_th != 0.0 else 1.0)
    return p_r * r_hat + p_th_new * theta_hat


def drift_to_plane(pos, mom, Q, normal):
    energy = math.sqrt(float(np.dot(mom, mom)) + 1.0)
    vel = mom / energy
    denom = float(np.dot(vel, normal))
    if denom <= 1e-14:
        raise RuntimeError("drift parallel to the target plane")
    time = float(np.dot(Q - pos, normal) / denom)
    if time < -1e-9:
        raise RuntimeError(f"plane is behind the particle, t={time}")
    time = max(time, 0.0)
    new_pos = pos + vel * time
    path = float(np.linalg.norm(vel) * time)
    return new_pos, time, path


def bend_field(_t, y, center):
    rel_x = y[0] - center[0]
    rel_z = y[1] - center[1]
    r2 = rel_x * rel_x + rel_z * rel_z
    energy = math.sqrt(y[2] * y[2] + y[3] * y[3] + 1.0)
    # E = -(FIELD / r^2) * r_vec
    scale = -FIELD / r2
    return [y[2] / energy, y[3] / energy, scale * rel_x, scale * rel_z]


def bend_magnetic(_t, y, b_field):
    energy = math.sqrt(y[2] * y[2] + y[3] * y[3] + 1.0)
    vx = y[2] / energy
    vz = y[3] / energy
    # B = (0, B, 0), F = q v x B = (-vz B, vx B)
    return [vx, vz, -vz * b_field, vx * b_field]


def track_bend(pos, mom, center, alpha0, alpha1, magnetic):
    # y = X, Z, PX, PZ, path
    y0 = [pos[0], pos[1], mom[0], mom[1], 0.0]

    def cross(_t, y):
        rel_x = y[0] - center[0]
        rel_z = y[1] - center[1]
        return rel_x * math.sin(alpha1) - rel_z * math.cos(alpha1)

    cross.terminal = True
    cross.direction = -1.0

    if magnetic:
        b_field = P_S / RHO

        def rhs(t, y):
            base = bend_magnetic(t, y, b_field)
            speed = math.hypot(base[0], base[1])
            return base + [speed]

    else:

        def rhs(t, y):
            base = bend_field(t, y, center)
            speed = math.hypot(base[0], base[1])
            return base + [speed]

    t_max = 8.0 * RHO * abs(alpha1 - alpha0) / V_S + 1.0
    sol = solve_ivp(
        rhs,
        (0.0, t_max),
        y0,
        events=cross,
        rtol=RTOL,
        atol=ATOL,
        dense_output=False,
        method="DOP853",
    )
    if sol.status != 1 or len(sol.t_events[0]) == 0:
        raise RuntimeError(f"bend event failed, status={sol.status}, message={sol.message}")
    yf = sol.y[:, -1]
    time = float(sol.t_events[0][0])
    return np.array([yf[0], yf[1]]), np.array([yf[2], yf[3]]), time, float(yf[4]), sol


def quad_kick(pos, mom, element, delta_for_scale):
    """Thin lens that preserves the momentum modulus.

    A magnetic quadrupole scales as 1/p. An electrostatic quadrupole scales
    as 1/(p v), the electric rigidity, and the particle leaves the lens with
    the same modulus it entered with because the quadrupole potential vanishes
    on axis and the lens is crossed completely. The sextupole is magnetic.
    """
    del delta_for_scale
    normal = element["n"]
    trans = element["t"]
    x = float(np.dot(pos - element["Q"], trans))
    p_n = float(np.dot(mom, normal))
    p_t = float(np.dot(mom, trans))
    p_mag = math.hypot(p_n, p_t)
    xprime = p_t / p_n
    if element["magnetic_quad"]:
        dk = element["k"] * (P_S / p_mag)
    else:
        energy = math.sqrt(p_mag * p_mag + 1.0)
        pv = p_mag * p_mag / energy
        dk = element["k"] * ((P_S * V_S) / pv)
    xprime_new = xprime - dk * x
    if element["S"] != 0.0:
        xprime_new += 0.5 * element["S"] * x * x * (P_S / p_mag)
    p_n_new = p_mag / math.sqrt(1.0 + xprime_new * xprime_new)
    p_t_new = xprime_new * p_n_new
    if p_n < 0.0:
        p_n_new = -p_n_new
        p_t_new = -p_t_new
    return p_n_new * normal + p_t_new * trans


def one_turn(pos, mom, elements, magnetic_bend, samples=None):
    time = 0.0
    path = 0.0
    # The particle starts on the first quadrupole plane, before the kick.
    # Every quadrupole is kicked once; the closing drift lands on that plane
    # again and does not kick it a second time.
    n = len(elements)
    for i, el in enumerate(elements):
        nxt = elements[(i + 1) % n]
        if el["kind"] == "quad":
            mom = quad_kick(pos, mom, el, 0.0)
            pos, dt, dp = drift_to_plane(pos, mom, nxt["Q_in"], nxt["n_in"])
            time += dt
            path += dp
        elif el["kind"] == "bend":
            if not magnetic_bend:
                mom = fringe_kick(pos, mom, el["center"], entering=True)
            pos, mom, dt, dp, sol = track_bend(
                pos, mom, el["center"], el["alpha0"], el["alpha1"], magnetic_bend
            )
            time += dt
            path += dp
            if samples is not None:
                samples.append((el, sol))
            if not magnetic_bend:
                mom = fringe_kick(pos, mom, el["center"], entering=False)
            # Fly to the next quad plane.
            pos, dt, dp = drift_to_plane(pos, mom, nxt["Q"], nxt["n"])
            time += dt
            path += dp
        else:
            raise RuntimeError(el["kind"])
    return pos, mom, time, path


def state_from_xt(x, angle, delta, element0):
    normal = element0["n"]
    trans = element0["t"]
    pos = element0["Q"] + x * trans
    direction = normal * math.cos(angle) + trans * math.sin(angle)
    mom = p_of_delta(delta) * direction
    return pos, mom


def xt_from_state(pos, mom, element0):
    x = float(np.dot(pos - element0["Q"], element0["t"]))
    p_n = float(np.dot(mom, element0["n"]))
    p_t = float(np.dot(mom, element0["t"]))
    angle = math.atan2(p_t, p_n)
    p_mag = math.sqrt(float(np.dot(mom, mom)))
    return x, angle, p_mag


def track_xt(x, angle, delta, elements, magnetic_bend):
    pos, mom = state_from_xt(x, angle, delta, elements[0])
    pos, mom, time, path = one_turn(pos, mom, elements, magnetic_bend)
    xf, af, _ = xt_from_state(pos, mom, elements[0])
    return xf, af, time, path


def closed_orbit(delta, elements, magnetic_bend):
    def residual(vec):
        xf, af, _, _ = track_xt(float(vec[0]), float(vec[1]), delta, elements, magnetic_bend)
        return np.array([xf - vec[0], af - vec[1]])

    # One Newton step from the origin gives the linear closed orbit. Further
    # steps pick up the quadratic dispersion. Steps that leave the aperture
    # are rejected.
    vec = np.zeros(2)
    err = 1.0
    for _ in range(8):
        r0 = residual(vec)
        err = float(np.linalg.norm(r0))
        if err < 1e-12:
            break
        eps = 1e-6
        jac = np.column_stack(
            [
                (residual(vec + np.array([eps, 0.0])) - r0) / eps,
                (residual(vec + np.array([0.0, eps])) - r0) / eps,
            ]
        )
        step = np.linalg.solve(jac, -r0)
        scale = 1.0
        accepted = False
        for _trial in range(8):
            trial = vec + scale * step
            try:
                r_try = residual(trial)
            except Exception:
                scale *= 0.5
                continue
            if np.linalg.norm(r_try) < err:
                vec = trial
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            break
    x, a = float(vec[0]), float(vec[1])
    xf, af, time, path = track_xt(x, a, delta, elements, magnetic_bend)
    return x, a, time, path, math.hypot(xf - x, af - a)


def one_turn_matrix(delta, elements, magnetic_bend, x0, a0):
    eps = 1e-7
    base = track_xt(x0, a0, delta, elements, magnetic_bend)
    dx = track_xt(x0 + eps, a0, delta, elements, magnetic_bend)
    da = track_xt(x0, a0 + eps, delta, elements, magnetic_bend)
    m00 = (dx[0] - base[0]) / eps
    m10 = (dx[1] - base[1]) / eps
    m01 = (da[0] - base[0]) / eps
    m11 = (da[1] - base[1]) / eps
    return np.array([[m00, m01], [m10, m11]])


def tune_twiss(matrix):
    cos_mu = 0.5 * (matrix[0, 0] + matrix[1, 1])
    cos_mu = float(np.clip(cos_mu, -1.0, 1.0))
    mu = math.acos(cos_mu)
    sin_mu = math.sin(mu)
    if abs(sin_mu) < 1e-12:
        raise RuntimeError("tune too close to an integer")
    beta = matrix[0, 1] / sin_mu
    alpha = (matrix[0, 0] - matrix[1, 1]) / (2.0 * sin_mu)
    # Fractional tune in [0, 1). The sign of sin follows beta>0.
    if beta < 0:
        mu = 2.0 * math.pi - mu
        beta = -beta
        alpha = -alpha
    nu = mu / (2.0 * math.pi)
    return nu, beta, alpha


def polyfit_eta(deltas, periods, t_ref):
    rel = np.array([(t - t_ref) / t_ref for t in periods])
    # Fit through the origin in the odd part plus even powers: a1 d + a2 d^2 + a3 d^3.
    design = np.column_stack([deltas, deltas**2, deltas**3])
    coef, *_ = np.linalg.lstsq(design, rel, rcond=None)
    pred = design @ coef
    err = float(np.max(np.abs(pred - rel)))
    return coef, err


def smooth_circle_times():
    """Continuous cylindrical field, no straights. Analytic benchmark."""
    lines = ["== smooth cylindrical ring, no straights =="]
    rho = RHO

    def circular_time(delta):
        e_gap = total_energy_of_delta(delta)
        eps = e_gap / GAMMA - 1.0
        r = rho * math.exp(eps / BETA2)
        # Kinetic energy on every circular orbit equals the reference.
        return 2.0 * math.pi * r / V_S

    t_ref = 2.0 * math.pi * rho / V_S
    deltas = np.array([-1e-3, -5e-4, -2e-4, 2e-4, 5e-4, 1e-3])
    times = np.array([circular_time(d) for d in deltas])
    # Also integrate the ODE for one value as a check of the field normalisation.
    delta = 1.0e-3
    e_gap = total_energy_of_delta(delta)
    eps = e_gap / GAMMA - 1.0
    r = rho * math.exp(eps / BETA2)
    # Placed at (r, 0), tangential momentum +z, kinetic energy = GAMMA.
    p_tang = math.sqrt(GAMMA * GAMMA - 1.0)
    y0 = [r, 0.0, 0.0, p_tang, 0.0]  # X, Z, PX, PZ, alpha

    def rhs(_t, y):
        energy = math.sqrt(y[2] * y[2] + y[3] * y[3] + 1.0)
        r2 = y[0] * y[0] + y[1] * y[1]
        scale = -FIELD / r2
        vx, vz = y[2] / energy, y[3] / energy
        alpha_dot = (y[0] * vz - y[1] * vx) / r2
        return [vx, vz, scale * y[0], scale * y[1], alpha_dot]

    def hit(_t, y):
        return y[4] - 2.0 * math.pi

    hit.terminal = True
    hit.direction = 1.0
    sol = solve_ivp(rhs, (0.0, 5.0 * t_ref), y0, events=hit, rtol=RTOL, atol=ATOL, method="DOP853")
    t_num = float(sol.t_events[0][0])
    h_end = math.sqrt(sol.y[2, -1] ** 2 + sol.y[3, -1] ** 2 + 1.0) + potential(
        math.hypot(sol.y[0, -1], sol.y[1, -1])
    )
    h_start = GAMMA + potential(r)
    coef, err = polyfit_eta(deltas, times, t_ref)
    # (T-Ts)/Ts = exp(eps/beta^2)-1, eps = sqrt(1+beta^2(2d+d^2))-1
    # series: delta + (1/2 + 1/(2*gamma^2)) delta^2
    eta0_exact = 1.0
    eta1_exact = 0.5 + 0.5 / GAMMA**2
    lines.append(f"gamma={GAMMA:.6f} beta2={BETA2:.6f} 1/gamma2={1/GAMMA**2:.6f}")
    lines.append(f"ODE one-turn time vs analytic: {t_num:.12e} / {circular_time(delta):.12e}")
    lines.append(f"energy drift over the turn: {h_end - h_start:.3e}")
    lines.append(f"fit eta0={coef[0]:.8f} exact={eta0_exact:.8f}")
    lines.append(f"fit eta1={coef[1]:.8f} exact={eta1_exact:.8f} resid={err:.3e}")
    lines.append(
        f"magnetic smooth with alpha0=1 would give eta0={BETA2:.6f}, not {eta0_exact:.6f}"
    )
    return pack_report(lines)


def dispersion_alpha(elements, magnetic_bend, circumference):
    """Geometric alpha0, alpha1 from the path length of closed orbits."""
    deltas = np.array([-1e-3, -5e-4, -2e-4, 2e-4, 5e-4, 1e-3])
    paths = []
    times = []
    for d in deltas:
        _, _, time, path, err = closed_orbit(d, elements, magnetic_bend)
        if err > 1e-8:
            raise RuntimeError(f"orbit residual {err}")
        paths.append(path)
        times.append(time)
    _, _, t_ref, c_ref, err0 = closed_orbit(0.0, elements, magnetic_bend)
    path_coef, path_err = polyfit_eta(deltas, paths, c_ref)
    time_coef, time_err = polyfit_eta(deltas, times, t_ref)
    return {
        "alpha0": path_coef[0],
        "alpha1": path_coef[1],
        "alpha_resid": path_err,
        "eta0": time_coef[0],
        "eta1": time_coef[1],
        "eta_resid": time_err,
        "T_s": t_ref,
        "C_s": c_ref,
        "C_geom": circumference,
        "orbit_err0": err0,
    }


def chromaticity_and_action(elements, magnetic_bend, info):
    deltas = np.array([-8e-4, -4e-4, -2e-4, 2e-4, 4e-4, 8e-4])
    nus = []
    for d in deltas:
        x, a, _, _, _ = closed_orbit(d, elements, magnetic_bend)
        matrix = one_turn_matrix(d, elements, magnetic_bend, x, a)
        nu, _, _ = tune_twiss(matrix)
        nus.append(nu)
    nus = np.array(nus)
    # nu = nu0 + xi * delta + xi2 * delta^2
    design = np.column_stack([np.ones_like(deltas), deltas, deltas**2])
    coef, *_ = np.linalg.lstsq(design, nus, rcond=None)
    nu0_fit = coef[0]
    xi = coef[1]

    x, a, _, _, _ = closed_orbit(0.0, elements, magnetic_bend)
    matrix = one_turn_matrix(0.0, elements, magnetic_bend, x, a)
    nu, beta, alpha = tune_twiss(matrix)

    actions = np.array([1e-6, 2e-6, 4e-6])
    rows = []
    for j in actions:
        # Average over initial phases. 6 phases, 1 turn each, on the linear ellipse.
        rel_t = []
        rel_c = []
        for phase in np.linspace(0.0, 2.0 * math.pi, 8, endpoint=False):
            amp_x = math.sqrt(2.0 * j * beta)
            amp_a = math.sqrt(2.0 * j / beta)
            xx = x + amp_x * math.cos(phase)
            aa = a - amp_a * (math.sin(phase) + alpha * math.cos(phase))
            _, _, time, path = track_xt(xx, aa, 0.0, elements, magnetic_bend)
            rel_t.append((time - info["T_s"]) / info["T_s"])
            rel_c.append((path - info["C_s"]) / info["C_s"])
        rows.append((j, float(np.mean(rel_t)), float(np.mean(rel_c))))
    # Slope versus J.
    js = np.array([r[0] for r in rows])
    lt = np.array([r[1] for r in rows])
    lc = np.array([r[2] for r in rows])
    slope_t = float(np.linalg.lstsq(js.reshape(-1, 1), lt, rcond=None)[0][0])
    slope_c = float(np.linalg.lstsq(js.reshape(-1, 1), lc, rcond=None)[0][0])
    predicted = -2.0 * math.pi * xi / info["C_s"]
    return {
        "nu": nu,
        "nu_fit": nu0_fit,
        "beta": beta,
        "alpha": alpha,
        "xi": xi,
        "det": float(np.linalg.det(matrix)),
        "slope_time": slope_t,
        "slope_path": slope_c,
        "predicted_time": predicted,
        "phases": rows,
        "matrix_trace": float(np.trace(matrix)),
    }


def scan_sextupole(base_k, magnetic_bend, magnetic_quad):
    """Find a thin-sextupole strength that zeroes horizontal chromaticity."""
    # Coarse scan, then a linear interpolation.
    strengths = np.linspace(-0.4, 0.4, 9)
    records = []
    for s in strengths:
        elements, length = build_ring(base_k, sextupole=float(s), magnetic_quad=magnetic_quad)
        try:
            info = dispersion_alpha(elements, magnetic_bend, length)
            chrom = chromaticity_and_action(elements, magnetic_bend, info)
        except Exception as exc:  # noqa: BLE001 — scan may cross a resonance
            records.append((float(s), None, str(exc)))
            continue
        records.append((float(s), chrom["xi"], chrom))
    return records


def formula_block(info, chrom, label):
    a0 = info["alpha0"]
    a1 = info["alpha1"]
    g2 = GAMMA**2
    g4 = g2 * g2
    eta0_mag = a0 - 1.0 / g2
    eta0_es = (2.0 - BETA2) * a0 - 1.0 / g2
    # equivalent: a0 - (1-a0)/g2
    eta1_mag = a1 - a0 / g2 + 1.0 / g4 + 1.5 * BETA2 / g2
    # Senichev (25): F_v chosen so that his eta0 matches the tracked eta0.
    # alpha0 - 1/(F g^2) = eta0_tracked => F = 1/(g^2 (alpha0 - eta0))
    denom = a0 - info["eta0"]
    if abs(denom) < 1e-14:
        f_v = float("inf")
    else:
        f_v = 1.0 / (g2 * denom)
    eta1_sen = a1 - a0 / (f_v * g2) + 1.0 / (f_v**2 * g4)
    lines = [f"== {label} =="]
    lines.append(
        f"C_track={info['C_s']:.6f} C_geom={info['C_geom']:.6f} T={info['T_s']:.6f} "
        f"orbit_err={info['orbit_err0']:.2e}"
    )
    lines.append(
        f"alpha0={a0:.8f} alpha1={a1:.8f} path_resid={info['alpha_resid']:.2e}"
    )
    lines.append(
        f"eta0_track={info['eta0']:.8f} eta1_track={info['eta1']:.8f} "
        f"time_resid={info['eta_resid']:.2e}"
    )
    lines.append(f"eta0_magnetic_formula={eta0_mag:.8f}  (alpha0 - 1/gamma^2)")
    lines.append(f"eta0_electrostatic={(eta0_es):.8f}  ((2-beta^2)*alpha0 - 1/gamma^2)")
    lines.append(f"eta0_es alt={a0 - (1.0 - a0) / g2:.8f}  (alpha0 - (1-alpha0)/gamma^2)")
    lines.append(f"eta1_magnetic_formula={eta1_mag:.8f}")
    lines.append(f"F_v fitted to eta0={f_v:.6f}")
    lines.append(f"eta1_Senichev25_with_that_Fv={eta1_sen:.8f}")
    lines.append(
        f"nu={chrom['nu']:.6f} beta_twiss={chrom['beta']:.4f} alpha_twiss={chrom['alpha']:.4f} "
        f"xi={chrom['xi']:.6f} detM={chrom['det']:.8f}"
    )
    lines.append(
        f"v_s*T={V_S * info['T_s']:.6f} versus C_track (equal when |p| is the reference)"
    )
    lines.append(
        f"d(dT/T)/dJ track={chrom['slope_time']:.6e} path={chrom['slope_path']:.6e} "
        f"symplectic -2*pi*xi/C={chrom['predicted_time']:.6e}"
    )
    if abs(chrom["predicted_time"]) > 0:
        lines.append(
            f"time/symplectic={chrom['slope_time'] / chrom['predicted_time']:.6f} "
            f"path/symplectic={chrom['slope_path'] / chrom['predicted_time']:.6f}"
        )
    for j, rt, rc in chrom["phases"]:
        lines.append(f"  J={j:.3e} dT/T={rt:.6e} dC/C={rc:.6e}")
    return pack_report(lines)


def sector_matrix_check():
    """One cylindrical sector, transfer of (x, x', time) against the linear model."""
    lines = ["== single cylindrical sector =="]
    alpha = THETA
    length = RHO * alpha
    k_x = (2.0 - BETA2) / RHO**2
    sqrt_k = math.sqrt(k_x)
    psi = sqrt_k * length
    # Linear map for a particle that enters the field already, coordinates
    # relative to the reference arc. Time factor v_s*dt/ds = 1 - delta/gamma^2
    # + (2-beta^2)*x/rho, integrated on x = R16-like orbit.
    # Numerical: entrance at angle 0, center at origin, reference at (rho, 0)
    # moving +Z. Entrance plane is the X axis (normal +Z). We start INSIDE,
    # so no entrance kick: the particle is created with the interior momentum.
    center = np.array([0.0, 0.0])

    def launch(x, angle, delta):
        r = RHO + x
        # Direction: tangent (0, 1) plus a small angle. angle>0 points outward? 
        # x' = dx/ds, s along the arc. Positive x outward, positive x' increases x.
        # Tangent CCW is +Z at (r, 0). Outward is +X. A positive angle toward +X
        # means velocity has a +X component.
        p_out = p_of_delta(delta)
        # Outside momentum, then entrance kick into the field.
        direction = np.array([math.sin(angle), math.cos(angle)])
        mom = p_out * direction
        pos = np.array([r, 0.0])
        mom = fringe_kick(pos, mom, center, entering=True)
        return pos, mom

    def transfer(x, angle, delta):
        pos, mom = launch(x, angle, delta)
        alpha1 = alpha
        pos, mom, time, path, _ = track_bend(pos, mom, center, 0.0, alpha1, magnetic=False)
        mom = fringe_kick(pos, mom, center, entering=False)
        # Exit frame: reference exit at angle `alpha`, position rho*(cos alpha, sin alpha),
        # tangent (-sin alpha, cos alpha), outward (cos alpha, sin alpha).
        outward = np.array([math.cos(alpha1), math.sin(alpha1)])
        tangent = np.array([-math.sin(alpha1), math.cos(alpha1)])
        rel = pos - RHO * outward
        x_out = float(np.dot(rel, outward))
        p_t = float(np.dot(mom, tangent))
        p_o = float(np.dot(mom, outward))
        angle_out = math.atan2(p_o, p_t)
        return x_out, angle_out, time, path

    # Reference time.
    _, _, t_ref, _ = transfer(0.0, 0.0, 0.0)
    # Columns by finite difference.
    eps_x = 1e-5
    eps_a = 1e-6
    eps_d = 1e-5
    base = transfer(0.0, 0.0, 0.0)
    dx = transfer(eps_x, 0.0, 0.0)
    da = transfer(0.0, eps_a, 0.0)
    dd = transfer(0.0, 0.0, eps_d)
    r11 = (dx[0] - base[0]) / eps_x
    r21 = (dx[1] - base[1]) / eps_x
    r12 = (da[0] - base[0]) / eps_a
    r22 = (da[1] - base[1]) / eps_a
    r16 = (dd[0] - base[0]) / eps_d
    r26 = (dd[1] - base[1]) / eps_d
    # Time relative to reference, in length v_s*Delta t. Positive means later.
    r51 = V_S * (dx[2] - base[2]) / eps_x
    r52 = V_S * (da[2] - base[2]) / eps_a
    r56 = V_S * (dd[2] - base[2]) / eps_d

    # Analytic linear map. Drive F = (2-beta^2)/rho, focus k.
    # x = (F/k) (1-cos psi) delta = rho (1-cos) delta, since F/k = rho.
    c, s = math.cos(psi), math.sin(psi)
    a_r11 = c
    a_r12 = s / sqrt_k
    a_r21 = -sqrt_k * s
    a_r22 = c
    a_r16 = RHO * (1.0 - c)
    a_r26 = RHO * sqrt_k * s
    # Time: v_s dt/ds = 1 - delta/gamma^2 + (2-beta^2) x/rho
    # along the first-order orbit. For entrance (x,x')=(0,0):
    # x(s) = (F/k)*(1-cos(sqrt(k) s)) delta = rho (1-cos) delta
    # integral_0^L (2-beta^2) x/rho ds = (2-beta^2) * rho * delta * (L - sin/sqrt_k)
    # wait x/rho = (1-cos) delta, integral (1-cos) ds = L - sin(psi)/sqrt_k
    # times (2-beta^2): 
    geo = (2.0 - BETA2) * (length - s / sqrt_k)
    kin = length / GAMMA**2
    a_r56 = -kin + geo  # v_s Delta t per unit delta: later if geo wins
    # R51: constant offset x, the homogeneous orbit x(s)= x0 cos
    # integral (2-beta^2)/rho * x0 cos(sqrt k s) ds = (2-beta^2)/rho * x0 * sin/sqrt_k
    a_r51 = (2.0 - BETA2) / RHO * s / sqrt_k
    lines.append(f"psi={psi:.6f} sqrt(k)L, k={k_x:.6e}, L={length:.4f}")
    lines.append(f"t_ref={t_ref:.8f} L/v={length / V_S:.8f}")
    lines.append(f"R11 num={r11:.6f} ana={a_r11:.6f}")
    lines.append(f"R12 num={r12:.6f} ana={a_r12:.6f}")
    lines.append(f"R21 num={r21:.6f} ana={a_r21:.6f}")
    lines.append(f"R22 num={r22:.6f} ana={a_r22:.6f}")
    lines.append(f"R16 num={r16:.6f} ana={a_r16:.6f}")
    lines.append(f"R26 num={r26:.6f} ana={a_r26:.6f}")
    lines.append(f"R51=v dt/dx num={r51:.6f} ana={a_r51:.6f}  (magnet would be {s / (sqrt_k * RHO) * 1:.6f} if h=1/rho only)")
    lines.append(f"R56=v dt/dd num={r56:.6f} ana={a_r56:.6f}")
    # Wollnik-style geometric R56 uses coefficient 1 rather than (2-beta^2), plus L/gamma^2
    # with the opposite sign convention (path ahead). Our later-time R56:
    wollnik_geo = 1.0 * (length - s / sqrt_k)
    lines.append(
        f"if the x-coefficient were 1 (geometry only): R56={-kin + wollnik_geo:.6f}"
    )
    lines.append(
        f"R51 geometry-only would be {s / (RHO * sqrt_k):.6f}; "
        f"ratio ana/geometry={(2-BETA2):.6f}"
    )
    return pack_report(lines)


def main():
    smooth_circle_times()
    print()
    sector_matrix_check()
    print()

    elements, length = build_ring(K_QUAD, sextupole=0.0, magnetic_quad=False)
    info_e = dispersion_alpha(elements, magnetic_bend=False, circumference=length)
    chrom_e = chromaticity_and_action(elements, False, info_e)
    formula_block(info_e, chrom_e, "electrostatic FODO, cylindrical deflectors, electric quads")
    print()

    elements_m, length_m = build_ring(K_QUAD, sextupole=0.0, magnetic_quad=True)
    info_m = dispersion_alpha(elements_m, magnetic_bend=True, circumference=length_m)
    chrom_m = chromaticity_and_action(elements_m, True, info_m)
    formula_block(info_m, chrom_m, "magnetic FODO, uniform dipoles, magnetic quads")
    print()

    # Sextupole scan on the electrostatic ring, a few points around the natural xi.
    print("== electrostatic sextupole scan (magnetic thin sextupole on F planes) ==")
    for s in (-0.15, -0.05, 0.0, 0.05, 0.15):
        el, leng = build_ring(K_QUAD, sextupole=s, magnetic_quad=False)
        try:
            info = dispersion_alpha(el, False, leng)
            chrom = chromaticity_and_action(el, False, info)
        except Exception as exc:  # noqa: BLE001
            print(f"S={s:.3f} FAILED {exc}")
            continue
        pred = -2.0 * math.pi * chrom["xi"] / info["C_s"]
        print(
            f"S={s:.3f} xi={chrom['xi']:.4f} eta0={info['eta0']:.5f} eta1={info['eta1']:.5f} "
            f"slopeT={chrom['slope_time']:.4e} slopeC={chrom['slope_path']:.4e} "
            f"sym={pred:.4e} Lphi_over_J={chrom['slope_time'] - chrom['slope_path']:.4e}"
        )


if __name__ == "__main__":
    main()
