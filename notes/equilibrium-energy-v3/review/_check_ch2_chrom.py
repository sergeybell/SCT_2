# Checks for chapter 2 chromaticity / path lengthening.
# (a) oint x_nl/rho = oint <F^x> D_1
# (b) Lambda from path length vs -2*pi*J*xi/C, xi = d(nu)/d(delta)
# Thin-lens FODO with thin dipoles, and a sector-dipole comparison.
import math

# ---------- linear algebra, 2x2 ----------

def mm(A, B):
    return (
        (A[0][0] * B[0][0] + A[0][1] * B[1][0],
         A[0][0] * B[0][1] + A[0][1] * B[1][1]),
        (A[1][0] * B[0][0] + A[1][1] * B[1][0],
         A[1][0] * B[0][1] + A[1][1] * B[1][1]),
    )

def mv(A, v):
    return (A[0][0] * v[0] + A[0][1] * v[1],
            A[1][0] * v[0] + A[1][1] * v[1])

def mdrift(L):
    return ((1.0, L), (0.0, 1.0))

def mquad(k):
    # Delta x' = -k x, k = integral K_1 ds; k>0 focuses horizontally
    return ((1.0, 0.0), (-k, 1.0))

def msector(rho, theta):
    # Linear body of a sector dipole, K = 1/rho**2, independent of delta.
    c = math.cos(theta)
    s = math.sin(theta)
    return ((c, rho * s), (-s / rho, c))

I2 = ((1.0, 0.0), (0.0, 1.0))

def tune_from_matrix(M):
    tr = M[0][0] + M[1][1]
    if tr > 2.0 and tr < 2.0 + 1e-7:
        return 0.0
    if tr < -2.0 and tr > -2.0 - 1e-7:
        return 0.5
    if abs(tr) > 2.0 + 1e-8:
        raise RuntimeError("unstable trace %s" % tr)
    tr = max(-2.0, min(2.0, tr))
    return math.acos(tr / 2.0) / (2.0 * math.pi)

def twiss_from_matrix(M):
    mu = 2.0 * math.pi * tune_from_matrix(M)
    sm = math.sin(mu)
    if abs(sm) < 1e-14:
        raise RuntimeError("integer tune, beta undefined")
    beta = M[0][1] / sm
    alpha = (M[0][0] - M[1][1]) / (2.0 * sm)
    gamma = (1.0 + alpha * alpha) / beta
    return beta, alpha, gamma, mu

# ============================================================
# Thin FODO: curvature only in zero-length dipoles (no 1/rho**2)
# ============================================================

def build_thin_cell(k_f, k_d, L, theta, m_f, m_d):
    # QF, drift, bend, drift, sext at QD?  sextupoles sit on the quads.
    # Order: QF+SF, drift L/2, bend, drift L/2, QD+SD, drift L/2, bend, drift L/2
    half = 0.5 * L
    return [
        ("quad", k_f),
        ("sext", m_f),
        ("drift", half),
        ("dip", theta),
        ("drift", half),
        ("quad", k_d),
        ("sext", m_d),
        ("drift", half),
        ("dip", theta),
        ("drift", half),
    ]


def oneturn_thin(cell, delta):
    """Linear map about the closed orbit, strengths scaled by 1/(1+delta)."""
    scale = 1.0 / (1.0 + delta)
    # First find the periodic dispersion-like orbit including the
    # quadratic sextupole kick, by iteration. Dipole kick = theta*delta*scale.
    x, xp = 0.0, 0.0
    for _ in range(25):
        x, xp, _ = track_thin_orbit(cell, x, xp, delta)
        # fixed-point residual is removed by Newton below
        break
    # Newton on the horizontal fixed point (y=0)
    def resid(v):
        x2, xp2, _ = track_thin_orbit(cell, v[0], v[1], delta)
        return (x2 - v[0], xp2 - v[1])
    v = (0.0, 0.0)
    for _ in range(12):
        r = resid(v)
        # Jacobian
        eps = 1e-8
        r1 = resid((v[0] + eps, v[1]))
        r2 = resid((v[0], v[1] + eps))
        a11 = (r1[0] - r[0]) / eps - 1.0  # wait, resid already subtracts v, so Jacobian of resid
        # resid(v+eps) - resid(v) is the derivative of (map - id)
        a11 = (r1[0] - r[0]) / eps
        a12 = (r2[0] - r[0]) / eps
        a21 = (r1[1] - r[1]) / eps
        a22 = (r2[1] - r[1]) / eps
        det = a11 * a22 - a12 * a21
        dx = (-r[0] * a22 + r[1] * a12) / det
        dxp = (-a11 * r[1] + a21 * r[0]) / det
        v = (v[0] + dx, v[1] + dxp)
        if abs(r[0]) + abs(r[1]) < 1e-14:
            break
    xco, xpco = v
    # One-turn matrices linearized about the orbit.
    # Quad k -> k/(1+delta); sextupole feed-down k_extra = m*x/(1+delta)
    Mx = I2
    My = I2
    # Rebuild while propagating the known orbit, accumulating Jacobian.
    x, xp = xco, xpco
    y_ignored = 0.0
    for el in cell:
        kind = el[0]
        if kind == "drift":
            L = el[1]
            Mx = mm(mdrift(L), Mx)
            My = mm(mdrift(L), My)
            x, xp = x + L * xp, xp
        elif kind == "quad":
            k = el[1] * scale
            Mx = mm(mquad(k), Mx)
            My = mm(mquad(-k), My)  # K_y = -K_1
            xp = xp - k * x
        elif kind == "sext":
            m = el[1] * scale
            # about y=0: Delta xp = -(m/2)*(x^2) ; d(Delta xp)/dx = -m*x
            k_extra = m * x
            Mx = mm(mquad(k_extra), Mx)
            # vertical: Delta yp = m*x*y => Delta yp' contribution k_y extra = -m*x
            # equation y'' = +K2*x*y, so y'' + K_y y = 0 with K_y = -K2*x
            # thin: Delta y' = + (m*x) * y, which is mquad(k) with k = -m*x
            My = mm(mquad(-m * x), My)
            xp = xp - 0.5 * m * (x * x)
        elif kind == "dip":
            th = el[1]
            xp = xp + th * delta * scale
        else:
            raise RuntimeError(kind)
    return Mx, My, xco, xpco


def track_thin_orbit(cell, x, xp, delta):
    """Nonlinear horizontal track, y=0. Returns x, xp, path extra over reference."""
    scale = 1.0 / (1.0 + delta)
    path = 0.0
    for el in cell:
        kind = el[0]
        if kind == "drift":
            L = el[1]
            path += L * (math.sqrt(1.0 + xp * xp) - 1.0)
            x = x + L * xp
        elif kind == "quad":
            xp = xp - el[1] * scale * x
        elif kind == "sext":
            m = el[1] * scale
            xp = xp - 0.5 * m * (x * x)
        elif kind == "dip":
            th = el[1]
            path += x * th
            xp = xp + th * delta * scale
        else:
            raise RuntimeError(kind)
    return x, xp, path


def track_thin_4d(cell, x, xp, y, yp, delta):
    scale = 1.0 / (1.0 + delta)
    path = 0.0
    for el in cell:
        kind = el[0]
        if kind == "drift":
            L = el[1]
            path += L * (math.sqrt(1.0 + xp * xp + yp * yp) - 1.0)
            x = x + L * xp
            y = y + L * yp
        elif kind == "quad":
            k = el[1] * scale
            xp = xp - k * x
            yp = yp + k * y
        elif kind == "sext":
            m = el[1] * scale
            x2 = x * x
            y2 = y * y
            xp = xp - 0.5 * m * (x2 - y2)
            yp = yp + m * x * y
        elif kind == "dip":
            th = el[1]
            path += x * th
            xp = xp + th * delta * scale
        else:
            raise RuntimeError(kind)
    return x, xp, y, yp, path


def linear_optics_thin(cell):
    """Twiss, dispersion and gamma integral at delta=0, sextupoles off
    (feed-down vanishes on x=0). Sextupole locations still recorded.
    """
    Mx = I2
    My = I2
    # dispersion transport: (D, Dp) -> M (D,Dp) + (0, theta) at dipoles
    # Accumulate the inhomogeneous vector by propagating kicks to the end.
    # We'll do a second pass once the one-turn matrix is known.
    for el in cell:
        kind = el[0]
        if kind == "drift":
            Md = mdrift(el[1])
            Mx = mm(Md, Mx)
            My = mm(Md, My)
        elif kind == "quad":
            Mx = mm(mquad(el[1]), Mx)
            My = mm(mquad(-el[1]), My)
        elif kind == "sext":
            pass
        elif kind == "dip":
            pass
        else:
            raise RuntimeError(kind)
    bx, ax, gx, mux = twiss_from_matrix(Mx)
    by, ay, gy, muy = twiss_from_matrix(My)

    # Propagate Twiss and dispersion from the start.
    def prop_twiss(beta, alpha, gamma, el_kind, val):
        if el_kind == "drift":
            L = val
            b2 = beta - 2 * alpha * L + gamma * L * L
            a2 = alpha - gamma * L
            g2 = gamma  # K=0
            return b2, a2, g2, gamma * L  # gamma constant, integral gamma*L
        if el_kind == "quad":
            # alpha jumps by +beta*k?  x' := x' - k x
            # beta continuous, alpha_after = alpha_before + beta*k? 
            # Transport of (beta, alpha): for a thin lens Delta x'=-k x,
            # alpha_+ = alpha_- + k*beta? Standard: alpha_exit = alpha_entrance + beta*k
            # Check: M = [[1,0],[-k,1]], Twiss transform.
            # Convention alpha = -beta'/2. A focusing lens (k>0) makes beta
            # start decreasing if we are at a maximum? At a waist before a
            # focusing lens alpha=0, after the lens the beam converges, beta
            # decreases, so beta'<0, alpha>0. alpha increases. Yes alpha += k*beta.
            return beta, alpha + val * beta, (1.0 + (alpha + val * beta) ** 2) / beta, 0.0
        if el_kind == "sext" or el_kind == "dip":
            return beta, alpha, gamma, 0.0
        raise RuntimeError(el_kind)

    # Dispersion: integrate through the cell with periodic closure.
    # Build one-turn M and the kick vector reduced to the start.
    # Propagate an initial (D,Dp)=(0,0), record the result V = M*0 + kicks = inhomogeneous.
    D, Dp = 0.0, 0.0
    Minhom = I2  # actually track explicitly
    for el in cell:
        kind = el[0]
        if kind == "drift":
            D = D + el[1] * Dp
        elif kind == "quad":
            Dp = Dp - el[1] * D
        elif kind == "sext":
            pass
        elif kind == "dip":
            Dp = Dp + el[1]
        else:
            raise RuntimeError(kind)
    V = (D, Dp)  # image of 0
    # (I - M) D0 = V, but V = M_partial kicks... 
    # Fixed point: D_end = M D_start + V_from_zero_start, and D_end=D_start,
    # so (I - M) D_start = V only if V is the result of propagating 0.
    # Yes: map(D) = M D + map(0), map(0)=V, so D = M D + V, (I-M)D = V.
    Mx_full = Mx
    a11 = 1.0 - Mx_full[0][0]
    a12 = -Mx_full[0][1]
    a21 = -Mx_full[1][0]
    a22 = 1.0 - Mx_full[1][1]
    det = a11 * a22 - a12 * a21
    D0 = (V[0] * a22 - V[1] * a12) / det
    Dp0 = (-V[0] * a21 + V[1] * a11) / det

    # Walk again, record beta, D at sextupoles and quads, integrate gamma and D/rho.
    bx_s, ax_s, gx_s = bx, ax, gx
    by_s, ay_s, gy_s = by, ay, gy
    D, Dp = D0, Dp0
    oint_gx = 0.0
    oint_gy = 0.0
    sum_k_bx = 0.0
    sum_k_by = 0.0
    sum_xnl_prep = []  # (beta_x, beta_y, D, m) at sextupoles; m filled by caller
    sext_idx = 0
    dip_x_weight = []  # (theta, D) and we'll need x_nl; store positions by re-propagating later
    # Store element-resolved D and beta by simulating and saving state at each element entrance
    states = []
    for el in cell:
        states.append((el, bx_s, ax_s, gx_s, by_s, ay_s, gy_s, D, Dp))
        kind = el[0]
        if kind == "drift":
            L = el[1]
            oint_gx += gx_s * L
            oint_gy += gy_s * L
            bx_s, ax_s, gx_s, _ = prop_twiss(bx_s, ax_s, gx_s, "drift", L)
            by_s, ay_s, gy_s, _ = prop_twiss(by_s, ay_s, gy_s, "drift", L)
            D = D + L * Dp
        elif kind == "quad":
            k = el[1]
            sum_k_bx += k * bx_s
            sum_k_by += k * by_s  # integral K_1 beta_y ; K_y = -K_1 so oint K_y beta_y = -sum k beta_y
            bx_s, ax_s, gx_s, _ = prop_twiss(bx_s, ax_s, gx_s, "quad", k)
            # vertical lens strength is -k (defocusing if k>0). alpha_y += k_y * beta = -k*beta
            by_s, ay_s, gy_s, _ = prop_twiss(by_s, ay_s, gy_s, "quad", -k)
            Dp = Dp - k * D
        elif kind == "sext":
            sum_xnl_prep.append((bx_s, by_s, D, el[1]))
            sext_idx += 1
        elif kind == "dip":
            dip_x_weight.append(el[1])  # theta; x_nl computed later at this D
            # D continuous, Dp jumps after we record entrance D. Record D now.
            # We'll pair below.
            Dp = Dp + el[1]
        else:
            raise RuntimeError(kind)

    # Recompute D at dipoles in order
    D, Dp = D0, Dp0
    D_at_dip = []
    for el in cell:
        kind = el[0]
        if kind == "drift":
            D = D + el[1] * Dp
        elif kind == "quad":
            Dp = Dp - el[1] * D
        elif kind == "sext":
            pass
        elif kind == "dip":
            D_at_dip.append(D)
            Dp = Dp + el[1]
    return {
        "nu_x": mux / (2 * math.pi),
        "nu_y": muy / (2 * math.pi),
        "beta_x0": bx,
        "alpha_x0": ax,
        "gamma_x0": gx,
        "beta_y0": by,
        "alpha_y0": ay,
        "gamma_y0": gy,
        "D0": D0,
        "Dp0": Dp0,
        "oint_gx": oint_gx,
        "oint_gy": oint_gy,
        "sum_k_bx": sum_k_bx,
        "sum_k_by": sum_k_by,
        "sext": sum_xnl_prep,  # (bx, by, D, m)
        "D_dip": D_at_dip,
        "theta_dip": dip_x_weight,
        "Mx": Mx,
        "My": My,
    }


def closed_orbit_kicks(cell, kicks_at_sext, n_sext):
    """Periodic orbit from thin horizontal kicks applied at sextupole slots.
    kicks_at_sext[i] is Delta x' at sextupole i (the averaged force integral).
    Returns x at each dipole and at each sextupole.
    """
    # map of zero + record
    # Same (I-M) method but kicks are at sextupoles instead of dipoles.
    x, xp = 0.0, 0.0
    isext = 0
    for el in cell:
        kind = el[0]
        if kind == "drift":
            x = x + el[1] * xp
        elif kind == "quad":
            xp = xp - el[1] * x
        elif kind == "sext":
            xp = xp + kicks_at_sext[isext]
            isext += 1
        elif kind == "dip":
            pass
    V = (x, xp)
    # one-turn M without kicks (quads + drifts)
    M = I2
    for el in cell:
        kind = el[0]
        if kind == "drift":
            M = mm(mdrift(el[1]), M)
        elif kind == "quad":
            M = mm(mquad(el[1]), M)
    a11 = 1.0 - M[0][0]
    a12 = -M[0][1]
    a21 = -M[1][0]
    a22 = 1.0 - M[1][1]
    det = a11 * a22 - a12 * a21
    x0 = (V[0] * a22 - V[1] * a12) / det
    xp0 = (-V[0] * a21 + V[1] * a11) / det
    # replay
    x, xp = x0, xp0
    isext = 0
    x_dip = []
    x_sext = []
    for el in cell:
        kind = el[0]
        if kind == "drift":
            x = x + el[1] * xp
        elif kind == "quad":
            xp = xp - el[1] * x
        elif kind == "sext":
            x_sext.append(x)
            xp = xp + kicks_at_sext[isext]
            isext += 1
        elif kind == "dip":
            x_dip.append(x)
    return x_dip, x_sext, x0, xp0


def xi_numeric_thin(cell, eps=1e-6):
    nup = tune_from_matrix(oneturn_thin(cell, eps)[0])
    num = tune_from_matrix(oneturn_thin(cell, -eps)[0])
    nup_y = tune_from_matrix(oneturn_thin(cell, eps)[1])
    num_y = tune_from_matrix(oneturn_thin(cell, -eps)[1])
    return (nup - num) / (2 * eps), (nup_y - num_y) / (2 * eps)


def lambda_direct_thin(cell, opt, Jx, Jy, nphase=16):
    """Path length averaged over the betatron phase.

    The stored beam is centred on the periodic orbit x_nl generated by the
    phase-averaged sextupole force. Particles are launched from that orbit.
    The sextupole is kept nonlinear; the average kick is already contained
    in x_nl, and the piece oscillating at 2*psi averages out of the length
    at order J.
    """
    bx, ax = opt["beta_x0"], opt["alpha_x0"]
    by, ay = opt["beta_y0"], opt["alpha_y0"]
    kicks = []
    for bsx, bsy, D, m in opt["sext"]:
        kicks.append(-0.5 * m * (Jx * bsx - Jy * bsy))
    _xd, _xs, x_nl, xp_nl = closed_orbit_kicks(cell, kicks, len(kicks))
    acc = 0.0
    for i in range(nphase):
        psi = 2.0 * math.pi * i / nphase
        psi_y = 2.0 * math.pi * ((i * 3) % nphase) / nphase
        x = x_nl + (math.sqrt(2 * Jx * bx) * math.cos(psi) if Jx > 0 else 0.0)
        xp = xp_nl + (
            -math.sqrt(2 * Jx / bx) * (math.sin(psi) + ax * math.cos(psi)) if Jx > 0 else 0.0
        )
        y = math.sqrt(2 * Jy * by) * math.cos(psi_y) if Jy > 0 else 0.0
        yp = (
            -math.sqrt(2 * Jy / by) * (math.sin(psi_y) + ay * math.cos(psi_y)) if Jy > 0 else 0.0
        )
        _x2, _xp2, _y2, _yp2, path = track_thin_4d(cell, x, xp, y, yp, 0.0)
        acc += path
    return acc / nphase


def run_thin():
    k_f = 0.45
    k_d = -0.40
    L = 2.0
    theta = 0.15
    C = 2.0 * L  # two half+half drifts of length L, quads and dips have zero length
    # C = 4 * (L/2) = 2 L = 4.0
    print("--- thin FODO, C =", 2 * L)
    cases = [
        ("no sext", 0.0, 0.0),
        ("sext on", 1.2, -1.5),
    ]
    for name, m_f, m_d in cases:
        cell = build_thin_cell(k_f, k_d, L, theta, m_f, m_d)
        opt = linear_optics_thin(cell)
        # (a) identity
        sext = opt["sext"]
        Jx, Jy = 1.0e-6, 2.0e-6
        kicks = []
        oint_F_D = 0.0
        for bx, by, D, m in sext:
            # integrated <F^x> = -(m/2)*(Jx*bx - Jy*by)
            Fint = -0.5 * m * (Jx * bx - Jy * by)
            kicks.append(Fint)
            oint_F_D += Fint * D
        x_dip, _, _, _ = closed_orbit_kicks(cell, kicks, len(sext))
        oint_x_rho = 0.0
        for xnl, th in zip(x_dip, opt["theta_dip"]):
            oint_x_rho += xnl * th
        print(name)
        print("  (a) oint x_nl/rho", oint_x_rho, "oint <F> D", oint_F_D,
              "diff", oint_x_rho - oint_F_D)

        # gamma vs K beta
        print("  oint gamma_x", opt["oint_gx"], "sum k beta_x", opt["sum_k_bx"],
              "diff", opt["oint_gx"] - opt["sum_k_bx"])
        print("  oint gamma_y", opt["oint_gy"], "-sum k beta_y", -opt["sum_k_by"],
              "diff", opt["oint_gy"] + opt["sum_k_by"])

        # xi from integral and from tune
        # xi_x = -1/4pi * (sum k bx - sum m D bx)
        # xi_y = +1/4pi * (sum k by - sum m D by)
        sum_mD_bx = 0.0
        sum_mD_by = 0.0
        for bx, by, D, m in sext:
            sum_mD_bx += m * D * bx
            sum_mD_by += m * D * by
        xi_x_int = -(opt["sum_k_bx"] - sum_mD_bx) / (4.0 * math.pi)
        xi_y_int = (opt["sum_k_by"] - sum_mD_by) / (4.0 * math.pi)
        xi_x_num, xi_y_num = xi_numeric_thin(cell)
        print("  xi_x int", xi_x_int, "num", xi_x_num, "diff", xi_x_int - xi_x_num)
        print("  xi_y int", xi_y_int, "num", xi_y_num, "diff", xi_y_int - xi_y_num)

        # Lambda from slopes + x_nl
        C_s = 2.0 * L
        Lam_slopes = (0.5 * Jx * opt["oint_gx"] + 0.5 * Jy * opt["oint_gy"] + oint_x_rho) / C_s
        Lam_xi = -(2.0 * math.pi / C_s) * (Jx * xi_x_num + Jy * xi_y_num)
        # direct tracking. Use small actions so O(J^2) is small; report Lambda/J scale.
        # For a clean O(J) comparison, subtract the quadratic contamination by
        # using two amplitudes.
        def lam_of(sx, sy):
            return lambda_direct_thin(cell, opt, sx, sy) / C_s
        # Richardson: Lambda = a1*J + a2*J^2, use J and J/4? actions scale.
        # path is even in amplitude only if no sextupole. With sextupole, odd powers
        # can appear? x^2 kick is even in (x -> -x) together with a phase shift of pi,
        # and path should be even in the amplitude. So Lambda = c2 A^2 + c4 A^4 = c1 J + c2 J^2.
        Jx1, Jy1 = 1.0e-6, 2.0e-6
        Jx2, Jy2 = 0.5e-6, 1.0e-6
        L1 = lam_of(Jx1, Jy1)
        L2 = lam_of(Jx2, Jy2)
        # L = a*J + b*J^2, with J scaled by s and s^2: let s2 = 1/2
        # L1 = a*J + b*J^2, L2 = a*(J/2) + b*(J/2)^2 = a J/2 + b J^2/4
        # 2*L2 = a J + b J^2/2
        # 4*L2 - L1 = a J + b J^2 =?  2*L2 - 0.5*L1 = a*(J/2)
        aJ = 2.0 * L2 - 0.5 * L1   # this is a*(J/2)? Let's solve properly.
        # L1 = a + b   in units where we set u = J-scale of case 1, but Jx and Jy
        # scale together by factor f: L(f) = A f + B f^2
        # L1 = A + B, L2 = A/2 + B/4
        # 4 L2 = 2 A + B
        # 4 L2 - L1 = A
        A = 4.0 * L2 - L1  # this is Lambda linear in the (Jx1,Jy1) actions
        print("  Lambda slopes+xnl", Lam_slopes)
        print("  Lambda -2pi J xi ", Lam_xi)
        print("  Lambda direct O(J)", A, "direct raw", L1)
        print("  rel (direct-xi)/xi", (A - Lam_xi) / Lam_xi if Lam_xi else None)
        print("  rel (slopes-xi)/xi", (Lam_slopes - Lam_xi) / Lam_xi if Lam_xi else None)
        print("  nu", opt["nu_x"], opt["nu_y"], "D0", opt["D0"])



# ============================================================
# Sector dipole: exact circular arc in a hard-edge sector.
# Body focusing 1/rho**2 is in the geometry; there is no fringe integral.
# ============================================================

def track_sector(x, xp, y, yp, delta, rho, theta):
    """Uniform-field sector. Horizontal motion is an exact circular arc.
    Vertical motion is a drift of length rho*theta (no edge focusing).
    Returns x, xp, y, yp and the absolute path length.
    """
    R = rho * (1.0 + delta)
    if R <= 0.0 or 1.0 + x / rho <= 0.0:
        raise RuntimeError("non-physical sector coordinates")
    phi = math.atan2(xp, 1.0 + x / rho)
    X = rho + x
    Z = 0.0
    vx = math.sin(phi)
    vz = math.cos(phi)
    cx = X - R * vz
    cz = Z + R * vx
    ang0 = math.atan2(-cz, X - cx)

    def lab_and_deriv(ang):
        X1 = cx + R * math.cos(ang)
        Z1 = cz + R * math.sin(ang)
        r2 = X1 * X1 + Z1 * Z1
        dX = -R * math.sin(ang)
        dZ = R * math.cos(ang)
        dlab = (X1 * dZ - Z1 * dX) / r2
        lab = math.atan2(Z1, X1)
        while lab - theta > math.pi:
            lab -= 2.0 * math.pi
        while theta - lab > math.pi:
            lab += 2.0 * math.pi
        return lab, dlab, X1, Z1

    ang = ang0 + theta
    for _ in range(30):
        lab, dlab, _, _ = lab_and_deriv(ang)
        dang = (theta - lab) / dlab
        ang += dang
        if abs(dang) < 1e-14:
            break
    else:
        raise RuntimeError("sector exit not found")
    lab, _, X2, Z2 = lab_and_deriv(ang)
    if abs(lab - theta) > 1e-9:
        raise RuntimeError("sector exit angle %s" % lab)
    r2 = math.hypot(X2, Z2)
    x2 = r2 - rho
    vx2 = -(Z2 - cz) / R
    vz2 = (X2 - cx) / R
    cth = math.cos(theta)
    sth = math.sin(theta)
    v_rad = vx2 * cth + vz2 * sth
    v_tan = -vx2 * sth + vz2 * cth
    xp2 = (1.0 + x2 / rho) * (v_rad / v_tan)
    path_h = R * (ang - ang0)
    y2 = y + (rho * theta) * yp
    yp2 = yp
    path = path_h + 0.5 * (yp * yp) * (rho * theta)
    return x2, xp2, y2, yp2, path


def sector_matrix_numeric(rho, theta):
    e = 1e-6
    def f(x, xp):
        x2, xp2, _, _, _ = track_sector(x, xp, 0.0, 0.0, 0.0, rho, theta)
        return x2, xp2
    xp0, xpp0 = f(0.0, 0.0)
    xa, xpa = f(e, 0.0)
    xb, xpb = f(-e, 0.0)
    xc, xpc = f(0.0, e)
    xd, xpd = f(0.0, -e)
    M = (
        ((xa - xb) / (2 * e), (xc - xd) / (2 * e)),
        ((xpa - xpb) / (2 * e), (xpc - xpd) / (2 * e)),
    )
    return M, xp0, xpp0


def check_sector_matrix():
    rho, theta = 5.0, 0.3
    M, x0, xp0 = sector_matrix_numeric(rho, theta)
    A = msector(rho, theta)
    diffs = (
        M[0][0] - A[0][0], M[0][1] - A[0][1],
        M[1][0] - A[1][0], M[1][1] - A[1][1],
    )
    dlt = 1e-6
    xp, xpp, _, _, _ = track_sector(0.0, 0.0, 0.0, 0.0, dlt, rho, theta)
    xm, xpm, _, _, _ = track_sector(0.0, 0.0, 0.0, 0.0, -dlt, rho, theta)
    Dx = (xp - xm) / (2 * dlt)
    Dxp = (xpp - xpm) / (2 * dlt)
    print("--- sector map vs linear body")
    print("  matrix diff", diffs, "exit of reference", x0, xp0)
    print("  dispersion transport", Dx, Dxp,
          "analytic", rho * (1 - math.cos(theta)), math.sin(theta))


def ring_uniform(n_sectors, rho, J):
    theta = 2.0 * math.pi / n_sectors
    C = 2.0 * math.pi * rho
    nphase = 16
    extras = []
    for i in range(nphase):
        psi = 2.0 * math.pi * i / nphase
        x = math.sqrt(2.0 * J * rho) * math.cos(psi)
        xp = -math.sqrt(2.0 * J / rho) * math.sin(psi)
        path = 0.0
        for _ in range(n_sectors):
            x, xp, _, _, p = track_sector(x, xp, 0.0, 0.0, 0.0, rho, theta)
            path += p
        extras.append(path - C)
    Lam = sum(extras) / nphase / C
    # fractional tune about the concentric closed orbit
    def matrix_at(delta):
        def once(x, xp):
            for _ in range(n_sectors):
                x, xp, _, _, _ = track_sector(x, xp, 0.0, 0.0, delta, rho, theta)
            return x, xp
        v = [rho * delta, 0.0]
        for _ in range(8):
            x2, xp2 = once(v[0], v[1])
            r0, r1 = x2 - v[0], xp2 - v[1]
            e = 1e-8
            xa, xpa = once(v[0] + e, v[1])
            xb, xpb = once(v[0], v[1] + e)
            a11 = (xa - x2) / e - 1.0
            a12 = (xb - x2) / e
            a21 = (xpa - xp2) / e
            a22 = (xpb - xp2) / e - 1.0
            det = a11 * a22 - a12 * a21
            v[0] += (-r0 * a22 + r1 * a12) / det
            v[1] += (-a11 * r1 + a21 * r0) / det
            if abs(r0) + abs(r1) < 1e-13:
                break
        e = 1e-7
        x0, xp0 = once(v[0], v[1])
        xa, xpa = once(v[0] + e, v[1])
        xb, xpb = once(v[0] - e, v[1])
        xc, xpc = once(v[0], v[1] + e)
        xd, xpd = once(v[0], v[1] - e)
        M = (
            ((xa - xb) / (2 * e), (xc - xd) / (2 * e)),
            ((xpa - xpb) / (2 * e), (xpc - xpd) / (2 * e)),
        )
        return M, v
    eps = 1e-4
    Mp, vp = matrix_at(eps)
    Mm, vm = matrix_at(-eps)
    # Trace is 2: the tune is an integer, and acos is unstable. Report |M-I|.
    def dev(M):
        return max(abs(M[0][0] - 1.0), abs(M[0][1]), abs(M[1][0]), abs(M[1][1] - 1.0))
    Lam_slope = (0.5 * J * (C / rho)) / C  # gamma = 1/rho, oint gamma = C/rho
    print("--- uniform sector ring rho", rho, "N", n_sectors, "J", J)
    print("  orbit x", vp[0], "expected", rho * eps)
    print("  |M-I| +delta", dev(Mp), "-delta", dev(Mm), "(xi = 0)")
    print("  Lambda direct", Lam, "slope-only", Lam_slope, "-2pi J xi/C", 0.0)
    print("  path extra min/max", min(extras), max(extras))


def twiss_transport(M, beta, alpha, gamma):
    C, S = M[0][0], M[0][1]
    Cp, Sp = M[1][0], M[1][1]
    b2 = C * C * beta - 2 * C * S * alpha + S * S * gamma
    a2 = -C * Cp * beta + (C * Sp + S * Cp) * alpha - S * Sp * gamma
    g2 = Cp * Cp * beta - 2 * Cp * Sp * alpha + Sp * Sp * gamma
    return b2, a2, g2


def sector_fodo_linear(rho, theta, k, m_f, m_d, n_cell, n_slice=400):
    """Linear model: quads scale as 1/(1+delta), sector body K=1/rho**2 does not.
    Sextupole feed-down is m*D*delta. Returns optics at delta=0 and xi.
    """
    def sector_pieces():
        dth = theta / n_slice
        return msector(rho, dth), rho * dth

    def oneturn(delta):
        scale = 1.0 / (1.0 + delta)
        # dispersion-like closed orbit, then feed-down. Iterate once: D from delta,
        # feed-down uses x_co.
        # First pass, no feed-down, get orbit.
        def propagate(x, xp, use_feed):
            Mx = I2
            My = I2
            for _c in range(n_cell):
                for plane_sign, m_sext in ((k, m_f), (-k, m_d)):
                    ks = plane_sign * scale
                    if use_feed:
                        ks = ks + (m_sext * scale) * x
                    Mx = mm(mquad(ks), Mx)
                    # vertical lens: quad strength -plane_sign, sextupole K_y extra = -m*x
                    ky = -plane_sign * scale
                    if use_feed:
                        ky = ky - (m_sext * scale) * x
                    My = mm(mquad(ky), My)
                    xp = xp - ks * x
                    Ms, ds = sector_pieces()
                    # inhomogeneous sector: apply slice by slice
                    h = 1.0 / rho
                    for _s in range(n_slice):
                        # x'' + x/rho**2 = delta/rho
                        c = math.cos(theta / n_slice)
                        s = math.sin(theta / n_slice)
                        # map (x, xp, delta)
                        x_new = c * x + (rho * s) * xp + rho * (1.0 - c) * delta
                        xp_new = -(s / rho) * x + c * xp + s * delta
                        x, xp = x_new, xp_new
                        Mx = mm(Ms, Mx)
                        My = mm(mdrift(ds), My)
            return x, xp, Mx, My
        # fixed point without using feed-down inside the unknown; solve by Newton
        # on the nonlinear (feed-down) map. Feed-down depends on x at the sextupole.
        def once(x, xp):
            scale = 1.0 / (1.0 + delta)
            for _c in range(n_cell):
                for plane_sign, m_sext in ((k, m_f), (-k, m_d)):
                    ks = plane_sign * scale + m_sext * scale * x
                    xp = xp - ks * x
                    c = math.cos(theta)
                    s = math.sin(theta)
                    x = c * x + rho * s * xp + rho * (1.0 - c) * delta
                    xp = -(s / rho) * x + c * xp + s * delta
                    # BUG: the second line uses the updated x. Fix below.
            return x, xp
        # The bug note is real: standard order is simultaneous use of entrance x.
        # Rewrite once() correctly.
        def once(x, xp):
            for _c in range(n_cell):
                for plane_sign, m_sext in ((k, m_f), (-k, m_d)):
                    ks = plane_sign * scale + m_sext * scale * x
                    xp = xp - ks * x
                    c = math.cos(theta)
                    s = math.sin(theta)
                    xn = c * x + rho * s * xp + rho * (1.0 - c) * delta
                    xpn = -(s / rho) * x + c * xp + s * delta
                    x, xp = xn, xpn
            return x, xp
        v = [0.0, 0.0]
        for _ in range(12):
            x2, xp2 = once(v[0], v[1])
            r0, r1 = x2 - v[0], xp2 - v[1]
            e = 1e-8
            xa, xpa = once(v[0] + e, v[1])
            xb, xpb = once(v[0], v[1] + e)
            a11 = (xa - x2) / e - 1.0
            a12 = (xb - x2) / e
            a21 = (xpa - xp2) / e
            a22 = (xpb - xp2) / e - 1.0
            det = a11 * a22 - a12 * a21
            v[0] += (-r0 * a22 + r1 * a12) / det
            v[1] += (-a11 * r1 + a21 * r0) / det
            if abs(r0) + abs(r1) < 1e-14:
                break
        # Jacobian at v, quad and sextupole linearized, sector homogeneous
        x, xp = v
        Mx = I2
        My = I2
        bx = by = ax = ay = None
        # We only need the matrix here; Twiss is computed at delta=0 separately.
        for _c in range(n_cell):
            for plane_sign, m_sext in ((k, m_f), (-k, m_d)):
                ks = plane_sign * scale + m_sext * scale * x
                ky = -plane_sign * scale - m_sext * scale * x
                Mx = mm(mquad(ks), Mx)
                My = mm(mquad(ky), My)
                xp = xp - ks * x
                c = math.cos(theta)
                s = math.sin(theta)
                xn = c * x + rho * s * xp + rho * (1.0 - c) * delta
                xpn = -(s / rho) * x + c * xp + s * delta
                x, xp = xn, xpn
                Mx = mm(msector(rho, theta), Mx)
                My = mm(mdrift(rho * theta), My)
        return Mx, My, v

    Mx, My, v = oneturn(0.0)
    nu_x = tune_from_matrix(Mx)
    nu_y = tune_from_matrix(My)
    bx, ax, gx, _ = twiss_from_matrix(Mx)
    by, ay, gy, _ = twiss_from_matrix(My)

    # Integrate gamma and beta/rho**2, and record beta, D at sextupoles, by walking delta=0.
    # Periodic dispersion:
    def disp_once(D, Dp):
        for _c in range(n_cell):
            for plane_sign, _m in ((k, m_f), (-k, m_d)):
                Dp = Dp - plane_sign * D
                c = math.cos(theta)
                s = math.sin(theta)
                Dn = c * D + rho * s * Dp + rho * (1.0 - c)
                Dpn = -(s / rho) * D + c * Dp + s
                D, Dp = Dn, Dpn
        return D, Dp
    # (map - id) fixed point. map(D) = M D + map(0)
    zD, zDp = disp_once(0.0, 0.0)
    # Build M by finite difference
    e = 1e-8
    dD, dDp = disp_once(e, 0.0)
    pD, pDp = disp_once(0.0, e)
    M11 = (dD - zD) / e
    M12 = (pD - zD) / e
    M21 = (dDp - zDp) / e
    M22 = (pDp - zDp) / e
    # (I-M) D = map(0) = (zD, zDp)
    a11, a12, a21, a22 = 1 - M11, -M12, -M21, 1 - M22
    det = a11 * a22 - a12 * a21
    D0 = (zD * a22 - zDp * a12) / det
    Dp0 = (-zD * a21 + zDp * a11) / det

    # Walk Twiss and dispersion, integrate.
    b, a, g = bx, ax, gx
    bY, aY, gY = by, ay, gy
    D, Dp = D0, Dp0
    oint_gx = oint_gy = oint_bx_over_rho2 = 0.0
    sum_k_bx = sum_k_by = 0.0
    sext_rec = []  # bx, by, D, m
    xnl_slots = []
    for _c in range(n_cell):
        for plane_sign, m_sext in ((k, m_f), (-k, m_d)):
            sum_k_bx += plane_sign * b
            sum_k_by += plane_sign * bY
            sext_rec.append((b, bY, D, m_sext))
            # quad transport of Twiss
            b, a, g = b, a + plane_sign * b, (1.0 + (a + plane_sign * b) ** 2) / b
            ky = -plane_sign
            bY, aY, gY = bY, aY + ky * bY, (1.0 + (aY + ky * bY) ** 2) / bY
            Dp = Dp - plane_sign * D
            # sector slices
            dth = theta / n_slice
            ds = rho * dth
            Ms = msector(rho, dth)
            Md = mdrift(ds)
            c = math.cos(dth)
            s = math.sin(dth)
            for _s in range(n_slice):
                oint_gx += g * ds
                oint_gy += gY * ds
                oint_bx_over_rho2 += b / (rho * rho) * ds
                b, a, g = twiss_transport(Ms, b, a, g)
                bY, aY, gY = twiss_transport(Md, bY, aY, gY)
                Dn = c * D + rho * s * Dp + rho * (1.0 - c)
                Dpn = -(s / rho) * D + c * Dp + s
                D, Dp = Dn, Dpn
    # xi numeric from the linear model
    eps = 1e-6
    nup = tune_from_matrix(oneturn(eps)[0])
    num = tune_from_matrix(oneturn(-eps)[0])
    nup_y = tune_from_matrix(oneturn(eps)[1])
    num_y = tune_from_matrix(oneturn(-eps)[1])
    return {
        "nu_x": nu_x,
        "nu_y": nu_y,
        "beta_x": bx,
        "alpha_x": ax,
        "beta_y": by,
        "alpha_y": ay,
        "D0": D0,
        "Dp0": Dp0,
        "oint_gx": oint_gx,
        "oint_gy": oint_gy,
        "oint_bx_rho2": oint_bx_over_rho2,
        "sum_k_bx": sum_k_bx,
        "sum_k_by": sum_k_by,
        "sext": sext_rec,
        "xi_x": (nup - num) / (2 * eps),
        "xi_y": (nup_y - num_y) / (2 * eps),
        "C": n_cell * 2.0 * rho * theta,
    }


def sector_identity(rho, theta, k, m_f, m_d):
    """(a) on the linear sector model: same operator for D and for x_nl."""
    n_cell = 2
    opt = sector_fodo_linear(rho, theta, k, m_f, m_d, n_cell)
    Jx, Jy = 1.0e-6, 1.7e-6
    # Periodic orbit from average sextupole kicks. Kicks do not change the
    # sector driving term; they are thin Delta x' at the sextupoles.
    kicks = [-0.5 * m * (Jx * bx - Jy * by) for bx, by, D, m in opt["sext"]]
    oint_FD = 0.0
    for (bx, by, D, m), F in zip(opt["sext"], kicks):
        oint_FD += F * D

    def once(x, xp, apply_kick):
        ik = 0
        for _c in range(n_cell):
            for plane_sign, _m in ((k, m_f), (-k, m_d)):
                xp = xp - plane_sign * x
                if apply_kick:
                    xp = xp + kicks[ik]
                    ik += 1
                c = math.cos(theta)
                s = math.sin(theta)
                xn = c * x + rho * s * xp
                xpn = -(s / rho) * x + c * xp
                x, xp = xn, xpn
        return x, xp
    z = once(0.0, 0.0, True)
    e = 1e-8
    dx = once(e, 0.0, False)  # homogeneous; kicks must not be double-counted
    # Homogeneous map by turning kicks off and using finite difference of the
    # kick map is messy. Build M from quads+sectors only.
    def hom(x, xp):
        for _c in range(n_cell):
            for plane_sign in (k, -k):
                xp = xp - plane_sign * x
                c = math.cos(theta)
                s = math.sin(theta)
                xn = c * x + rho * s * xp
                xpn = -(s / rho) * x + c * xp
                x, xp = xn, xpn
        return x, xp
    z0 = hom(0.0, 0.0)
    # z0 should be 0
    dx = hom(e, 0.0)
    dp = hom(0.0, e)
    M11 = dx[0] / e
    M21 = dx[1] / e
    M12 = dp[0] / e
    M22 = dp[1] / e
    # map with kicks of 0 is hom. map(x) = M x + particular(kicks).
    part = once(0.0, 0.0, True)
    a11, a12 = 1 - M11, -M12
    a21, a22 = -M21, 1 - M22
    det = a11 * a22 - a12 * a21
    x0 = (part[0] * a22 - part[1] * a12) / det
    xp0 = (-part[0] * a21 + part[1] * a11) / det
    # Integrate x/rho along the periodic orbit. x continuous through slices.
    x, xp = x0, xp0
    ik = 0
    oint = 0.0
    n_slice = 400
    dth = theta / n_slice
    ds = rho * dth
    c = math.cos(dth)
    s = math.sin(dth)
    for _c in range(n_cell):
        for plane_sign in (k, -k):
            xp = xp - plane_sign * x
            xp = xp + kicks[ik]
            ik += 1
            for _s in range(n_slice):
                oint += x / rho * ds
                xn = c * x + rho * s * xp
                xpn = -(s / rho) * x + c * xp
                x, xp = xn, xpn
    print("--- sector linear identity (a)")
    print("  oint x_nl/rho", oint, "oint <F> D", oint_FD, "diff", oint - oint_FD)
    print("  nu", opt["nu_x"], opt["nu_y"], "xi linear", opt["xi_x"], opt["xi_y"])
    print("  oint gamma_x", opt["oint_gx"], "sum K1 beta_x", opt["sum_k_bx"],
          "oint beta/rho^2", opt["oint_bx_rho2"])
    print("  gamma_x - K1 beta - beta/rho^2",
          opt["oint_gx"] - opt["sum_k_bx"] - opt["oint_bx_rho2"])
    C = opt["C"]
    # naive xi treats 1/rho**2 as if it scaled like K1
    xi_naive_x = -(opt["sum_k_bx"] + opt["oint_bx_rho2"]) / (4.0 * math.pi)
    xi_quad_x = -opt["sum_k_bx"] / (4.0 * math.pi)
    print("  xi from d(nu)/d(delta)", opt["xi_x"], "quad-only integral", xi_quad_x,
          "naive with 1/rho^2", xi_naive_x)
    return opt


def tune_unfolded(M):
    tr = max(-2.0, min(2.0, M[0][0] + M[1][1]))
    mu = math.acos(tr / 2.0)
    if M[0][1] < 0.0:
        mu = 2.0 * math.pi - mu
    return mu / (2.0 * math.pi)


def twiss_unfolded(M):
    nu = tune_unfolded(M)
    sm = math.sin(2.0 * math.pi * nu)
    beta = M[0][1] / sm
    alpha = (M[0][0] - M[1][1]) / (2.0 * sm)
    return beta, alpha, nu


def sector_fodo_exact(rho, theta, k, Jx, Jy, n_cell=2):
    """Exact hard-edge ring. Compare direct path length with -2*pi*J*xi/C.
    No sextupoles: the only nonlinearities are geometric.
    """
    C = n_cell * 2.0 * rho * theta

    def track_ring(x, xp, y, yp, delta):
        path = 0.0
        for _c in range(n_cell):
            xp = xp - k / (1.0 + delta) * x
            yp = yp + k / (1.0 + delta) * y
            x, xp, y, yp, p = track_sector(x, xp, y, yp, delta, rho, theta)
            path += p
            xp = xp - (-k) / (1.0 + delta) * x
            yp = yp + (-k) / (1.0 + delta) * y
            x, xp, y, yp, p = track_sector(x, xp, y, yp, delta, rho, theta)
            path += p
        return x, xp, y, yp, path

    def fixed(delta):
        v = [0.0, 0.0, 0.0, 0.0]
        for _ in range(12):
            x2, xp2, y2, yp2, _ = track_ring(*v, delta)
            r = [x2 - v[0], xp2 - v[1], y2 - v[2], yp2 - v[3]]
            if sum(abs(t) for t in r) < 1e-13:
                break
            e = 1e-8
            cols = []
            for j in range(4):
                w = v[:]
                w[j] += e
                q = track_ring(*w, delta)
                rj = [q[i] - w[i] for i in range(4)]
                cols.append([(rj[i] - r[i]) / e for i in range(4)])
            A = [[cols[j][i] for j in range(4)] for i in range(4)]
            step = solve4(A, [-r[i] for i in range(4)])
            for j in range(4):
                v[j] += step[j]
        return v

    def matrix(delta, v):
        e = 1e-7
        base = track_ring(*v, delta)
        cols = []
        for j in range(4):
            w = v[:]
            w[j] += e
            q = track_ring(*w, delta)
            cols.append([(q[i] - base[i]) / e for i in range(4)])
        Mx = ((cols[0][0], cols[1][0]), (cols[0][1], cols[1][1]))
        My = ((cols[2][2], cols[3][2]), (cols[2][3], cols[3][3]))
        return Mx, My

    v0 = fixed(0.0)
    Mx, My = matrix(0.0, v0)
    bx, ax, nu_x = twiss_unfolded(Mx)
    by, ay, nu_y = twiss_unfolded(My)
    eps = 1e-4
    vp, vm = fixed(eps), fixed(-eps)
    Mxp, Myp = matrix(eps, vp)
    Mxm, Mym = matrix(-eps, vm)
    xi_x = (tune_unfolded(Mxp) - tune_unfolded(Mxm)) / (2 * eps)
    xi_y = (tune_unfolded(Myp) - tune_unfolded(Mym)) / (2 * eps)

    def lam(sx, sy):
        nphase = 12
        # closed-orbit path
        _, _, _, _, p0 = track_ring(*v0, 0.0)
        acc = 0.0
        for i in range(nphase):
            psi = 2.0 * math.pi * i / nphase
            psy = 2.0 * math.pi * ((i * 5) % nphase) / nphase
            x = v0[0] + (math.sqrt(2 * sx * bx) * math.cos(psi) if sx else 0.0)
            xp = v0[1] + (
                -math.sqrt(2 * sx / bx) * (math.sin(psi) + ax * math.cos(psi)) if sx else 0.0
            )
            y = v0[2] + (math.sqrt(2 * sy * by) * math.cos(psy) if sy else 0.0)
            yp = v0[3] + (
                -math.sqrt(2 * sy / by) * (math.sin(psy) + ay * math.cos(psy)) if sy else 0.0
            )
            _, _, _, _, p = track_ring(x, xp, y, yp, 0.0)
            acc += p - p0
        return (acc / nphase) / C

    def lin_part(sx, sy):
        # L(f)=A f + B f^2 with f=1 and f=1/4
        L1 = lam(sx, sy)
        L2 = lam(sx / 4.0, sy / 4.0)
        A = (16.0 * L2 - L1) / 3.0
        return A, L1

    print("--- exact sector FODO rho", rho, "theta", theta, "k", k, "n_cell", n_cell, "C", C)
    print("  nu", nu_x, nu_y, "beta", bx, by, "v0", v0)
    print("  xi exact", xi_x, xi_y)
    for label, sx, sy in (("Jx", Jx, 0.0), ("Jy", 0.0, Jy), ("both", Jx, Jy)):
        A, raw = lin_part(sx, sy)
        pred = -(2.0 * math.pi / C) * (sx * xi_x + sy * xi_y)
        print(" ", label, "Lambda direct", A, "raw", raw, "Lambda xi", pred,
              "rel", (A - pred) / pred if pred else A)


def solve4(A, b):
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(4):
        piv = max(range(col, 4), key=lambda r: abs(M[r][col]))
        M[col], M[piv] = M[piv], M[col]
        div = M[col][col]
        if abs(div) < 1e-18:
            raise RuntimeError("singular")
        for j in range(col, 5):
            M[col][j] /= div
        for r in range(4):
            if r == col:
                continue
            factor = M[r][col]
            for j in range(col, 5):
                M[r][j] -= factor * M[col][j]
    return [M[i][4] for i in range(4)]


def check_green_identity_constant():
    K = 0.25
    C = 40.0
    k = math.sqrt(K)
    nu = C * k / (2.0 * math.pi)
    beta = 1.0 / k
    n = 4000
    acc = 0.0
    for i in range(n):
        s = (i + 0.5) * C / n
        psi = k * s
        G = beta / (2.0 * math.sin(math.pi * nu)) * math.cos(math.pi * nu - abs(psi))
        acc += G * (C / n)
    print("--- Green constant-K", acc, "1/K", 1.0 / K, "diff", acc - 1.0 / K)


if __name__ == "__main__":
    check_green_identity_constant()
    run_thin()
    check_sector_matrix()
    ring_uniform(8, 10.0, 1.0e-6)
    ring_uniform(8, 10.0, 1.0e-4)
    sector_identity(8.0, 0.25, 0.30, 0.0, 0.0)
    sector_identity(8.0, 0.25, 0.30, 0.8, -0.6)
    sector_fodo_exact(8.0, 0.05, 0.30, 1.0e-6, 1.5e-6, n_cell=2)
    sector_fodo_exact(8.0, 0.25, 0.30, 1.0e-6, 1.5e-6, n_cell=2)
