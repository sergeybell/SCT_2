"""Fixed-number checks for chapter 3 (phase cubic term). No symbols, no fits."""
from fractions import Fraction

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import ellipk, ellipe


def frac_identities():
    assert Fraction(1, 16) - Fraction(1, 48) == Fraction(1, 24)
    assert Fraction(1, 48) - Fraction(1, 192) == Fraction(1, 64)
    assert 2 * Fraction(1, 64) == Fraction(1, 32)
    assert Fraction(1, 16) == Fraction(12, 192)
    assert Fraction(12, 192) - Fraction(1, 192) == Fraction(11, 192)
    assert Fraction(1, 64) == Fraction(3, 192)
    assert Fraction(11, 192) - Fraction(3, 192) == Fraction(8, 192)
    assert Fraction(8, 192) == Fraction(1, 24)
    # a = phi + phi^3/192 reproduces phi_hat = a - a^3/192 through phi^5
    # checked by expanding (phi + phi**3/192) - (phi + phi**3/192)**3/192
    # = phi + phi**3/192 - phi**3/192 + O(phi**5) = phi
    print("fraction identities ok")


def elliptic_ratio(phi):
    """<delta^2> / delta_m^2 for the pendulum, k = sin(phi/2)."""
    m = np.sin(phi / 2.0) ** 2
    K = ellipk(m)
    E = ellipe(m)
    return (E / K - 1.0 + m) / m


def theta_quadrature(phi):
    """Same average from the regular theta-integrand, without ellipk."""
    m = float(np.sin(phi / 2.0) ** 2)

    def integrand(theta):
        s = np.sin(theta)
        return (np.cos(theta) ** 2) / np.sqrt(1.0 - m * s * s)

    I, err = quad(integrand, 0.0, 0.5 * np.pi, epsabs=1e-14, limit=500)
    # K = int_0^{pi/2} (1 - m sin^2)^{-1/2} d theta
    def k_int(theta):
        s = np.sin(theta)
        return 1.0 / np.sqrt(1.0 - m * s * s)

    K, errk = quad(k_int, 0.0, 0.5 * np.pi, epsabs=1e-14, limit=500)
    return I / K, err, errk


def expansion_check():
    print("--- elliptic series vs (1/2)*(1 - phi^2/32) ---")
    target = -1.0 / 32.0
    for phi in (0.4, 0.2, 0.1, 0.05, 0.025):
        r = elliptic_ratio(phi)
        # r = 1/2 * (1 - phi^2/32) + O(phi^4) = 1/2 - phi^2/64 + O(phi^4)
        corr = (2.0 * r - 1.0) / phi**2
        series = 0.5 * (1.0 - phi**2 / 32.0)
        print(
            f"phi={phi:.3f}  ratio={r:.12f}  series={series:.12f}  "
            f"(2*ratio-1)/phi^2={corr:.8f}  target={target:.8f}  "
            f"diff={corr - target:.3e}"
        )


def quadrature_vs_closed():
    print("--- theta quadrature vs closed elliptic ratio ---")
    for phi in (0.3, 1.0, 2.0, 3.0):
        q, err, errk = theta_quadrature(phi)
        c = float(elliptic_ratio(phi))
        print(
            f"phi={phi:.1f}  quad={q:.14f}  closed={c:.14f}  "
            f"diff={q - c:.3e}  quad_err={err:.1e}"
        )


def _one_period(phi_of_t_rhs, y0, t_max):
    def event_back(t, y):
        return y[0]

    event_back.direction = 1  # phi crosses 0 upward
    event_back.terminal = False
    sol = solve_ivp(
        phi_of_t_rhs,
        (0.0, t_max),
        y0,
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
        events=event_back,
        method="DOP853",
    )
    times = sol.t_events[0]
    # t=0 is the initial upward crossing; the next one is the period
    period_candidates = [t for t in times if t > 1e-8]
    if not period_candidates:
        raise RuntimeError("period not found")
    T = float(period_candidates[0])
    ts = np.linspace(0.0, T, 8001)
    phi, phidot = sol.sol(ts)
    return T, ts, phi, phidot


def fourier_sin(ts, phi, n, T):
    """Coefficient of sin(2 pi n t / T)."""
    w = 2.0 * np.pi * n / T
    # trapz of phi * sin * 2/T
    return (2.0 / T) * np.trapezoid(phi * np.sin(w * ts), ts)


def pendulum_ode_check():
    print("--- full pendulum ODE vs energy and elliptic average ---")
    # phi'' + sin phi = 0, omega_s = 1, delta = phi_dot
    def rhs(_t, y):
        return [y[1], -np.sin(y[0])]

    for phi_hat in (0.2, 1.2, 2.5):
        v0 = 2.0 * np.sin(phi_hat / 2.0)
        T, ts, phi, phidot = _one_period(rhs, [0.0, v0], t_max=50.0)
        delta = phidot  # omega_rf = |eta0| = 1
        dmax = float(np.max(delta))
        dmin = float(np.min(delta))
        dm = 0.5 * (dmax - dmin)
        d2 = float(np.trapezoid(delta**2, ts) / T)
        pred_dm = 2.0 * np.sin(phi_hat / 2.0)
        pred_d2 = dm**2 * float(elliptic_ratio(phi_hat))
        # also the series relation and series average, for comparison
        series_dm = phi_hat * (1.0 - phi_hat**2 / 24.0)
        series_d2 = (dm**2 / 2.0) * (1.0 - phi_hat**2 / 32.0)
        m = np.sin(phi_hat / 2.0) ** 2
        Omega_exact = np.pi / (2.0 * ellipk(m))
        Omega_num = 2.0 * np.pi / T
        print(
            f"phi={phi_hat:.1f}  dm={dm:.10f}  energy={pred_dm:.10f}  "
            f"diff_dm={dm - pred_dm:.3e}"
        )
        print(
            f"         <d^2>={d2:.10f}  elliptic={pred_d2:.10f}  "
            f"diff={d2 - pred_d2:.3e}  series={series_d2:.10f}"
        )
        print(
            f"         Omega={Omega_num:.10f}  pi/(2K)={Omega_exact:.10f}  "
            f"series_dm={series_dm:.10f}"
        )


def cubic_ode_check():
    print("--- cubic truncation ODE vs BM ---")
    # phi'' + (phi - phi**3/6) = 0
    def rhs(_t, y):
        p = y[0]
        return [y[1], -(p - p**3 / 6.0)]

    for phi_hat in (0.15, 0.3):
        # energy: v_max = phi_hat * sqrt(1 - phi_hat**2/12)
        v0 = phi_hat * np.sqrt(1.0 - phi_hat**2 / 12.0)
        T, ts, phi, phidot = _one_period(rhs, [0.0, v0], t_max=40.0)
        delta = phidot
        dmax = float(np.max(delta))
        dmin = float(np.min(delta))
        dm = 0.5 * (dmax - dmin)
        d2 = float(np.trapezoid(delta**2, ts) / T)
        mean_d = float(np.trapezoid(delta, ts) / T)
        pmax = float(np.max(phi))
        pmin = float(np.min(phi))
        ph = 0.5 * (pmax - pmin)
        s1 = fourier_sin(ts, phi, 1, T)
        s3 = fourier_sin(ts, phi, 3, T)
        # cosine check: delta should be c1 cos + c3 cos 3
        w = 2.0 * np.pi / T
        c1 = (2.0 / T) * np.trapezoid(delta * np.cos(w * ts), ts)
        c3 = (2.0 / T) * np.trapezoid(delta * np.cos(3.0 * w * ts), ts)
        pred_ph_coeff_s1 = phi_hat + phi_hat**3 / 192.0
        pred_ph_coeff_s3 = phi_hat**3 / 192.0
        pred_dm = phi_hat * (1.0 - phi_hat**2 / 24.0)
        # shape at fixed numerical half-span
        pred_c1 = dm * (1.0 - phi_hat**2 / 64.0)
        pred_c3 = dm * (phi_hat**2 / 64.0)
        pred_d2 = (dm**2 / 2.0) * (1.0 - phi_hat**2 / 32.0)
        Omega = w
        pred_Omega = 1.0 - phi_hat**2 / 16.0
        # also in terms of first-harmonic coefficient a of the cosine form.
        # sine coeff s1 = a_sine; for the shifted solution s1 = a_cos * (something).
        # Our formula: s1 = phi + phi**3/192, and a_cosine = phi + phi**3/192
        # because the shift preserves the first-harmonic amplitude.
        a = s1  # first harmonic amplitude after the convention phi ~ a sin
        # BM says a_method where phi_hat = a_method * (1 - a_method**2/192)
        # and s1 should equal a_method (the cosine coefficient equals the sine coefficient)
        corr = (2.0 * d2 / dm**2 - 1.0) / phi_hat**2
        print(f"phi_target={phi_hat:.2f}  phi_halfspan={ph:.10f}  mean_delta={mean_d:.3e}")
        print(
            f"  sin1={s1:.10f} pred={pred_ph_coeff_s1:.10f} diff={s1 - pred_ph_coeff_s1:.3e}"
        )
        print(
            f"  sin3={s3:.10f} pred={pred_ph_coeff_s3:.10f} diff={s3 - pred_ph_coeff_s3:.3e}"
        )
        print(f"  dm={dm:.10f} pred={pred_dm:.10f} diff={dm - pred_dm:.3e}")
        print(f"  cos1={c1:.10f} pred={pred_c1:.10f} diff={c1 - pred_c1:.3e}")
        print(f"  cos3={c3:.10f} pred={pred_c3:.10f} diff={c3 - pred_c3:.3e}")
        print(
            f"  <d^2>={d2:.10f} pred={pred_d2:.10f} diff={d2 - pred_d2:.3e}  "
            f"(2<d^2>/dm^2-1)/phi^2={corr:.8f} target={-1/32:.8f}"
        )
        print(
            f"  Omega={Omega:.10f} pred={pred_Omega:.10f} diff={Omega - pred_Omega:.3e}"
        )
        # if delta_m were replaced by the first-harmonic coefficient, the
        # relative correction to <d^2> would start at phi^4
        d2_over_c1 = d2 / (0.5 * c1**2) - 1.0
        print(f"  <d^2>/(c1^2/2) - 1 = {d2_over_c1:.3e}  (should be O(phi^4))")


def direct_phi_quadrature():
    """Average by integrating dphi/phidot, singularity subtracted via theta map
    is already above. Here: integrate in phi with a sqrt-endpoint substitution
    that does not use the elliptic modulus as the integration variable name
    in the closed form. Independent of ellipk."""
    print("--- phi-integral quadrature vs closed form ---")

    def ratio(phi_hat):
        # u from 0 to 1, phi = phi_hat - s^2, s = sqrt(phi_hat - phi)
        # near the end phidot ~ const * s
        k2 = np.sin(phi_hat / 2.0) ** 2

        def speed(phi):
            return np.sqrt(np.maximum(k2 - np.sin(phi / 2.0) ** 2, 0.0))

        # split: 0 to phi_hat - delta is regular; the last piece uses s
        # Use s over the WHOLE range: phi = phi_hat * (1 - s^2)? better Chebyshev.
        # phi = phi_hat * sin(alpha), alpha in 0..pi/2, phi=0 at alpha=0? 
        # That vanishes speed at alpha=pi/2 only.
        # phi = phi_hat * sin(alpha), dphi = phi_hat cos(alpha) d alpha
        # speed(phi_hat)=0, speed finite at 0.
        def num(alpha):
            phi = phi_hat * np.sin(alpha)
            sp = speed(phi)
            # int speed dphi = int speed * phi_hat * cos(alpha) d alpha
            return sp * phi_hat * np.cos(alpha)

        def den(alpha):
            phi = phi_hat * np.sin(alpha)
            sp = speed(phi)
            # int dphi/speed
            return phi_hat * np.cos(alpha) / sp

        # At alpha=pi/2, cos=0 and speed=0. The ratio cos/speed stays finite.
        # Avoid the exact endpoint.
        A, e1 = quad(num, 0.0, 0.5 * np.pi - 1e-12, epsabs=1e-13, limit=500)
        B, e2 = quad(den, 1e-14, 0.5 * np.pi - 1e-10, epsabs=1e-12, limit=800)
        # <phi_dot^2>/phi_dot_max^2 = (int speed dphi) / (k^2 * int dphi/speed)
        return A / (k2 * B), e1, e2

    for phi in (0.5, 1.2, 2.2):
        r, e1, e2 = ratio(phi)
        c = float(elliptic_ratio(phi))
        print(f"phi={phi:.1f}  phi_quad={r:.12f}  closed={c:.12f}  diff={r - c:.3e}")


if __name__ == "__main__":
    frac_identities()
    expansion_check()
    quadrature_vs_closed()
    direct_phi_quadrature()
    pendulum_ode_check()
    cubic_ode_check()
