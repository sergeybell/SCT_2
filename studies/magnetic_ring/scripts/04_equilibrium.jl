# Step 4. Equilibrium momentum of the closed magnetic ring, sextupoles off.
#
#   julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/04_equilibrium.jl
#
# Chapters 1–3 of notes/equilibrium-energy-v3. The cavity phase stored in
# BELL is π. The notes' phase is late-positive, so φ^N = −φ^B, and below
# transition the notes' stationary phase is 0.
#
# Columns, for each case: the note formula, the normal-form torus average,
# iteration of the degree-3 map of ring_lie_factor, and exact_track_element
# with rf_kick on the thin cavity. Nothing is imported from
# Julia/equilibrium_energy. --turns N replaces the default of three
# synchrotron periods. That flag is the explicit longer or shorter run.

import Pkg

const STUDY_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(STUDY_DIR)

using LinearAlgebra
using Printf
using Bell
using Bell.NormalForm: dragt_finn_nf6d
using MagneticRing

include(joinpath(STUDY_DIR, "src", "generator_chromaticity.jl"))
include(joinpath(STUDY_DIR, "src", "readings.jl"))
using .GeneratorChromaticity
using .LieReadings

const OUT_DIR = joinpath(STUDY_DIR, "output")
const REGIME_JSON = joinpath(OUT_DIR, "02_rf_regimes.json")
const TRACK_DS = 0.1
const PERIODS = 3
const BELL_PHI_S = Float64(π)
const XI_FACTOR = -1 / (2π)
# Bulk δ_eq. The next omitted slip term is O(δ_m) relative to η1 when η2 ~ η1
# (δ_m = 1e-4). Half a percent is past that scale and past φ̂²/32.
const BULK_REL = 0.005
# Chapter 3 drops an O(φ̂²) piece of the second harmonic. Linear φ̂² is 0.0025
# and cubic φ̂² is 0.09. Losing the velocity piece −2/γ² moves the coefficient
# by about 40%, which is outside both bands.
const A2_TOL_LINEAR = 0.05
const A2_TOL_CUBIC = 0.20

const STOPS = String[]

function stop!(msg)
    push!(STOPS, msg)
    println("STOP: ", msg)
    flush(stdout)
    return nothing
end

struct SeriesMeas
    kind::String
    delta_eq::Float64
    delta_m::Float64
    a1::Float64
    a2::Float64
    a2_cos::Float64
    a2_sin::Float64
    a3::Float64
    Q::Float64
    delta2::Float64
    phi_hat::Float64
    rms::Float64
    unc_delta_eq::Float64
    unc_Q::Float64
    unc_a2::Float64
    gamma_eff::Float64
    dE_over_E::Float64
    n::Int
    aborted::Bool
end

struct Targets
    delta_eq_ch12::Float64
    delta_eq_ch3::Float64
    delta_eq_pendulum::Float64
    a2::Float64
    omega_ch3::Float64
    omega_pendulum::Float64
    delta2_ch12::Float64
    delta2_ch3::Float64
    delta2_pendulum::Float64
    gamma_formula::Float64
    dE_formula::Float64
    Lambda::Float64
    phi_hat::Float64
    cubic::Bool
end

# --- note formulas -----------------------------------------------------------

"Complete elliptic integrals K(k) and E(k) on a uniform θ grid."
function elliptic_KE(k2::Float64; n::Int = 8000)
    k2 < 1 || error("elliptic modulus must be below 1")
    sK = 0.0
    sE = 0.0
    for i in 1:n
        θ = (i - 0.5) * π / (2n)
        w = 1 - k2 * sin(θ)^2
        sK += 1 / sqrt(w)
        sE += sqrt(w)
    end
    dθ = π / (2n)
    return sK * dθ, sE * dθ
end

"""⟨δ²⟩ / (δ_m²/2) for the full pendulum. Equals 1 − φ̂²/32 through O(φ̂²)."""
function pendulum_factor(phi_hat::Float64)
    phi_hat > 0 || return 1.0
    k = sin(phi_hat / 2)
    K, E = elliptic_KE(k^2)
    ratio = (E / K - 1 + k^2) / k^2
    return 2 * ratio
end

"Ω/ω_s of the full pendulum, π / (2 K(sin(φ̂/2)))."
function pendulum_omega_ratio(phi_hat::Float64)
    phi_hat > 0 || return 1.0
    k = sin(phi_hat / 2)
    K, _ = elliptic_KE(k^2)
    return π / (2K)
end

function slip_mean_square(delta_m, factor)
    return delta_m^2 / 2 * factor
end

function delta_eq_of(eta0, eta1, delta2, Lambda)
    return -(eta1 * delta2 + Lambda) / eta0
end

"Coefficient of cos 2ψ in chapter 1. The cubic chapter does not correct it."
function second_harmonic(eta0, eta1, gamma, delta_m)
    return (delta_m^2 / 6) * (eta1 / eta0 - 2 / gamma^2)
end

function kinematics(gamma, beta, delta_eq, delta2)
    b2 = beta^2
    g_eff = gamma + b2 * gamma * delta_eq + b2 * delta2 / (2 * gamma)
    dE = b2 * (delta_eq + delta2 / (2 * gamma^2))
    return g_eff, dE
end

function make_targets(eta0, eta1, gamma, beta, delta_m, phi_hat, Lambda, cubic)
    f_ch3 = 1 - phi_hat^2 / 32
    f_pend = pendulum_factor(phi_hat)
    d2_12 = slip_mean_square(delta_m, 1.0)
    d2_3 = slip_mean_square(delta_m, f_ch3)
    d2_p = slip_mean_square(delta_m, f_pend)
    de12 = delta_eq_of(eta0, eta1, d2_12, Lambda)
    de3 = delta_eq_of(eta0, eta1, d2_3, Lambda)
    dep = delta_eq_of(eta0, eta1, d2_p, Lambda)
    d2_f = cubic ? d2_3 : d2_12
    de_f = cubic ? de3 : de12
    ge, dE = kinematics(gamma, beta, de_f, d2_f)
    return Targets(de12, de3, dep, second_harmonic(eta0, eta1, gamma, delta_m),
                   1 - phi_hat^2 / 16, pendulum_omega_ratio(phi_hat),
                   d2_12, d2_3, d2_p, ge, dE, Lambda, phi_hat, cubic)
end

# --- harmonic least squares (rewritten; not the old study) ------------------

function _design(times, Q)
    θ = (2π * Q) .* times
    return hcat(ones(length(times)), cos.(θ), sin.(θ), cos.(2 .* θ), sin.(2 .* θ),
                cos.(3 .* θ), sin.(3 .* θ))
end

function _coef_at(deltas, times, Q)
    A = _design(times, Q)
    return A \ deltas
end

function _rms_at(deltas, times, Q)
    coef = _coef_at(deltas, times, Q)
    resid = deltas .- _design(times, Q) * coef
    return coef, sqrt(sum(abs2, resid) / length(deltas))
end

"Minimise the harmonic residual in a ±2% window around Q0, then zoom."
function _refine_Q(deltas, Q0)
    times = float.(0:length(deltas)-1)
    center = Q0
    half = 0.02 * Q0
    bestQ = Q0
    best = Inf
    for _ in 1:6
        step = 2half / 40
        for i in 0:40
            Q = center - half + i * step
            Q > 0 || continue
            _, rms = _rms_at(deltas, times, Q)
            if rms < best
                best = rms
                bestQ = Q
            end
        end
        center = bestQ
        half = step
    end
    return bestQ
end

function _signed_second(coef)
    a1 = hypot(coef[2], coef[3])
    a1 > 0 || return (a2 = 0.0, a2_cos = 0.0, a2_sin = 0.0, a1 = 0.0, a3 = hypot(coef[6], coef[7]))
    cφ = coef[2] / a1
    sφ = coef[3] / a1
    c2 = cφ^2 - sφ^2
    s2 = 2 * sφ * cφ
    C, D = coef[4], coef[5]
    return (
        a2 = hypot(C, D),
        a2_cos = C * c2 + D * s2,
        a2_sin = -C * s2 + D * c2,
        a1 = a1,
        a3 = hypot(coef[6], coef[7]),
    )
end

"""
    harmonic_synchrotron(deltas, phis, Q0; Qfixed=nothing) -> named tuple

δ_n = δ_eq + a1 cos ψ + a2 cos 2ψ + a3 cos 3ψ, with ψ measured from the
crest of the fundamental. Q is refined on a grid unless `Qfixed` is set.
The split of the record is the uncertainty of δ_eq, Q, and the second harmonic.
"""
function harmonic_synchrotron(deltas, phis, Q0; Qfixed = nothing)
    n = length(deltas)
    n >= 16 || error("harmonic fit needs at least 16 samples, got $n")
    times = float.(0:n-1)
    Q = Qfixed === nothing ? _refine_Q(deltas, Q0) : Float64(Qfixed)
    coef, rms = _rms_at(deltas, times, Q)
    harm = _signed_second(coef)
    deq = coef[1]
    # Period average of the fitted waveform. A few synchrotron periods at the
    # shifted tune are not an integer window, so the raw second moment leaks.
    d2 = (harm.a1^2 + harm.a2^2 + harm.a3^2) / 2
    half = (maximum(deltas) - minimum(deltas)) / 2
    ph = length(phis) == n ? (maximum(phis) - minimum(phis)) / 2 : NaN
    n2 = n ÷ 2
    # Global time on each half, so the constant is not a leaked fundamental.
    c1 = _coef_at(deltas[1:n2], times[1:n2], Q)
    c2 = _coef_at(deltas[n2+1:end], times[n2+1:end], Q)
    h1 = _signed_second(c1)
    h2 = _signed_second(c2)
    unc_deq = abs(c1[1] - c2[1]) / 2
    unc_a2 = abs(h1.a2 - h2.a2) / 2
    unc_Q = 0.0
    if Qfixed === nothing && n2 >= 16
        unc_Q = abs(_refine_Q(deltas[1:n2], Q0) - _refine_Q(deltas[n2+1:end], Q0)) / 2
    end
    return (
        delta_eq = deq,
        Q = Q,
        a1 = harm.a1,
        a2 = harm.a2,
        a2_cos = harm.a2_cos,
        a2_sin = harm.a2_sin,
        a3 = harm.a3,
        delta_m_half = half,
        delta2 = d2,
        phi_hat = ph,
        rms = rms,
        unc_delta_eq = unc_deq,
        unc_Q = unc_Q,
        unc_a2 = unc_a2,
        n = n,
    )
end

function _meas(kind, fit, gamma, beta; aborted = false, Q = fit.Q)
    ge, dE = kinematics(gamma, beta, fit.delta_eq, fit.delta2)
    return SeriesMeas(kind, fit.delta_eq, fit.delta_m_half, fit.a1, fit.a2, fit.a2_cos,
                      fit.a2_sin, fit.a3, Q, fit.delta2, fit.phi_hat, fit.rms,
                      fit.unc_delta_eq, fit.unc_Q, fit.unc_a2, ge, dE, fit.n, aborted)
end

function formula_meas(targets, gamma, beta, Q_note, delta_m)
    de = targets.cubic ? targets.delta_eq_ch3 : targets.delta_eq_ch12
    d2 = targets.cubic ? targets.delta2_ch3 : targets.delta2_ch12
    # Chapter 1 leaves the frequency at ω_s. Chapter 3 multiplies by 1−φ̂²/16.
    ratio = targets.cubic ? targets.omega_ch3 : 1.0
    ge, dE = kinematics(gamma, beta, de, d2)
    return SeriesMeas("formula", de, delta_m, delta_m, abs(targets.a2), targets.a2, 0.0, NaN,
                      Q_note * ratio, d2, targets.phi_hat, 0.0, 0.0, 0.0, 0.0,
                      ge, dE, 0, false)
end

# --- self check, before any lattice work ------------------------------------

function _self_check()
    Q = 0.001
    n = 4000
    deq = 1.23e-8
    dm = 1.0e-4
    a2 = -5.7e-9
    φ = 0.2
    t = float.(0:n-1)
    θ = 2π * Q .* t .- φ
    delta = deq .+ dm .* cos.(θ) .+ a2 .* cos.(2 .* θ)
    ph = 0.05 .* sin.(θ)
    fit = harmonic_synchrotron(delta, ph, 0.00102)
    abs(fit.delta_eq - deq) < 1e-12 ||
        error("harmonic fit missed delta_eq: $(fit.delta_eq)")
    abs(fit.a1 - dm) < 1e-12 || error("harmonic fit missed the fundamental")
    abs(fit.a2_cos - a2) < 1e-12 ||
        error("signed second harmonic is $(fit.a2_cos), expected $a2")
    abs(fit.a2_sin) < 1e-12 || error("second-harmonic sine leaked: $(fit.a2_sin)")
    abs(fit.Q - Q) < 1e-8 || error("refined Q is $(fit.Q)")
    abs(fit.delta2 - (dm^2 + a2^2) / 2) < 1e-18 ||
        error("period-average mean square is $(fit.delta2)")
    eta0, eta1, dm0 = -0.5, 1.0, 1e-4
    de = delta_eq_of(eta0, eta1, dm0^2 / 2, 0.0)
    abs(de - 1e-8) < 1e-18 || error("chapter-1 delta_eq sign or factor")
    a = second_harmonic(eta0, eta1, 1.0, dm0)
    abs(a - (dm0^2 / 6) * (eta1 / eta0 - 2)) < 1e-18 || error("second-harmonic formula")
    Λ = betatron_lambda(100.0, 2.0, 3.0, -1.5, -0.5)
    abs(Λ - (-2π / 100) * (2 * -1.5 + 3 * -0.5)) < 1e-12 || error("Lambda formula")
    abs(pendulum_factor(0.3) - 0.997176908) < 1e-7 ||
        error("pendulum factor at phi_hat=0.3 is $(pendulum_factor(0.3))")
    abs(pendulum_omega_ratio(0.3) - 0.994377614) < 1e-7 ||
        error("pendulum tune ratio is $(pendulum_omega_ratio(0.3))")
    # ⟨δ²⟩ = δ_m²/2 * factor must reproduce both kinematic terms.
    b2, γ, deq2, fac = 0.25, 1.2, 2e-8, 0.99
    d2 = slip_mean_square(1e-3, fac)
    g1, e1 = kinematics(γ, sqrt(b2), deq2, d2)
    g2 = γ + b2 * γ * deq2 + b2 * (1e-3)^2 / (4γ) * fac
    e2 = b2 * (deq2 + (1e-3)^2 / (4γ^2) * fac)
    abs(g1 - g2) < 1e-18 || error("gamma_eff does not match the note")
    abs(e1 - e2) < 1e-18 || error("Delta E/E does not match the note")
    println("Formula and harmonic-fit self-check passed")
    flush(stdout)
    return nothing
end

# --- JSON, without a new dependency ------------------------------------------

function _json_object(text, key)
    m = findfirst(Regex("\"$(key)\"\\s*:\\s*\\{"), text)
    m === nothing && error("JSON object \"$key\" was not found")
    start = last(m)
    depth = 0
    stop = start
    for k in start:lastindex(text)
        c = text[k]
        if c == '{'
            depth += 1
        elseif c == '}'
            depth -= 1
            if depth == 0
                stop = k
                break
            end
        end
    end
    depth == 0 || error("JSON object \"$key\" does not close")
    return text[start:stop]
end

function _json_number(block, key)
    m = match(Regex("\"$(key)\"\\s*:\\s*([-+0-9.eE]+)"), block)
    m === nothing && error("JSON object has no number \"$key\"")
    return parse(Float64, m.captures[1])
end

function _load_inputs()
    isfile(REGIME_JSON) || error("regime file is missing: $REGIME_JSON")
    text = read(REGIME_JSON, String)
    optics = _json_object(text, "optics")
    placeholder = _json_object(text, "placeholder_action")
    function regime(name)
        block = _json_object(text, name)
        return (
            name = name,
            voltage_MV = _json_number(block, "voltage_MV"),
            harmonic = Int(round(_json_number(block, "harmonic"))),
            delta_m = _json_number(block, "delta_m"),
            phi_hat = _json_number(block, "phi_hat"),
            Qs = _json_number(block, "Qs"),
            phi_s_rad = _json_number(block, "phi_s_rad"),
        )
    end
    return (
        text = text,
        eta0_json = _json_number(optics, "eta0"),
        eta1_json = _json_number(optics, "eta1"),
        Jx = _json_number(placeholder, "Jx"),
        Jy = _json_number(placeholder, "Jy"),
        Lambda_json = _json_number(placeholder, "Lambda"),
        linear = regime("linear"),
        cubic = regime("cubic_phase"),
    )
end

function _require_bare_ring(lattice)
    ref = lattice.ref
    ref.particle.name == "deuteron" || error("reference particle is $(ref.particle.name)")
    ref.Ekin_MeV == 270.0 || error("kinetic energy is $(ref.Ekin_MeV) MeV")
    abs(total_bend_angle(lattice) - 2π) < 1.0e-9 ||
        error("total bend is not 2π: $(total_bend_angle(lattice))")
    any(e -> e isa RFCavity, lattice.elements) &&
        error("the closed lattice already contains an RF cavity")
    for e in lattice.elements
        e isa Sextupole && e.k2 != 0 &&
            error("sextupole $(e.name) has k2 = $(e.k2); this step reads k2 = 0")
    end
    return nothing
end

function _requested_turns()
    i = findfirst(==("--turns"), ARGS)
    i === nothing && return nothing
    i < length(ARGS) || error("--turns needs an integer")
    n = parse(Int, ARGS[i + 1])
    n >= 16 || error("--turns must be at least 16, got $n")
    return n
end

function _case_turns(Qs)
    requested = _requested_turns()
    requested !== nothing && return requested
    Qs > 0 || error("synchrotron tune is not positive")
    return max(80, round(Int, PERIODS / Qs))
end

# --- tracking ----------------------------------------------------------------

function _phi_note(z, ref, C, h)
    return -canonical_to_laboratory(z, ref, C, h)[5]
end

function _launch(po, delta, Jx, Jy)
    ax = Jx > 0 ? sqrt(2 * Jx * po.x.beta) : 0.0
    ay = Jy > 0 ? sqrt(2 * Jy * po.y.beta) : 0.0
    Dx, Dpx, Dy, Dpy = po.D4
    x = ax + Dx * delta
    px = -po.x.alpha / po.x.beta * ax + Dpx * delta
    y = ay + Dy * delta
    py = -po.y.alpha / po.y.beta * ay + Dpy * delta
    z = [x, px, y, py, 0.0, delta]
    return z, ax, ay
end

function _run_series(step, z0, nturns, ref, C, h; abort_at, label)
    deltas = Vector{Float64}(undef, nturns + 1)
    phis = Vector{Float64}(undef, nturns + 1)
    z = collect(Float64, z0)
    t0 = time()
    used = 0
    aborted = false
    for i in 0:nturns
        if i > 0
            z = step(z)
        end
        used = i + 1
        if !all(isfinite, z)
            deltas[used] = NaN
            phis[used] = NaN
            aborted = true
            println("  ", label, " non-finite state at turn ", i)
            flush(stdout)
            break
        end
        deltas[used] = z[6]
        phis[used] = _phi_note(z, ref, C, h)
        if abs(z[6]) > abort_at
            aborted = true
            @printf("  %s |delta| %.3e exceeds %.3e at turn %d\n", label, abs(z[6]), abort_at, i)
            flush(stdout)
            break
        end
        if i == 0 || i == nturns || i % 500 == 0
            @printf("  %s turn %6d / %d   elapsed %7.1f s   delta % .6e\n",
                    label, i, nturns, time() - t0, z[6])
            flush(stdout)
        end
    end
    return (
        delta = deltas[1:used],
        phi_note = phis[1:used],
        aborted = aborted,
        elapsed = time() - t0,
    )
end

function _map_step(M)
    return z -> Float64[M[i](z) for i in 1:6]
end

function _exact_step(lattice, ref, C)
    return z -> track_turn_exact(z, lattice, ref, C; ds = TRACK_DS)
end

function _fixed_point_delta(M)
    z = zeros(6)
    return Float64[M[i](z) for i in 1:6]
end

# --- normal-form torus -------------------------------------------------------

function _torus(nf, J; nth = 24)
    r = sqrt.(2 .* Float64.(collect(J)))
    ni = r[1] == 0 ? 1 : nth
    nj = r[2] == 0 ? 1 : nth
    series = zeros(nth)
    counts = zeros(Int, nth)
    acc = 0.0
    n = 0
    for i in 0:ni-1, j in 0:nj-1, k in 0:nth-1
        θx = 2π * i / ni
        θy = 2π * j / nj
        θs = 2π * k / nth
        w = [r[1] * cos(θx), -r[1] * sin(θx), r[2] * cos(θy), -r[2] * sin(θy),
             r[3] * cos(θs), -r[3] * sin(θs)]
        d = Float64(nf.U6(w)) + nf.zstar[6]
        acc += d
        n += 1
        series[k + 1] += d
        counts[k + 1] += 1
    end
    series ./= counts
    mean = acc / n
    # The angle grid is one exact synchrotron period. Do not refit Q.
    fit = harmonic_synchrotron(series, zeros(nth), 1 / nth; Qfixed = 1 / nth)
    actions = collect(Float64, J)
    d_dot = nf.zstar[6] + dot(nf.d, actions)
    g_dot = nf.zstar[6] + dot(nf.d_grid, actions)
    return (mean = mean, fit = fit, d_dot = d_dot, grid_dot = g_dot, samples = n)
end

function _js_for_delta_m(nf, delta_m)
    σ = nf.sigma_delta_s
    σ > 0 || error("sigma_delta_s is not positive: $σ")
    return delta_m^2 / (2 * σ^2), σ
end

# --- one regime: factor, normal form, map -----------------------------------

function _build_cavity_ring(lattice, voltage, harmonic)
    lat, cav = append_thin_cavity(lattice, voltage, harmonic, BELL_PHI_S)
    lat.harmonic == cav.harmonic ||
        error("ring harmonic $(lat.harmonic) differs from the cavity harmonic $(cav.harmonic)")
    println("  ring_lie_factor of the ring with one thin cavity")
    flush(stdout)
    F = ring_lie_factor(lat)
    println("  degree-3 polynomial, factor_to_map order 3")
    flush(stdout)
    M = factor_to_map(F; order = 3)
    println("  dragt_finn_nf6d")
    flush(stdout)
    nf = dragt_finn_nf6d(F)
    return (lattice = lat, cav = cav, F = F, M = M, nf = nf)
end

# --- judgment ----------------------------------------------------------------

function _judge_bulk!(tag, meas, formula)
    scale = max(abs(formula), 1e-30)
    gap = abs(meas - formula)
    rel = gap / scale
    if rel > BULK_REL
        stop!("$tag bulk delta_eq relative gap $(rel) exceeds $(BULK_REL). measured $(meas), formula $(formula)")
    else
        @printf("  agree  %-22s  rel %.3e\n", tag, rel)
        flush(stdout)
    end
    return rel
end

function _judge_cubic_factor!(tag, meas, unc, targets)
    correction = abs(targets.delta_eq_ch3 - targets.delta_eq_ch12)
    series = abs(targets.delta_eq_pendulum - targets.delta_eq_ch3)
    gap3 = abs(meas - targets.delta_eq_ch3)
    gap12 = abs(meas - targets.delta_eq_ch12)
    if unc >= correction / 4
        @printf("  %-22s  cubic factor not resolved (unc %.3e, correction %.3e)\n",
                tag, unc, correction)
        flush(stdout)
        if min(gap3, gap12) > max(BULK_REL * abs(targets.delta_eq_ch12), 8 * unc)
            stop!("$tag is far from both the chapter 1-2 and chapter 3 values. measured $meas")
        else
            stop!("$tag did not resolve 1 - phi_hat^2/32. unc $unc, correction $correction, measured $meas")
        end
        return
    end
    if gap12 + correction / 4 < gap3
        stop!(string(
            "$tag δ_eq is closer to chapter 1-2 than to 1 - phi_hat^2/32. ",
            "measured $meas, chapter 1-2 $(targets.delta_eq_ch12), chapter 3 $(targets.delta_eq_ch3). ",
            "The tune and the third harmonic of this series do carry phi_hat^2; ",
            "the mean momentum does not.",
        ))
        return
    end
    band = max(8 * unc, 5 * series, 0.15 * correction)
    if gap3 > band
        stop!("$tag misses the cubic delta_eq. measured $meas, formula $(targets.delta_eq_ch3), gap $gap3, band $band")
    else
        @printf("  agree  %-22s  cubic factor, gap %.3e, band %.3e\n", tag, gap3, band)
        flush(stdout)
    end
    return nothing
end

function _judge_a2!(tag, meas::SeriesMeas, targets, tol)
    formula = targets.a2
    scale = max(abs(formula), 1e-30)
    if meas.rms > abs(formula) && meas.unc_a2 > abs(formula)
        stop!("$tag second harmonic is below the residual. rms $(meas.rms), formula $formula")
        return
    end
    gap = abs(meas.a2_cos - formula)
    if gap > max(tol * scale, 5 * meas.unc_a2)
        stop!("$tag second-harmonic coefficient $(meas.a2_cos) vs formula $formula (tol $tol)")
    elseif abs(meas.a2_sin) > max(tol * scale, 5 * meas.unc_a2)
        stop!("$tag second harmonic has a quadrature piece $(meas.a2_sin) beside coefficient $(meas.a2_cos)")
    else
        @printf("  agree  %-22s  second harmonic rel %.3e\n", tag, gap / scale)
        flush(stdout)
    end
    return nothing
end

function _judge_tune!(tag, Q, unc_Q, Q_note, targets)
    ratio = Q / Q_note
    unc = unc_Q / Q_note
    if !targets.cubic
        if abs(ratio - 1) > 0.01
            stop!("$tag linear Q/Q_note = $ratio, off by more than 1%")
        else
            @printf("  agree  %-22s  Q/Q_note %.8f\n", tag, ratio)
            flush(stdout)
        end
        return
    end
    predicted = targets.omega_ch3
    series = abs(targets.omega_pendulum - predicted)
    if unc >= abs(predicted - 1) / 4
        stop!("$tag did not resolve Omega/omega_s. ratio $ratio, unc $unc, predicted $predicted")
        return
    end
    if abs(ratio - predicted) > max(0.0015, 5 * unc, 5 * series)
        stop!("$tag Omega/omega_s = $ratio, chapter 3 predicts $predicted")
    else
        @printf("  agree  %-22s  Omega/omega_s %.8f vs %.8f\n", tag, ratio, predicted)
        flush(stdout)
    end
    return nothing
end

function _judge_nf!(tag, torus, targets)
    formula = targets.delta_eq_ch12
    scale = max(abs(formula), 1e-30)
    gap = abs(torus.mean - formula)
    spread = maximum(abs, (torus.mean - torus.d_dot, torus.mean - torus.grid_dot))
    if spread > max(1e-3 * max(abs(torus.d_dot), 1e-30), 1e-14)
        stop!("$tag torus mean $(torus.mean), d·J $(torus.d_dot), d_grid·J $(torus.grid_dot)")
    end
    # The normal form stops at the quadratic action term, so φ̂² is absent.
    # 0.1% is tighter than φ̂²/32 at φ̂ = 0.3 (0.28%) and looser than the
    # generator-versus-orders chromaticity gap.
    if gap / scale > 1e-3
        stop!("$tag torus $(torus.mean) vs chapter 1-2 $formula, relative $(gap / scale)")
    else
        @printf("  agree  %-22s  NF vs ch.1-2 rel %.3e\n", tag, gap / scale)
        flush(stdout)
    end
    return nothing
end

function _judge_phihat!(tag, measured, predicted)
    predicted > 0 || return nothing
    rel = abs(measured - predicted) / predicted
    # The linear φ̂ and the half-range differ by about φ̂²/24.
    if rel > 0.02
        stop!("$tag phi_hat $(measured) vs regime $(predicted), relative $rel")
    else
        @printf("  agree  %-22s  phi_hat rel %.3e\n", tag, rel)
        flush(stdout)
    end
    return nothing
end

function _judge_identity!(tag, meas, eta0, eta1, Lambda)
    # ⟨δ²⟩ in the note is the mean of δ². The fit stores the variance about δ_eq;
    # δ_eq² is negligible next to that variance at these amplitudes.
    predicted = delta_eq_of(eta0, eta1, meas.delta2 + meas.delta_eq^2, Lambda)
    _judge_bulk!(tag * " identity", meas.delta_eq, predicted)
    return predicted
end

function _judge_case!(name, targets, map_meas, exact_meas, Q_note, torus, eta0, eta1)
    println("Judgement ", name)
    flush(stdout)
    a2tol = targets.cubic ? A2_TOL_CUBIC : A2_TOL_LINEAR
    _judge_nf!(name * " NF", torus, targets)
    for (tag, meas) in ((name * " map", map_meas), (name * " exact", exact_meas))
        meas.aborted && stop!("$tag left the momentum window")
        meas.aborted && continue
        if targets.cubic
            _judge_cubic_factor!(tag, meas.delta_eq, meas.unc_delta_eq, targets)
        else
            _judge_bulk!(tag, meas.delta_eq, targets.delta_eq_ch12)
        end
        _judge_identity!(tag, meas, eta0, eta1, targets.Lambda)
        _judge_a2!(tag, meas, targets, a2tol)
        _judge_tune!(tag, meas.Q, meas.unc_Q, Q_note, targets)
        _judge_phihat!(tag, meas.phi_hat, targets.phi_hat)
    end
    return nothing
end

# --- phase dictionary --------------------------------------------------------

function _phase_probe(lattice, ref, C, po, regime, phi_s, Qs)
    println("Phase probe phi_s = ", phi_s)
    flush(stdout)
    lat, _ = append_thin_cavity(lattice, regime.voltage_MV, regime.harmonic, phi_s)
    z0, _, _ = _launch(po, regime.delta_m, 0.0, 0.0)
    nturns = max(50, round(Int, 1 / Qs))
    local series
    try
        series = _run_series(_exact_step(lat, ref, C), z0, nturns, ref, C, regime.harmonic;
                             abort_at = 0.01, label = "phi_s=$(round(phi_s; digits=3))")
    catch err
        println("Phase probe phi_s = ", phi_s, " failed: ", sprint(showerror, err))
        flush(stdout)
        return (
            phi_s = phi_s, n = 0, nturns = nturns, aborted = true,
            minimum = NaN, maximum = NaN, max_abs = NaN, final = NaN, elapsed = 0.0,
            oscillates = false, ran_away = true,
        )
    end
    d = series.delta
    finite = filter(isfinite, d)
    mx = maximum(abs, finite)
    return (
        phi_s = phi_s,
        n = length(d),
        nturns = nturns,
        aborted = series.aborted,
        minimum = minimum(finite),
        maximum = maximum(finite),
        max_abs = mx,
        final = finite[end],
        elapsed = series.elapsed,
        oscillates = minimum(finite) < -0.3 * regime.delta_m && maximum(finite) > 0.3 * regime.delta_m,
        ran_away = series.aborted || mx > 10 * regime.delta_m,
    )
end

# --- report ------------------------------------------------------------------

function _fmt(x)
    return isfinite(x) ? @sprintf("%.8e", float(x)) : "nan"
end

function _fmt_ratio(x)
    return isfinite(x) ? @sprintf("%.8f", float(x)) : "nan"
end

function _csv_escape(s)
    t = string(s)
    if occursin(',', t) || occursin('"', t)
        return "\"" * replace(t, '"' => "\"\"") * "\""
    end
    return t
end

function _push!(rows, case, column, quantity, value)
    push!(rows, (case, column, quantity, value))
    return nothing
end

function _push_meas!(rows, case, column, m::SeriesMeas)
    _push!(rows, case, column, "delta_eq", m.delta_eq)
    _push!(rows, case, column, "delta_m_half", m.delta_m)
    _push!(rows, case, column, "a1", m.a1)
    _push!(rows, case, column, "a2", m.a2)
    _push!(rows, case, column, "a2_cos", m.a2_cos)
    _push!(rows, case, column, "a2_sin", m.a2_sin)
    _push!(rows, case, column, "a3", m.a3)
    _push!(rows, case, column, "Q", m.Q)
    _push!(rows, case, column, "delta2", m.delta2)
    _push!(rows, case, column, "phi_hat", m.phi_hat)
    _push!(rows, case, column, "rms", m.rms)
    _push!(rows, case, column, "unc_delta_eq", m.unc_delta_eq)
    _push!(rows, case, column, "unc_Q", m.unc_Q)
    _push!(rows, case, column, "unc_a2", m.unc_a2)
    _push!(rows, case, column, "gamma_eff", m.gamma_eff)
    _push!(rows, case, column, "delta_E_over_E", m.dE_over_E)
    return nothing
end

function _md_row(cells)
    return "| " * join(cells, " | ") * " |"
end

function _meas_cells(m::SeriesMeas, Q_note)
    ratio = Q_note > 0 ? m.Q / Q_note : NaN
    return [
        _fmt(m.delta_eq),
        _fmt(m.delta_m),
        _fmt(m.a2_cos),
        _fmt(m.a2),
        _fmt_ratio(ratio),
        _fmt(m.delta2),
        _fmt(m.phi_hat),
        _fmt(m.gamma_eff),
        _fmt(m.dE_over_E),
    ]
end

function _write_numbers(path, rows)
    open(path, "w") do io
        println(io, "case,column,quantity,value")
        for (case, column, quantity, value) in rows
            println(io, case, ",", column, ",", quantity, ",", _fmt(value))
        end
    end
    return nothing
end

function _write_report(path, ctx, probes, cases, rows_header)
    lines = String[]
    push!(lines, "# Equilibrium momentum of the closed magnetic ring")
    push!(lines, "")
    if !isempty(STOPS)
        push!(lines, "## STOP")
        push!(lines, "")
        for s in STOPS
            push!(lines, "- " * s)
        end
        push!(lines, "")
    end
    push!(lines, "Sextupoles are at `k2 = 0`. One thin cavity is appended. BELL stores `phi_s = π` because `l > 0` is an early arrival. The notes' phase is late-positive, so `φ^N = −φ^B`. Below transition the notes' stationary phase is 0, and `η0 cos φ_s^B > 0` at `φ_s^B = π`.")
    push!(lines, "")
    push!(lines, "The normal-form torus is the average of `U6` on the action torus. It is the quadratic piece `z*_6 + d · J` (`d_grid` is the same average at action `1e-8`). That polynomial stops before `φ̂²`, so the torus is compared with chapters 1–2 even in the cubic-phase regime.")
    push!(lines, "")
    push!(lines, "`δ_eq`, `δ_m`, and the second harmonic of a tracked series come from a harmonic least-squares fit, harmonics 1–3, with the tune refined on a grid. `δ_m` in the table is the half-range `(max − min) / 2`. The signed second harmonic is the coefficient of `cos 2ψ` after `ψ` is measured from the crest of the fundamental. `⟨δ²⟩` is `(a1² + a2² + a3²) / 2`, the average of that waveform over a full period. `γ_eff` and `ΔE/E` use that column's `δ_eq` and `⟨δ²⟩`. Agreement is judged on `δ_eq`.")
    push!(lines, "")
    push!(lines, "## Phase dictionary")
    push!(lines, "")
    push!(lines, "Linear-regime voltage, `J = 0`, launch at `δ = δ_m` and `l = 0`. The control uses the same voltage with `phi_s = 0`. A run that reaches `|δ| = 0.01` stops early.")
    push!(lines, "")
    push!(lines, _md_row(["phi_s BELL", "turns asked", "samples", "aborted", "min delta", "max delta", "max |delta|", "oscillates", "left bucket"]))
    push!(lines, _md_row(["---", "---", "---", "---", "---", "---", "---", "---", "---"]))
    for p in probes
        push!(lines, _md_row([
            _fmt(p.phi_s), string(p.nturns), string(p.n), string(p.aborted),
            _fmt(p.minimum), _fmt(p.maximum), _fmt(p.max_abs),
            string(p.oscillates), string(p.ran_away),
        ]))
    end
    push!(lines, "")
    push!(lines, "## Slip, chromaticity, inequalities")
    push!(lines, "")
    push!(lines, _md_row(["quantity", "value"]))
    push!(lines, _md_row(["---", "---"]))
    for (q, v) in rows_header
        push!(lines, _md_row([q, _fmt(v)]))
    end
    push!(lines, "")
    push!(lines, "Generator `ξ` is `−c / (2π)` for the coefficient of `J δ`. `Λ = −2π (J_x ξ_x + J_y ξ_y) / C` uses that `ξ` and the placeholder actions of `02_rf_regimes.json`.")
    push!(lines, "")
    push!(lines, "`|η1 δ_m / η0|`, `δ_m`, `|Λ / (η0 δ_m)|`, the linear `φ̂`, the cubic `φ̂²`, and `ω_s T_s` (also divided by the betatron phase) are much less than 1. The cubic `φ̂` is not required to be much less than 1.")
    push!(lines, "")
    for case in cases
        t = case.targets
        push!(lines, "## " * case.title)
        push!(lines, "")
        push!(lines, case.blurb)
        push!(lines, "")
        push!(lines, _md_row(["target", "value"]))
        push!(lines, _md_row(["---", "---"]))
        push!(lines, _md_row(["delta_eq chapter 1-2", _fmt(t.delta_eq_ch12)]))
        push!(lines, _md_row(["delta_eq chapter 3", _fmt(t.delta_eq_ch3)]))
        push!(lines, _md_row(["delta_eq full pendulum", _fmt(t.delta_eq_pendulum)]))
        push!(lines, _md_row(["second harmonic cos 2psi", _fmt(t.a2)]))
        push!(lines, _md_row(["Omega/omega_s chapter 3", _fmt_ratio(t.omega_ch3)]))
        push!(lines, _md_row(["Omega/omega_s pendulum", _fmt_ratio(t.omega_pendulum)]))
        push!(lines, _md_row(["Lambda", _fmt(t.Lambda)]))
        push!(lines, "")
        push!(lines, "NF `d · J` = $(_fmt(case.torus.d_dot)), `d_grid · J` = $(_fmt(case.torus.grid_dot)), torus mean = $(_fmt(case.torus.mean)). `J_s` = $(_fmt(case.Js)) from `δ_m² / (2 σ_δ²)` with `σ_δ` = $(_fmt(case.sigma)).")
        push!(lines, "")
        header = ["quantity", "formula", "NF torus", "degree-3 map", "exact track"]
        push!(lines, _md_row(header))
        push!(lines, _md_row(["---", "---", "---", "---", "---"]))
        labels = [
            "delta_eq",
            "delta_m half-range",
            "second harmonic cos 2psi",
            "second harmonic amplitude",
            "Q / Q_note",
            "mean square of delta",
            "phi_hat",
            "gamma_eff",
            "Delta E / E",
        ]
        cols = [_meas_cells(case.formula, case.Q_note),
                _meas_cells(case.nf, case.Q_note),
                _meas_cells(case.map, case.Q_note),
                _meas_cells(case.exact, case.Q_note)]
        for i in eachindex(labels)
            push!(lines, _md_row([labels[i], cols[1][i], cols[2][i], cols[3][i], cols[4][i]]))
        end
        push!(lines, "")
        push!(lines, "Map residual rms $(_fmt(case.map.rms)), exact residual rms $(_fmt(case.exact.rms)). Map `δ_eq` split $(_fmt(case.map.unc_delta_eq)), exact split $(_fmt(case.exact.unc_delta_eq)). Turns: map $(max(case.map.n - 1, 0)), exact $(max(case.exact.n - 1, 0)). Degree-3 image of the origin, delta component $(_fmt(case.fixed[6])).")
        push!(lines, "")
    end
    if !isempty(STOPS)
        push!(lines, "The run stopped on a physical disagreement. Numbers were not adjusted.")
        push!(lines, "")
    else
        push!(lines, "No physical STOP. The torus matches chapters 1–2. Where the cubic factor is resolved, the map and the exact track match `1 − φ̂²/32`.")
        push!(lines, "")
    end
    write(path, join(lines, "\n"))
    return nothing
end

# --- cases -------------------------------------------------------------------

function _track_column(kind, step, z0, nturns, ref, C, h, Q0, gamma, beta, label)
    local series
    try
        series = _run_series(step, z0, nturns, ref, C, h; abort_at = 0.05, label = label)
    catch err
        println(label, " failed: ", sprint(showerror, err))
        flush(stdout)
        ge, dE = kinematics(gamma, beta, NaN, NaN)
        return SeriesMeas(kind, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                          NaN, NaN, NaN, ge, dE, 0, true)
    end
    if series.aborted || !all(isfinite, series.delta) || length(series.delta) < 16
        ge, dE = kinematics(gamma, beta, NaN, NaN)
        return SeriesMeas(kind, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                          NaN, NaN, NaN, ge, dE, length(series.delta), true)
    end
    fit = harmonic_synchrotron(series.delta, series.phi_note, Q0)
    return _meas(kind, fit, gamma, beta; aborted = false)
end

function _nf_column(torus, gamma, beta, Q_nf)
    # δ_m of the torus is the half-range of δ along the synchrotron angle.
    # Q stays the small-amplitude normal-form tune: the torus has no φ̂² shift.
    fit = torus.fit
    de = torus.mean
    ge, dE = kinematics(gamma, beta, de, fit.delta2)
    return SeriesMeas("nf", de, fit.delta_m_half, fit.a1, fit.a2, fit.a2_cos, fit.a2_sin,
                      fit.a3, Q_nf, fit.delta2, NaN, fit.rms, abs(de - torus.d_dot), 0.0,
                      fit.unc_a2, ge, dE, fit.n, false)
end

function _run_case(title, blurb, targets, built, po, Jx, Jy, delta_m, Q_note, gamma, beta, C, eta0, eta1)
    nf = built.nf
    Js, σ = _js_for_delta_m(nf, delta_m)
    J = (Jx, Jy, Js)
    println("Torus ", title, "  Js = ", Js)
    flush(stdout)
    torus = _torus(nf, J)
    z0, ax, ay = _launch(po, delta_m, Jx, Jy)
    xb, pxb = dispersion_free_transverse(z0[1], z0[2], po.D4[1], po.D4[2], z0[6])
    yb, pyb = dispersion_free_transverse(z0[3], z0[4], po.D4[3], po.D4[4], z0[6])
    Jx_cs = courant_snyder_action(xb, pxb, po.x.beta, po.x.alpha)
    Jy_cs = courant_snyder_action(yb, pyb, po.y.beta, po.y.alpha)
    @printf("  launch ax %.6e ay %.6e  Jx_cs %.6e Jy_cs %.6e  (asked %.6e %.6e)\n",
            ax, ay, Jx_cs, Jy_cs, Jx, Jy)
    flush(stdout)
    if Jx == 0 && Jy == 0
        max(Jx_cs, Jy_cs) < 1e-18 ||
            error("dispersion-free launch action is not zero: $Jx_cs, $Jy_cs")
    else
        abs(Jx_cs - Jx) <= 1e-8 * Jx || error("launch Jx $Jx_cs is not the requested $Jx")
        abs(Jy_cs - Jy) <= 1e-8 * Jy || error("launch Jy $Jy_cs is not the requested $Jy")
    end
    nturns = _case_turns(Q_note)
    println("Map ", title, "  turns ", nturns)
    flush(stdout)
    map_m = _track_column("map", _map_step(built.M), z0, nturns, built.lattice.ref, C,
                          built.cav.harmonic, Q_note, gamma, beta, title * " map")
    println("Exact ", title, "  turns ", nturns)
    flush(stdout)
    ex_m = _track_column("exact", _exact_step(built.lattice, built.lattice.ref, C), z0, nturns,
                         built.lattice.ref, C, built.cav.harmonic, Q_note, gamma, beta,
                         title * " exact")
    Q_nf = abs(nf.mu[3]) / (2π)
    nf_m = _nf_column(torus, gamma, beta, Q_nf)
    formula = formula_meas(targets, gamma, beta, Q_note, delta_m)
    @printf("  delta_eq formula %.8e   NF %.8e   map %.8e   exact %.8e\n",
            formula.delta_eq, nf_m.delta_eq, map_m.delta_eq, ex_m.delta_eq)
    @printf("  a2_cos   formula %.8e   NF %.8e   map %.8e   exact %.8e\n",
            formula.a2_cos, nf_m.a2_cos, map_m.a2_cos, ex_m.a2_cos)
    flush(stdout)
    _judge_case!(title, targets, map_m, ex_m, Q_note, torus, eta0, eta1)
    return (
        title = title,
        blurb = blurb,
        targets = targets,
        formula = formula,
        nf = nf_m,
        map = map_m,
        exact = ex_m,
        torus = torus,
        Js = Js,
        sigma = σ,
        Q_note = Q_note,
        Q_nf = Q_nf,
        Jx_cs = Jx_cs,
        Jy_cs = Jy_cs,
        fixed = _fixed_point_delta(built.M),
    )
end

function _d_formula(xi_x, xi_y, eta0, eta1, C, sigma)
    return [
        2π * xi_x / (C * eta0),
        2π * xi_y / (C * eta0),
        -eta1 * sigma^2 / eta0,
    ]
end

function _judge_d!(tag, nf, formula)
    for (i, name) in enumerate(("x", "y", "s"))
        rel = abs(nf.d[i] - formula[i]) / max(abs(formula[i]), 1e-30)
        _push_ok = rel <= 1e-3
        @printf("  d_%s  nf %.8e  grid %.8e  formula %.8e  rel %.3e\n",
                name, nf.d[i], nf.d_grid[i], formula[i], rel)
        flush(stdout)
        _push_ok || stop!("$tag d_$name relative gap $rel. nf $(nf.d[i]), formula $(formula[i])")
        rel_g = abs(nf.d[i] - nf.d_grid[i]) / max(abs(nf.d[i]), 1e-30)
        rel_g <= 1e-3 ||
            stop!("$tag d_$name and d_grid differ by $rel_g")
    end
    return nothing
end

function main()
    _self_check()
    MAP_ORDER == 3 ||
        error("this step evaluates a degree-3 map; MAP_ORDER is $(MAP_ORDER)")
    path = closed_lattice_path()
    println("Loading ", path)
    flush(stdout)
    lattice = load_closed_lattice()
    _require_bare_ring(lattice)
    ref = lattice.ref
    C = path_length(lattice)
    inputs = _load_inputs()
    inputs.linear.phi_s_rad == 0.0 ||
        error("linear regime file phase is $(inputs.linear.phi_s_rad), expected the notes' 0")
    inputs.cubic.phi_s_rad == 0.0 ||
        error("cubic regime file phase is $(inputs.cubic.phi_s_rad), expected the notes' 0")
    inputs.linear.delta_m == inputs.cubic.delta_m ||
        error("the two regimes do not share delta_m")

    println("Magnetic one-turn factor")
    flush(stdout)
    maglat = magnetic_lattice(lattice)
    Fmag = ring_lie_factor(maglat)
    Mmag = factor_to_map(Fmag)
    pe = path_expansion(Mmag, ref, C)
    po = periodic_optics(Fmag.R, ref, C)
    po.x.stable && po.y.stable || error("magnetic periodic optics is unstable")
    gen = generator_chromaticity(Fmag, po)
    xi_x = XI_FACTOR * gen.x1.c
    xi_y = XI_FACTOR * gen.y1.c
    pe.eta0 < 0 || error("eta0 = $(pe.eta0) is not below transition")
    pe.eta0 * cos(BELL_PHI_S) > 0 || error("phi_s = π is not stable for this eta0")
    abs(pe.eta0 - inputs.eta0_json) / abs(inputs.eta0_json) < 1e-6 ||
        error("eta0 moved relative to 02_rf_regimes.json")
    abs(pe.eta1 - inputs.eta1_json) / abs(inputs.eta1_json) < 1e-6 ||
        error("eta1 moved relative to 02_rf_regimes.json")

    gamma = ref.gamma
    beta = ref.beta
    Λ = betatron_lambda(C, inputs.Jx, inputs.Jy, xi_x, xi_y)
    dm = inputs.linear.delta_m

    println("Inequalities")
    flush(stdout)
    header = Tuple{String,Float64}[]
    function note_ineq(q, v)
        push!(header, (q, v))
        println("  ", q, " = ", v)
        flush(stdout)
        return nothing
    end
    note_ineq("C_m", C)
    note_ineq("gamma", gamma)
    note_ineq("beta", beta)
    note_ineq("eta0", pe.eta0)
    note_ineq("eta1", pe.eta1)
    note_ineq("xi_x_generator", xi_x)
    note_ineq("xi_y_generator", xi_y)
    note_ineq("Jx", inputs.Jx)
    note_ineq("Jy", inputs.Jy)
    note_ineq("Lambda", Λ)
    note_ineq("Lambda_json_placeholder", inputs.Lambda_json)
    note_ineq("abs_eta1_delta_m_over_eta0", abs(pe.eta1 * dm / pe.eta0))
    note_ineq("delta_m", dm)
    note_ineq("abs_Lambda_over_eta0_delta_m", abs(Λ / (pe.eta0 * dm)))
    for reg in (inputs.linear, inputs.cubic)
        w = 2π * reg.Qs
        note_ineq(reg.name * "_phi_hat", reg.phi_hat)
        note_ineq(reg.name * "_phi_hat_squared", reg.phi_hat^2)
        note_ineq(reg.name * "_omega_s_Ts", w)
        note_ineq(reg.name * "_omega_over_betatron_x", w / po.x.mu)
        note_ineq(reg.name * "_omega_over_betatron_y", w / po.y.mu)
        note_ineq(reg.name * "_Qs", reg.Qs)
        Qs_bell = synchrotron_tune(ref, pe.eta0, reg.voltage_MV, reg.harmonic, BELL_PHI_S)
        note_ineq(reg.name * "_Qs_smooth_phi_pi", Qs_bell)
        note_ineq(reg.name * "_voltage_MV", reg.voltage_MV)
    end
    println("  |eta1 delta_m / eta0| = ", abs(pe.eta1 * dm / pe.eta0), " is much less than 1")
    println("  delta_m = ", dm, " is much less than 1")
    println("  |Lambda / (eta0 delta_m)| = ", abs(Λ / (pe.eta0 * dm)), " is much less than 1")
    println("  linear phi_hat = ", inputs.linear.phi_hat, " is much less than 1")
    println("  cubic phi_hat = ", inputs.cubic.phi_hat, " is not the linear condition")
    println("  cubic phi_hat^2 = ", inputs.cubic.phi_hat^2, " is much less than 1")
    for reg in (inputs.linear, inputs.cubic)
        w = 2π * reg.Qs
        println("  ", reg.name, " omega_s Ts = ", w, " is much less than 1")
        println("  ", reg.name, " omega_s Ts / betatron phase x = ", w / po.x.mu, " is much less than 1")
    end
    flush(stdout)

    probe_pi = _phase_probe(lattice, ref, C, po, inputs.linear, BELL_PHI_S, inputs.linear.Qs)
    probe_0 = _phase_probe(lattice, ref, C, po, inputs.linear, 0.0, inputs.linear.Qs)
    probes = (probe_pi, probe_0)
    if !(probe_pi.oscillates && !probe_pi.ran_away)
        stop!("phi_s = π did not oscillate inside the bucket. min $(probe_pi.minimum), max $(probe_pi.maximum), aborted $(probe_pi.aborted)")
    else
        println("Phase dictionary: phi_s = π oscillates and stays bounded")
        flush(stdout)
    end
    if !probe_0.ran_away
        stop!("phi_s = 0 did not leave the bucket. min $(probe_0.minimum), max $(probe_0.maximum)")
    else
        println("Phase dictionary: phi_s = 0 leaves the bucket")
        flush(stdout)
    end

    cases = []
    rows = Tuple{String,String,String,Float64}[]
    phase_ok = probe_pi.oscillates && !probe_pi.ran_away && probe_0.ran_away
    if phase_ok
        println("Linear-regime factor")
        flush(stdout)
        linear = _build_cavity_ring(lattice, inputs.linear.voltage_MV, inputs.linear.harmonic)
        dlin = _d_formula(xi_x, xi_y, pe.eta0, pe.eta1, C, linear.nf.sigma_delta_s)
        println("Normal-form d versus chapters 1-2, linear voltage")
        flush(stdout)
        _judge_d!("linear", linear.nf, dlin)
        fixed = _fixed_point_delta(linear.M)
        println("  |M(0)| delta = ", fixed[6])
        flush(stdout)
        abs(fixed[6]) > 1e-9 &&
            stop!("linear degree-3 map does not fix the origin in delta: $(fixed[6])")

        tA = make_targets(pe.eta0, pe.eta1, gamma, beta, dm, inputs.linear.phi_hat, 0.0, false)
        caseA = _run_case(
            "A linear J=0",
            "Linear voltage, `J_x = J_y = 0`. The formula is chapter 1: `⟨δ²⟩ = δ_m²/2`, `Λ = 0`.",
            tA, linear, po, 0.0, 0.0, dm, inputs.linear.Qs, gamma, beta, C, pe.eta0, pe.eta1)
        tB = make_targets(pe.eta0, pe.eta1, gamma, beta, dm, inputs.linear.phi_hat, Λ, false)
        caseB = _run_case(
            "B linear finite J",
            "Same linear voltage. `J_x` and `J_y` are the placeholder actions. `Λ` uses generator chromaticity.",
            tB, linear, po, inputs.Jx, inputs.Jy, dm, inputs.linear.Qs, gamma, beta, C, pe.eta0, pe.eta1)

        println("Cubic-phase factor")
        flush(stdout)
        cubic = _build_cavity_ring(lattice, inputs.cubic.voltage_MV, inputs.cubic.harmonic)
        dcub = _d_formula(xi_x, xi_y, pe.eta0, pe.eta1, C, cubic.nf.sigma_delta_s)
        println("Normal-form d versus chapters 1-2, cubic voltage")
        flush(stdout)
        _judge_d!("cubic", cubic.nf, dcub)
        fixed_c = _fixed_point_delta(cubic.M)
        abs(fixed_c[6]) > 1e-9 &&
            stop!("cubic degree-3 map does not fix the origin in delta: $(fixed_c[6])")
        tC = make_targets(pe.eta0, pe.eta1, gamma, beta, dm, inputs.cubic.phi_hat, 0.0, true)
        caseC = _run_case(
            "C cubic J=0",
            "Cubic-phase voltage, `J = 0`. The formula column uses `⟨δ²⟩ = (δ_m²/2) (1 − φ̂²/32)` and `Ω/ω_s = 1 − φ̂²/16`. The torus is still the chapter 1–2 value. The pendulum column is the full sine, kept as the size of the dropped `φ^5` term.",
            tC, cubic, po, 0.0, 0.0, dm, inputs.cubic.Qs, gamma, beta, C, pe.eta0, pe.eta1)
        cases = [caseA, caseB, caseC]
        for (tag, nf, formula) in (("linear", linear.nf, dlin), ("cubic", cubic.nf, dcub))
            for (i, name) in enumerate(("x", "y", "s"))
                _push!(rows, tag, "nf.d", "d_" * name, nf.d[i])
                _push!(rows, tag, "nf.d_grid", "d_" * name, nf.d_grid[i])
                _push!(rows, tag, "formula_ch12", "d_" * name, formula[i])
            end
            _push!(rows, tag, "nf", "Qs", abs(nf.mu[3]) / (2π))
            _push!(rows, tag, "nf", "sigma_delta", nf.sigma_delta_s)
        end
        for case in cases
            for (col, m) in (("formula", case.formula), ("nf", case.nf),
                             ("map", case.map), ("exact", case.exact))
                _push_meas!(rows, case.title, col, m)
            end
            _push!(rows, case.title, "target", "delta_eq_ch12", case.targets.delta_eq_ch12)
            _push!(rows, case.title, "target", "delta_eq_ch3", case.targets.delta_eq_ch3)
            _push!(rows, case.title, "target", "delta_eq_pendulum", case.targets.delta_eq_pendulum)
            _push!(rows, case.title, "nf", "d_dot_J", case.torus.d_dot)
            _push!(rows, case.title, "nf", "d_grid_dot_J", case.torus.grid_dot)
        end
    end

    mkpath(OUT_DIR)
    report = joinpath(OUT_DIR, "04_equilibrium.md")
    numbers = joinpath(OUT_DIR, "04_numbers.csv")
    _write_numbers(numbers, rows)
    _write_report(report, nothing, probes, cases, header)
    println("Wrote ", report)
    println("Wrote ", numbers)
    if isempty(STOPS)
        println("No physical STOP.")
    else
        println("STOP count: ", length(STOPS))
        for s in STOPS
            println("  ", s)
        end
    end
    flush(stdout)
    isempty(STOPS) || exit(2)
    return nothing
end

main()
