# Magnetic ring study on the current BELL API.
# The BELL tree is read-only. This module loads examples/lattices/magnetic_2.bell
# and writes a closed copy under this study. Later steps read that copy:
# linear optics, slip factor, and the stationary RF regimes chosen from them.

module MagneticRing

using Bell

"""Absolute path of the BELL example lattice. Read-only."""
const BELL_LATTICE = "C:/GIT_REPOS/BELL/examples/lattices/magnetic_2.bell"

"""Sector bending angle that closes sixteen BM magnets on 2π."""
const TARGET_BM_ANGLE = Float64(π) / 8

_study_root() = normpath(joinpath(@__DIR__, ".."))

"""Path of the closed lattice written by [`write_closed_lattice`](@ref)."""
closed_lattice_path() = joinpath(_study_root(), "lattices", "magnetic_2_closed.bell")

"""Load the BELL magnetic_2 example."""
load_magnetic_2() = load_structure(BELL_LATTICE)

"""
    close_bm_angles(lattice; angle = TARGET_BM_ANGLE) -> Lattice

Copy `lattice` with every sector bend named `BM` set to `angle`.
Drifts named `OBM` and every other element are unchanged.
"""
function close_bm_angles(lattice::Lattice; angle::Real = TARGET_BM_ANGLE)
    ang = Float64(angle)
    elements = AbstractElement[_retarget_bm(e, ang) for e in lattice.elements]
    return Lattice("magnetic_2_closed", elements, lattice.ref, lattice.harmonic)
end

function _retarget_bm(e::SBend, angle::Float64)
    e.name == "BM" || return e
    return SBend(e.name, e.L, angle, e.k1, e.e1, e.e2, e.gap, e.fint)
end
_retarget_bm(e::AbstractElement, ::Float64) = e

"""
    element_parameters(e) -> NamedTuple

`type`, `length`, `angle`, `k1`, and `k2`. A field the element does not carry
is `nothing` (a sector bend with `k1 = 0` keeps that zero).
"""
function element_parameters(e::AbstractElement)
    return (
        name = e.name,
        type = kind_string(e),
        length = elem_length(e),
        angle = e isa SBend ? e.angle : nothing,
        k1 = (e isa SBend || e isa Quadrupole) ? e.k1 : nothing,
        k2 = e isa Sextupole ? e.k2 : nothing,
    )
end

"""Distinct elements in first-appearance order, with a multiplicity."""
function unique_rows(lattice::Lattice)
    rows = NamedTuple[]
    for e in unique_elements(lattice)
        n = count(isequal(e), lattice.elements)
        push!(rows, merge(element_parameters(e), (count = n,)))
    end
    return rows
end

"""
    source_layout(lattice) -> NamedTuple

Require 16 sector bends named BM and 4 drifts named OBM, all of one length,
and a single BM angle. Throws if the example does not have that layout.
"""
function source_layout(lattice::Lattice)
    bms = [e for e in lattice.elements if e.name == "BM"]
    obms = [e for e in lattice.elements if e.name == "OBM"]
    length(bms) == 16 || error("expected 16 elements named BM, found $(length(bms))")
    length(obms) == 4 || error("expected 4 elements named OBM, found $(length(obms))")
    all(e -> e isa SBend, bms) || error("an element named BM is not a sector bend")
    all(e -> e isa Drift, obms) || error("an element named OBM is not a drift")
    count(e -> e isa SBend, lattice.elements) == 16 ||
        error("the ring has a sector bend that is not named BM")
    L = bms[1].L
    all(e -> e.L == L, bms) || error("BM lengths are not all equal")
    all(e -> e.L == L, obms) || error("OBM length differs from the BM arc length")
    angles = unique(e.angle for e in bms)
    length(angles) == 1 || error("BM angles are not all equal: $(collect(angles))")
    return (
        n_bm = 16,
        n_obm = 4,
        n_elements = length(lattice.elements),
        length_m = L,
        bm_angle = angles[1],
    )
end

"""True when `closed` differs from `source` only by the BM bending angle."""
function only_bm_angle_changed(source::Lattice, closed::Lattice)
    length(source.elements) == length(closed.elements) || return false
    source.harmonic == closed.harmonic || return false
    isequal(source.ref, closed.ref) || return false
    for (a, b) in zip(source.elements, closed.elements)
        if a.name == "BM"
            (b.name == "BM" && a isa SBend && b isa SBend) || return false
            b.L == a.L && b.k1 == a.k1 && b.e1 == a.e1 && b.e2 == a.e2 &&
                b.gap == a.gap && b.fint == a.fint || return false
        else
            isequal(a, b) || return false
        end
    end
    return true
end

"""Element-wise equality, including the reference particle and the harmonic."""
function same_elements(a::Lattice, b::Lattice)
    length(a.elements) == length(b.elements) || return false
    all(isequal(x, y) for (x, y) in zip(a.elements, b.elements)) || return false
    return a.harmonic == b.harmonic && isequal(a.ref, b.ref)
end

"""
    ring_readings(lattice; ds = 0.02) -> NamedTuple

Path length, total bend, and tunes.

`Qx` and `Qy` are the accumulated tunes from `optics_along` (integer part
included). `Qx_frac` and `Qy_frac` are μ/2π from `periodic_optics` of
`ring_matrix`. `closure_Qx` and `closure_Qy` are the fractional-tune residuals
of that transport against the one-turn matrix.
"""
function ring_readings(lattice::Lattice; ds::Real = 0.02)
    C = path_length(lattice)
    bend = total_bend_angle(lattice)
    R = ring_matrix(lattice)
    po = periodic_optics(R, lattice.ref, C)
    along = optics_along(lattice, po; ds = ds)
    po.x.stable && po.y.stable || error("periodic optics of $(lattice.name) is unstable")
    return (
        path_length = C,
        total_bend_angle = bend,
        bend_minus_2pi = bend - 2π,
        Qx = along.Qx,
        Qy = along.Qy,
        Qx_frac = po.x.mu / (2π),
        Qy_frac = po.y.mu / (2π),
        closure_Qx = along.closure.Qx_frac,
        closure_Qy = along.closure.Qy_frac,
        symplectic_error = symplectic_error(R),
    )
end

"""
    write_closed_lattice(path = closed_lattice_path()) -> path

Copy `BELL_LATTICE` and set the single BM card's angle to a decimal that
parses as [`TARGET_BM_ANGLE`](@ref). The destination must stay inside this
study. The BELL file is not opened for writing.
"""
function write_closed_lattice(path::AbstractString = closed_lattice_path())
    dst = normpath(path)
    _refuse_outside_study(dst)
    literal = _angle_literal(TARGET_BM_ANGLE)
    text = _with_bm_angle(read(BELL_LATTICE, String), literal)
    mkpath(dirname(dst))
    write(dst, text)
    loaded = load_structure(dst)
    got = source_layout(loaded).bm_angle
    got == TARGET_BM_ANGLE ||
        error("reloaded BM angle $got is not $(TARGET_BM_ANGLE)")
    return dst
end

"""Load the closed lattice written under this study."""
function load_closed_lattice()
    path = closed_lattice_path()
    isfile(path) || error("closed lattice is missing: $path. Run scripts/01_lattice.jl")
    return load_structure(path)
end

function _refuse_outside_study(dst::AbstractString)
    root = lowercase(normpath(_study_root()))
    bell = lowercase(normpath("C:/GIT_REPOS/BELL"))
    place = lowercase(dst)
    startswith(place, root) || error("refusing to write outside the magnetic_ring study: $dst")
    startswith(place, bell) && error("refusing to write into BELL: $dst")
    return nothing
end

function _angle_literal(angle::Float64)
    literal = string(angle)
    parse(Float64, literal) == angle || error("angle literal does not round-trip: $literal")
    return literal
end

_code(line::AbstractString) = begin
    i = findfirst('#', line)
    return i === nothing ? line : line[1:prevind(line, i)]
end

function _with_bm_angle(text::AbstractString, literal::AbstractString)
    nl = occursin("\r\n", text) ? "\r\n" : "\n"
    lines = split(text, r"\r?\n", keepempty = true)
    n = 0
    kept = String[]
    for line in lines
        tokens = split(strip(_code(line)))
        if length(tokens) >= 2 && tokens[1] == "BM" && lowercase(tokens[2]) == "sbend"
            line = _replace_angle_token(line, literal)
            n += 1
        end
        # Keep the lattice cards and the author line. Drop other header comments
        # from the example file so the study copy stays a plain structure file.
        stripped = lstrip(line)
        if startswith(stripped, "#") && !occursin("Author:", line)
            continue
        end
        push!(kept, line)
    end
    n == 1 || error("expected one BM sector-bend card in the source lattice, found $n")
    header = [
        "# Bell structure file",
        "# SCT_2 magnetic ring: every BM angle is set to pi/8.",
        "# OBM drifts are unchanged. Copied from BELL examples/lattices/magnetic_2.bell.",
    ]
    return join(vcat(header, kept), nl)
end

function _replace_angle_token(line::AbstractString, literal::AbstractString)
    re = r"angle\s*=\s*[^\s#]+"
    count(re, line) == 1 || error("BM card does not have one angle: $line")
    return replace(line, re => "angle = " * literal; count = 1)
end

"""Angle token on the BM card of a structure file, for the record."""
function source_bm_angle_token(path::AbstractString = BELL_LATTICE)
    text = read(path, String)
    for line in split(text, r"\r?\n")
        tokens = split(strip(_code(line)))
        if length(tokens) >= 2 && tokens[1] == "BM" && lowercase(tokens[2]) == "sbend"
            m = match(r"angle\s*=\s*([^\s#]+)", line)
            m === nothing && error("BM card has no angle token")
            return m.captures[1]
        end
    end
    error("no BM sector-bend card in $path")
end

"""
    sample_dispersion(lattice, po; ds = 0.02) -> NamedTuple

Horizontal and vertical dispersion on the same slices as `optics_along`.

`optics_along` keeps only `D_x` and `D_x'`. The slice count is the same rule:
`n = ceil(L/ds)` pieces, each advanced by `linear_matrix` of length `L/n`.
`po` is a `periodic_optics` result (the start value `D4`).
"""
function sample_dispersion(lattice::Lattice, po; ds::Real = 0.02)
    ref = lattice.ref
    eta = vcat(Float64.(po.D4), 0.0, 1.0)
    Dx = [eta[1]]
    Dpx = [eta[2]]
    Dy = [eta[3]]
    Dpy = [eta[4]]
    for e in lattice.elements
        L = elem_length(e)
        n = L == 0 ? 0 : max(1, ceil(Int, L / ds))
        for _ in 1:n
            eta = linear_matrix(e, ref, L / n) * eta
            push!(Dx, eta[1])
            push!(Dpx, eta[2])
            push!(Dy, eta[3])
            push!(Dpy, eta[4])
        end
    end
    return (Dx = Dx, Dpx = Dpx, Dy = Dy, Dpy = Dpy)
end

"""
    revolution_period(ref, path_length) -> T_s [s]

Revolution period `C / (β c)` of the reference particle.
"""
revolution_period(ref::RefParticle, path_length::Real) = path_length / (ref.beta * c_light)

"""
    stationary_phi_s(eta0) -> phi_s [rad]

Synchronous phase of a stationary bucket. The reference particle is not
accelerated, so `sin(phi_s) = 0`. Small oscillations are stable for a positive
voltage when `cos(phi_s) = -sign(eta0)`: `0` below transition, `π` above it.
"""
function stationary_phi_s(eta0::Real)
    iszero(eta0) && error("eta0 is zero; the stationary-bucket phase is not defined")
    return eta0 < 0 ? 0.0 : float(π)
end

"""
    synchrotron_tune_note(ref, eta0, voltage_MV, harmonic) -> Q_s

Small-amplitude synchrotron tune from the notes,

    Q_s = sqrt( h |η0| Z V / (2 π β² E) ),

with `V` in MV and `E = γ m` in MeV. `Q_s = ω_s T_s / 2π`. The revolution
period cancels, and the stable-phase sign is already taken as `|η0|`.
"""
function synchrotron_tune_note(ref::RefParticle, eta0::Real, voltage_MV::Real, harmonic::Integer)
    voltage_MV > 0 || error("voltage must be positive, got $voltage_MV MV")
    harmonic > 0 || error("harmonic must be positive, got $harmonic")
    iszero(eta0) && error("eta0 is zero")
    Z = ref.particle.charge
    Z > 0 || error("charge number must be positive, got $Z")
    E = ref.gamma * ref.mass_MeV
    a = harmonic * abs(eta0) * Z * voltage_MV / (2π * ref.beta^2 * E)
    return sqrt(a)
end

"""
    voltage_for_phi_hat(ref, eta0, delta_m, harmonic, phi_hat) -> V [MV]

Voltage that realises a phase half-range `phi_hat` at momentum half-range
`delta_m` for a stationary bucket:

    φ̂ = h |η0| δ_m / Q_s.
"""
function voltage_for_phi_hat(ref::RefParticle, eta0::Real, delta_m::Real, harmonic::Integer, phi_hat::Real)
    delta_m > 0 || error("delta_m must be positive")
    phi_hat > 0 || error("phi_hat must be positive")
    Z = ref.particle.charge
    E = ref.gamma * ref.mass_MeV
    return harmonic * abs(eta0) * delta_m^2 * 2π * ref.beta^2 * E / (phi_hat^2 * Z)
end

"""Phase half-range `φ̂ = h |η0| δ_m / Q_s`."""
phase_half_range(eta0::Real, delta_m::Real, harmonic::Integer, Qs::Real) =
    harmonic * abs(eta0) * delta_m / Qs

"""`ω_s T_s = 2π Q_s`."""
omega_s_Ts(Qs::Real) = 2π * Qs

"""
    action_for_amplitude(beta, amplitude) -> J

Courant–Snyder action whose peak displacement is `amplitude` at this beta:
`J = a² / (2 β)`.
"""
function action_for_amplitude(beta::Real, amplitude::Real)
    beta > 0 || error("beta must be positive, got $beta")
    amplitude >= 0 || error("amplitude must be non-negative")
    return amplitude^2 / (2 * beta)
end

"""
    betatron_lambda(path_length, Jx, Jy, xi_x, xi_y) -> Λ

Relative betatron path lengthening at fixed actions,

    Λ = -(2π / C) (J_x ξ_x + J_y ξ_y).

`ξ` is `dQ/dδ`. The actions are chosen by the caller; this function does not
pick an amplitude.
"""
function betatron_lambda(path_length::Real, Jx::Real, Jy::Real, xi_x::Real, xi_y::Real)
    path_length > 0 || error("path length must be positive")
    return -2π / path_length * (Jx * xi_x + Jy * xi_y)
end

# Concrete stand-ins for the notes' "much less than one". Recorded with the
# regimes. Not a fit to a measured voltage.
const RF_CUTS = (
    eta_ratio_max = 0.02,
    lambda_ratio_max = 0.02,
    delta_m_max = 1.0e-2,
    omega_s_Ts_max = 0.02,
    omega_over_betatron_max = 0.05,
    linear_phi_hat = 0.05,
    cubic_phi_hat = 0.30,
)
const DELTA_M_CANDIDATES = (1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5)
const AMPLITUDE_CANDIDATES = (1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5)

"""`2π{ν}` and `2π` times the distance from `{ν}` to the nearest integer."""
function _mixing_scales(tune::Real)
    f = mod(float(tune), 1.0)
    return (two_pi_frac = 2π * f, two_pi_nearest = 2π * min(f, 1.0 - f))
end

function _select_delta_m(eta0, eta1, harmonic, phi_hat, mix_min)
    abs_ratio = abs(eta1 / eta0)
    for dm in DELTA_M_CANDIDATES
        dm <= RF_CUTS.delta_m_max || continue
        abs_ratio * dm <= RF_CUTS.eta_ratio_max || continue
        Qs = harmonic * abs(eta0) * dm / phi_hat
        w = omega_s_Ts(Qs)
        w <= RF_CUTS.omega_s_Ts_max || continue
        w <= RF_CUTS.omega_over_betatron_max * mix_min || continue
        return dm
    end
    error("no momentum half-range meets the linear-regime cuts at h=$harmonic")
end

function _select_amplitude(path_length, betx_max, bety_max, xi_x, xi_y, eta0, delta_m)
    for a in AMPLITUDE_CANDIDATES
        Jx = action_for_amplitude(betx_max, a)
        Jy = action_for_amplitude(bety_max, a)
        Λ = betatron_lambda(path_length, Jx, Jy, xi_x, xi_y)
        if abs(Λ / (eta0 * delta_m)) <= RF_CUTS.lambda_ratio_max
            return (amplitude_m = a, Jx = Jx, Jy = Jy, Lambda = Λ,
                    beta_x_m = betx_max, beta_y_m = bety_max,
                    xi_x = xi_x, xi_y = xi_y)
        end
    end
    error("no placeholder amplitude keeps |Lambda/(eta0*delta_m)| within the cut")
end

function _regime(name, ref, eta0, eta1, harmonic, delta_m, phi_hat, phi_s, T_s, Λ, Qx, Qy)
    V = voltage_for_phi_hat(ref, eta0, delta_m, harmonic, phi_hat)
    Qs = synchrotron_tune_note(ref, eta0, V, harmonic)
    Qs_from_phase = harmonic * abs(eta0) * delta_m / phi_hat
    abs(Qs - Qs_from_phase) <= 1e-12 * max(Qs, 1e-16) ||
        error("voltage for $name does not reproduce the target phase half-range")
    w = omega_s_Ts(Qs)
    sx = _mixing_scales(Qx)
    sy = _mixing_scales(Qy)
    ph = phase_half_range(eta0, delta_m, harmonic, Qs)
    cubic = name == "cubic_phase" ? (1 - ph^2 / 32) : 1.0
    delta_eq = -(eta1 * delta_m^2 / 2 * cubic + Λ) / eta0
    return (
        name = name,
        harmonic = harmonic,
        voltage_MV = V,
        delta_m = delta_m,
        phi_s_rad = phi_s,
        sin_phi_s = sin(phi_s),
        cos_phi_s = cos(phi_s),
        phi_hat = ph,
        phi_hat_squared = ph^2,
        Qs = Qs,
        omega_s_Ts = w,
        omega_rf_rad_s = 2π * harmonic / T_s,
        abs_eta1_delta_m_over_eta0 = abs(eta1 * delta_m / eta0),
        abs_Lambda_over_eta0_delta_m = abs(Λ / (eta0 * delta_m)),
        two_pi_frac_nu_x = sx.two_pi_frac,
        two_pi_frac_nu_y = sy.two_pi_frac,
        two_pi_nearest_nu_x = sx.two_pi_nearest,
        two_pi_nearest_nu_y = sy.two_pi_nearest,
        omega_s_Ts_over_two_pi_frac_nu_x = w / sx.two_pi_frac,
        omega_s_Ts_over_two_pi_frac_nu_y = w / sy.two_pi_frac,
        omega_s_Ts_over_two_pi_nearest_nu_x = w / sx.two_pi_nearest,
        omega_s_Ts_over_two_pi_nearest_nu_y = w / sy.two_pi_nearest,
        formula_delta_eq = delta_eq,
    )
end

"""
    build_rf_regimes(ref, path_length, eta0, eta1, Qx, Qy, betx_max, bety_max, xi_x, xi_y, harmonic)

Two stationary buckets that share `δ_m` and the lattice harmonic.

The linear regime targets `φ̂ = 0.05`. The cubic-phase regime targets
`φ̂ = 0.30` (`φ̂² = 0.09`), so the small-phase condition is not the one in
force and the cubic-chapter condition `φ̂² ≪ 1` is. `δ_m` is the largest
entry of a fixed list that meets the linear-regime cuts. The betatron action
is a placeholder: peak displacement 1 mm at `β_max`, reduced only if
`|Λ/(η0 δ_m)|` misses its cut. `ξ` is an input, not computed here.
"""
function build_rf_regimes(ref::RefParticle, path_length::Real, eta0::Real, eta1::Real,
                          Qx::Real, Qy::Real, betx_max::Real, bety_max::Real,
                          xi_x::Real, xi_y::Real, harmonic::Integer)
    iszero(eta0) && error("eta0 is zero")
    harmonic > 0 || error("harmonic must be positive")
    T_s = revolution_period(ref, path_length)
    sx = _mixing_scales(Qx)
    sy = _mixing_scales(Qy)
    mix_min = minimum((sx.two_pi_frac, sy.two_pi_frac, sx.two_pi_nearest, sy.two_pi_nearest))
    mix_min > 0 || error("a fractional tune is an integer; betatron-phase averaging is not justified")
    dm = _select_delta_m(eta0, eta1, harmonic, RF_CUTS.linear_phi_hat, mix_min)
    hold = _select_amplitude(path_length, betx_max, bety_max, xi_x, xi_y, eta0, dm)
    phi_s = stationary_phi_s(eta0)
    linear = _regime("linear", ref, eta0, eta1, harmonic, dm, RF_CUTS.linear_phi_hat,
                     phi_s, T_s, hold.Lambda, Qx, Qy)
    cubic = _regime("cubic_phase", ref, eta0, eta1, harmonic, dm, RF_CUTS.cubic_phi_hat,
                    phi_s, T_s, hold.Lambda, Qx, Qy)
    linear.phi_hat <= 0.10 || error("linear phase half-range is not small: $(linear.phi_hat)")
    cubic.phi_hat_squared <= 0.10 || error("cubic phase half-range squared is not small: $(cubic.phi_hat_squared)")
    cubic.phi_hat > linear.phi_hat * 4 || error("cubic and linear phase half-ranges are not separated")
    for reg in (linear, cubic)
        reg.omega_s_Ts <= RF_CUTS.omega_s_Ts_max ||
            error("$(reg.name): omega_s T_s = $(reg.omega_s_Ts)")
        reg.abs_eta1_delta_m_over_eta0 <= RF_CUTS.eta_ratio_max ||
            error("$(reg.name): |eta1 delta_m/eta0| misses the cut")
        reg.abs_Lambda_over_eta0_delta_m <= RF_CUTS.lambda_ratio_max ||
            error("$(reg.name): |Lambda/(eta0 delta_m)| misses the cut")
        reg.delta_m <= RF_CUTS.delta_m_max || error("$(reg.name): delta_m misses the cut")
        abs(reg.sin_phi_s) <= 1e-14 || error("$(reg.name): stationary bucket has sin phi_s != 0")
    end
    return (
        cuts = RF_CUTS,
        harmonic = harmonic,
        phi_s_rad = phi_s,
        delta_m = dm,
        T_s = T_s,
        placeholder = hold,
        linear = linear,
        cubic = cubic,
    )
end

export BELL_LATTICE, TARGET_BM_ANGLE, closed_lattice_path,
       load_magnetic_2, close_bm_angles, write_closed_lattice, load_closed_lattice,
       element_parameters, unique_rows, source_layout, only_bm_angle_changed,
       same_elements, ring_readings, source_bm_angle_token,
       sample_dispersion, revolution_period, stationary_phi_s, synchrotron_tune_note,
       voltage_for_phi_hat, phase_half_range, omega_s_Ts, action_for_amplitude,
       betatron_lambda, build_rf_regimes, RF_CUTS

end # module
