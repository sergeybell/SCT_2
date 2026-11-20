# Step 6. Dipole fringe on the own sextupole solution.
#
#   julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/06_fringe.jl
#
# Rebuild the step-5 placement (optics_along at ds = 0.02 m, then
# place_sextupole_families) and set k2 from output/05_sext_strengths.json,
# object k2_per_m3.own. Strengths are not solved again.
#
# Three FringeModel rows, edges on dipoles only:
#   NoEdge                         FringeModel()
#   HardEdge()                     FringeModel(dipole = HardEdge())
#   SoftEdge(LinearRamp(0.05))     FringeModel(dipole = SoftEdge(profile = LinearRamp(0.05)))
# The ramp width 0.05 m is an assumption of this step.
#
# Q, ξ, η1, and the linear-regime δ_eq formula come from ring_lie_factor(lat, model).
# ξ is generator_chromaticity of that factor. The normalising matrix is built
# from periodic_optics of the same factor, including when the model changes
# the linear optics. δ_m, Jx, and Jy stay the linear-regime values of steps 4 and 5.
#
# exact_track_element has no EdgeMap method. A short track walks fringe_sequence:
# EdgeMap → apply_factor of the Lie factor of edge_map; magnets → exact_track_element;
# thin cavity → rf_kick. --turns N replaces that track.

import Pkg

const STUDY_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(STUDY_DIR)

using LinearAlgebra
using Printf
using Bell
using MagneticRing

include(joinpath(STUDY_DIR, "src", "generator_chromaticity.jl"))
include(joinpath(STUDY_DIR, "src", "readings.jl"))
using .GeneratorChromaticity
using .LieReadings

const OUT_DIR = joinpath(STUDY_DIR, "output")
const REGIME_JSON = joinpath(OUT_DIR, "02_rf_regimes.json")
const STRENGTH_JSON = joinpath(OUT_DIR, "05_sext_strengths.json")
const DS = 0.02
const TRACK_DS = 0.1
const BELL_PHI_S = Float64(π)
# Step 2b: ξ1 = −c / (2π), c the coefficient of J δ.
const XI_FACTOR = -1 / (2π)
# LinearRamp width [m]. An assumption: the lattice card does not carry it.
const RAMP_WIDTH = 0.05
const MATRIX_TOL = 1.0e-8
const NORM_TOL = 1.0e-8
const SLIP_TOL = 1.0e-8
# Generator ξ and chromaticity_orders already differ in the sixth digit on
# this ring (step 5). A gap above this is a new disagreement.
const XI_READ_TOL = 1.0e-3
const TRACK_REL = 0.005
const TRACK_ABS = 1.0e-11
# A fringe shift larger than this leaves the integer part of Q ambiguous.
const TUNE_SHIFT_MAX = 0.4

const STOPS = String[]

function stop!(msg)
    push!(STOPS, msg)
    println("STOP: ", msg)
    flush(stdout)
    return nothing
end

# --- files -------------------------------------------------------------------

function _object_after(text::AbstractString, key::AbstractString)
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

function _json_number(block::AbstractString, key::AbstractString)
    m = match(Regex("\"$(key)\"\\s*:\\s*([-+0-9.eE]+)"), block)
    m === nothing && error("JSON object has no number \"$key\"")
    return parse(Float64, m.captures[1])
end

function _json_string(block::AbstractString, key::AbstractString)
    m = match(Regex("\"$(key)\"\\s*:\\s*\"([^\"]*)\""), block)
    m === nothing && error("JSON object has no string \"$key\"")
    return m.captures[1]
end

function _load_linear(path)
    isfile(path) || error("regime file is missing: $path")
    text = read(path, String)
    lin = _object_after(text, "linear")
    hold = _object_after(text, "placeholder_action")
    return (
        voltage_MV = _json_number(lin, "voltage_MV"),
        harmonic = Int(round(_json_number(lin, "harmonic"))),
        delta_m = _json_number(lin, "delta_m"),
        Qs = _json_number(lin, "Qs"),
        phi_hat = _json_number(lin, "phi_hat"),
        notes_phi_s = _json_number(lin, "phi_s_rad"),
        Jx = _json_number(hold, "Jx"),
        Jy = _json_number(hold, "Jy"),
    )
end

function _family_indices(text::AbstractString, family::AbstractString)
    m = match(Regex("\\{[^{}]*\"name\"\\s*:\\s*\"$(family)\"[^{}]*\\}"), text)
    m === nothing && error("strength file has no family $family")
    nums = match(r"\"element_index\"\s*:\s*\[([^\]]*)\]", m.match)
    nums === nothing && error("family $family has no element_index")
    out = Int[]
    for part in split(nums.captures[1], ",")
        s = strip(part)
        isempty(s) && continue
        push!(out, parse(Int, s))
    end
    return out
end

function _load_own(path)
    isfile(path) || error("strength file is missing: $path")
    text = read(path, String)
    _json_string(text, "primary") == "own" ||
        error("step 5 primary is not own")
    abs(_json_number(text, "optics_ds_m") - DS) > 1.0e-12 &&
        error("step 5 placement step is not $(DS) m")
    own = match(r"\"own\"\s*:\s*\{([^{}]*)\}", text)
    own === nothing && error("strength file has no k2_per_m3.own object")
    k = Dict{String,Float64}()
    idx = Dict{String,Vector{Int}}()
    for f in SEXT_FAMILIES
        mk = match(Regex("\"$(f)\"\\s*:\\s*([-+0-9.eE]+)"), own.captures[1])
        mk === nothing && error("own strengths have no $f")
        k[f] = parse(Float64, mk.captures[1])
        idx[f] = _family_indices(text, f)
    end
    return (k2 = k, index = idx)
end

# --- checks and formula ------------------------------------------------------

function _self_check()
    abs(XI_FACTOR + 1 / (2π)) < 1.0e-15 || error("chromaticity factor is not -1/(2π)")
    abs(RAMP_WIDTH - 0.05) < 1.0e-15 || error("ramp width is not the stated assumption")
    LinearRamp(RAMP_WIDTH).width == RAMP_WIDTH || error("LinearRamp did not keep the width")
    # δ_eq = −(η1 δ_m²/2 + Λ) / η0
    abs(_delta_eq(-0.5, 1.0, 0.1, 0.0) - 0.01) < 1.0e-14 || error("delta_eq formula")
    abs(betatron_lambda(2π, 1.0, 0.0, 1.0, 0.0) + 1.0) < 1.0e-14 || error("Lambda formula")
    abs(_tune_shift(0.2, 0.1) - 0.1 / 2π) < 1.0e-14 || error("tune shift")
    abs(_tune_shift(0.1, 2π - 0.1) - 0.2 / 2π) < 1.0e-14 || error("tune shift across 0")
    return nothing
end

function _delta_eq(eta0, eta1, delta_m, Lambda)
    eta0 == 0 && error("eta0 is zero")
    return -(eta1 * delta_m^2 / 2 + Lambda) / eta0
end

function _gamma_shift(gamma, beta, delta_eq, delta2)
    b2 = float(beta)^2
    g = float(gamma)
    g_eff = g + b2 * g * delta_eq + b2 * delta2 / (2 * g)
    dE = b2 * (delta_eq + delta2 / (2 * g^2))
    return (gamma_eff = g_eff, dE_over_E = dE)
end

"Shortest fringe shift of a tune, in turns. Both phases are in radians."
function _tune_shift(mu, mu_ref)
    d = float(mu) - float(mu_ref)
    if d > π
        d -= 2π
    elseif d < -π
        d += 2π
    end
    return d / 2π
end

function _frac_diff(a, b)
    d = float(a) - float(b)
    return d - round(d)
end

function _turns(Qs)
    i = findfirst(==("--turns"), ARGS)
    if i !== nothing
        i < length(ARGS) || error("--turns needs an integer")
        n = parse(Int, ARGS[i + 1])
        n >= 16 || error("--turns must be at least 16, got $n")
        return n
    end
    Qs > 0 || error("synchrotron tune is not positive")
    return max(64, round(Int, 1 / Qs))
end

_turns_fixed() = findfirst(==("--turns"), ARGS) !== nothing

function _require_bare(lattice)
    ref = lattice.ref
    ref.particle.name == "deuteron" || error("reference particle is $(ref.particle.name)")
    ref.Ekin_MeV == 270.0 || error("kinetic energy is $(ref.Ekin_MeV) MeV")
    abs(total_bend_angle(lattice) - 2π) < 1.0e-9 ||
        error("total bend is not 2π: $(total_bend_angle(lattice))")
    any(e -> e isa RFCavity, lattice.elements) &&
        error("the closed lattice already contains an RF cavity")
    for e in lattice.elements
        e isa Sextupole && e.k2 != 0 &&
            error("sextupole $(e.name) has k2 = $(e.k2); placement starts from k2 = 0")
    end
    return nothing
end

function _norm_ok(g)
    return g.symplectic_error_A <= NORM_TOL &&
           g.symplectic_error_shear <= NORM_TOL &&
           g.dispersion_column <= NORM_TOL &&
           g.shear_dispersion_column <= NORM_TOL &&
           g.rotation_residual <= NORM_TOL &&
           g.transverse_coupling <= NORM_TOL &&
           g.delta_row <= NORM_TOL &&
           abs(g.delta_delta - 1) <= NORM_TOL &&
           g.homological_residual <= 1.0e-6 &&
           g.max_l_f3 <= 1.0e-8 &&
           g.max_l_f4 <= 1.0e-8
end

# --- lattice -----------------------------------------------------------------

function _place(lattice)
    ref = lattice.ref
    C = path_length(lattice)
    R = ring_matrix(lattice)
    po = periodic_optics(R, ref, C)
    po.x.stable && po.y.stable || error("periodic optics of the closed lattice is unstable")
    opt = optics_along(lattice, po; ds = DS)
    opt.closure.betx < 1.0e-8 || error("optics_along beta_x did not close")
    opt.closure.D < 1.0e-8 || error("optics_along dispersion did not close")
    println("place_sextupole_families")
    flush(stdout)
    placed, fam_idx, rows, notes = place_sextupole_families(lattice, opt)
    return (lattice = placed, fam_idx = fam_idx, rows = rows, notes = notes, po = po, along = opt)
end

function _check_placement(placed, own)
    for f in SEXT_FAMILIES
        got = placed.fam_idx[f]
        exp = own.index[f]
        got == exp ||
            error("family $f indices $got differ from step 5 $exp")
        for i in got
            e = placed.lattice.elements[i]
            (e isa Sextupole && e.name == f) ||
                error("index $i is $(e.name), expected sextupole $f")
        end
    end
    return nothing
end

function _models()
    ramp = LinearRamp(RAMP_WIDTH)
    return (
        (key = "NoEdge", dipole = "NoEdge", model = FringeModel()),
        (key = "HardEdge", dipole = "HardEdge()", model = FringeModel(dipole = HardEdge())),
        (key = "SoftEdge", dipole = "SoftEdge(profile=LinearRamp(0.05))",
         model = FringeModel(dipole = SoftEdge(profile = ramp))),
    )
end

function _ramp_record()
    profile = LinearRamp(RAMP_WIDTH)
    _, W, shift = profile_function(profile)
    soft = SoftEdge(profile = profile)
    return (
        width_m = RAMP_WIDTH,
        assumption = true,
        half_width_m = W,
        profile_shift_m = shift,
        fringe_integral_m = fringe_integral(profile),
        method = string(soft.method),
    )
end

function _dipole_card(lattice)
    bm = first(e for e in lattice.elements if e isa SBend && e.name == "BM")
    return (
        name = bm.name,
        count = count(e -> e isa SBend && e.name == "BM", lattice.elements),
        L_m = bm.L,
        angle_rad = bm.angle,
        e1 = bm.e1,
        e2 = bm.e2,
        gap_m = bm.gap,
        fint = bm.fint,
        k1 = bm.k1,
    )
end

# --- one model ---------------------------------------------------------------

function _read_model(lattice, spec, body, ramp)
    ref = lattice.ref
    C = path_length(lattice)
    mag = magnetic_lattice(lattice)
    length(mag.elements) == length(lattice.elements) ||
        error("$(spec.key): magnetic_lattice dropped elements")
    any(e -> e isa RFCavity, mag.elements) &&
        error("$(spec.key): readings include a cavity")
    println("ring_lie_factor  ", spec.key)
    flush(stdout)
    t0 = time()
    F = ring_lie_factor(mag, spec.model)
    M = factor_to_map(F)
    pe = path_expansion(M, ref, C)
    po = periodic_optics(F.R, ref, C)
    po.x.stable && po.y.stable || error("$(spec.key): periodic optics is unstable")
    # A is built inside generator_chromaticity from this factor's optics.
    gen = generator_chromaticity(F, po)
    ch = chromaticity_orders(M)
    xi_x = XI_FACTOR * gen.x1.c
    xi_y = XI_FACTOR * gen.y1.c
    dQx = _tune_shift(po.x.mu, body.po.x.mu)
    dQy = _tune_shift(po.y.mu, body.po.y.mu)
    if max(abs(dQx), abs(dQy)) > TUNE_SHIFT_MAX
        stop!("$(spec.key) fringe tune shift $(dQx), $(dQy) exceeds $TUNE_SHIFT_MAX; the integer part of Q is ambiguous")
    end
    Qx = body.along.Qx + dQx
    Qy = body.along.Qy + dQy
    slip = eta1_formula(pe.alpha0, pe.alpha1, ref.gamma)
    if abs(pe.eta1 - slip) > SLIP_TOL
        stop!(string(
            spec.key, " path_expansion eta1 ", pe.eta1,
            " and the slip identity ", slip, " differ by ", pe.eta1 - slip,
        ))
    end
    q_gap_x = _frac_diff(ch.x[1], po.x.mu / 2π)
    q_gap_y = _frac_diff(ch.y[1], po.y.mu / 2π)
    if max(abs(q_gap_x), abs(q_gap_y)) > MATRIX_TOL
        stop!("$(spec.key) fractional tune from chromaticity_orders and μ/2π differ by $q_gap_x, $q_gap_y")
    end
    xi_gap_x = xi_x - ch.x[2]
    xi_gap_y = xi_y - ch.y[2]
    if max(abs(xi_gap_x), abs(xi_gap_y)) > XI_READ_TOL
        stop!("$(spec.key) generator chromaticity and chromaticity_orders differ by $xi_gap_x, $xi_gap_y")
    end
    if !_norm_ok(gen)
        stop!("$(spec.key) normalising matrix of this factor's optics is outside the step-2b tolerance")
    end
    dR = maximum(abs, F.R .- body.R)
    optics_changed = dR > MATRIX_TOL || max(abs(dQx), abs(dQy)) > MATRIX_TOL
    sym = symplectic_error(F.R)
    @printf("  %-8s %7.2f s   Q % .8f % .8f   xi % .6e % .6e   eta1 % .6e   |dR| %.3e\n",
            spec.key, time() - t0, Qx, Qy, xi_x, xi_y, pe.eta1, dR)
    flush(stdout)
    return (
        key = spec.key,
        dipole = spec.dipole,
        model = spec.model,
        Qx = Qx,
        Qy = Qy,
        Qx_frac = po.x.mu / 2π,
        Qy_frac = po.y.mu / 2π,
        dQx = dQx,
        dQy = dQy,
        xi_x = xi_x,
        xi_y = xi_y,
        xi_orders_x = ch.x[2],
        xi_orders_y = ch.y[2],
        xi_gap_x = xi_gap_x,
        xi_gap_y = xi_gap_y,
        eta0 = pe.eta0,
        eta1 = pe.eta1,
        alpha0 = pe.alpha0,
        alpha1 = pe.alpha1,
        slip = slip,
        po = po,
        gen = gen,
        norm_ok = _norm_ok(gen),
        matrix_gap = dR,
        optics_changed = optics_changed,
        symplectic_error = sym,
        seconds = time() - t0,
        ramp = spec.key == "SoftEdge" ? ramp : nothing,
        C = C,
    )
end

# --- tracking ----------------------------------------------------------------

function _launch(po, delta, Jx, Jy)
    ax = Jx > 0 ? sqrt(2 * Jx * po.x.beta) : 0.0
    ay = Jy > 0 ? sqrt(2 * Jy * po.y.beta) : 0.0
    Dx, Dpx, Dy, Dpy = po.D4
    x = ax + Dx * delta
    px = -po.x.alpha / po.x.beta * ax + Dpx * delta
    y = ay + Dy * delta
    py = -po.y.alpha / po.y.beta * ay + Dpy * delta
    return [x, px, y, py, 0.0, delta]
end

function _design(times, Q)
    θ = (2π * Q) .* times
    return hcat(ones(length(times)), cos.(θ), sin.(θ), cos.(2 .* θ), sin.(2 .* θ),
                cos.(3 .* θ), sin.(3 .* θ))
end

function _coef(deltas, times, Q)
    A = _design(times, Q)
    coef = A \ deltas
    rms = sqrt(sum(abs2, deltas .- A * coef) / length(deltas))
    return coef, rms
end

"Minimise the harmonic residual in a ±2% window around Q0, then zoom."
function _refine_Q(deltas, Q0)
    times = float.(0:length(deltas)-1)
    center = Float64(Q0)
    half = 0.02 * center
    bestQ = center
    best = Inf
    for _ in 1:6
        step = 2half / 40
        for i in 0:40
            Q = center - half + i * step
            Q > 0 || continue
            _, rms = _coef(deltas, times, Q)
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

function _fit_delta(deltas, Q0)
    n = length(deltas)
    n >= 16 || error("harmonic fit needs at least 16 samples, got $n")
    times = float.(0:n-1)
    Q = _refine_Q(deltas, Q0)
    coef, rms = _coef(deltas, times, Q)
    a1 = hypot(coef[2], coef[3])
    a2 = hypot(coef[4], coef[5])
    a3 = hypot(coef[6], coef[7])
    n2 = n ÷ 2
    unc = NaN
    if n2 >= 16
        c1, _ = _coef(deltas[1:n2], times[1:n2], Q)
        c2, _ = _coef(deltas[n2+1:end], times[n2+1:end], Q)
        unc = abs(c1[1] - c2[1]) / 2
    end
    return (
        delta_eq = coef[1],
        delta2 = (a1^2 + a2^2 + a3^2) / 2,
        a1 = a1,
        rms = rms,
        unc = unc,
        n = n,
        Q = Q,
    )
end

function _edge_factors(items, ref)
    cache = Dict{Any,LieFactor}()
    factors = Union{Nothing,LieFactor}[nothing for _ in items]
    n_edge = 0
    for (i, it) in enumerate(items)
        it isa EdgeMap || continue
        it isa AbstractElement &&
            error("EdgeMap is an AbstractElement; the tracking split would be wrong")
        n_edge += 1
        factors[i] = get!(cache, (it.parent, it.side, it.model)) do
            dragt_finn_factor(edge_map(it, ref); name = it.name)
        end
    end
    return factors, n_edge
end

"""
One turn of `fringe_sequence`. An `EdgeMap` is the Lie factor of `edge_map`.
A magnet is `exact_track_element`. A thin cavity is `rf_kick`.
"""
function _track_turn(z, items, factors, ref, path_length; ds::Real = TRACK_DS)
    w = collect(Float64, z)
    for (it, F) in zip(items, factors)
        if it isa EdgeMap
            w = apply_factor(F, w)
        elseif it isa RFCavity
            it.L == 0 || error("track expected a thin cavity, L = $(it.L)")
            w = rf_kick(w, it, ref, path_length)
        else
            w = exact_track_element(w, it, ref, path_length; ds = ds)
        end
    end
    return w
end

function _probe_exact_track(lattice, model)
    items = fringe_sequence(lattice, model)
    edge = findfirst(it -> it isa EdgeMap, items)
    edge === nothing && error("HardEdge sequence has no EdgeMap")
    it = items[edge]
    try
        exact_track_element(zeros(6), it, lattice.ref, path_length(lattice); ds = TRACK_DS)
        return (
            accepted = true,
            method_error = false,
            exception = "",
            message = "exact_track_element returned on an EdgeMap",
            edge_name = it.name,
            parent = it.parent.name,
            side = string(it.side),
        )
    catch err
        full = sprint(showerror, err)
        kept = String[]
        for line in split(full, '\n')
            startswith(line, "Closest candidates") && break
            push!(kept, line)
            length(kept) >= 6 && break
        end
        return (
            accepted = false,
            method_error = err isa MethodError,
            exception = string(typeof(err)),
            message = join(kept, "\n"),
            edge_name = it.name,
            parent = it.parent.name,
            side = string(it.side),
        )
    end
end

function _short_track(lattice, reading, regime, nturns)
    label = reading.key
    ref = lattice.ref
    C = reading.C
    println("Fringe track ", label, "  turns ", nturns, "  ds ", TRACK_DS)
    flush(stdout)
    lat, cav = append_thin_cavity(lattice, regime.voltage_MV, regime.harmonic, BELL_PHI_S)
    cav.L == 0 || error("cavity is not thin")
    abs(path_length(lat) - C) < 1.0e-12 || error("the thin cavity changed the path length")
    items = fringe_sequence(lat, reading.model)
    factors, n_edge = _edge_factors(items, ref)
    z = _launch(reading.po, regime.delta_m, regime.Jx, regime.Jy)
    Qs_note = synchrotron_tune_note(ref, reading.eta0, regime.voltage_MV, regime.harmonic)
    # The file Qs belongs to the step-2 slip. A model that moves η0 is fitted
    # from the note tune of that η0, at the same voltage.
    Q0 = abs(Qs_note - regime.Qs) / regime.Qs < 0.02 ? regime.Qs : Qs_note
    deltas = Vector{Float64}(undef, nturns + 1)
    aborted = false
    used = 0
    t0 = time()
    for i in 0:nturns
        if i > 0
            z = _track_turn(z, items, factors, ref, C; ds = TRACK_DS)
        end
        used = i + 1
        if !all(isfinite, z) || abs(z[6]) > 0.01
            aborted = true
            deltas[used] = all(isfinite, z) ? z[6] : NaN
            println("  ", label, " left the window at turn ", i)
            flush(stdout)
            break
        end
        deltas[used] = z[6]
        if i == 0 || i == nturns || i % 100 == 0
            @printf("  %s turn %5d / %d   elapsed %7.1f s   delta % .6e\n",
                    label, i, nturns, time() - t0, z[6])
            flush(stdout)
        end
    end
    if used < 16
        stop!("$label fringe track kept $used samples")
        fit = (delta_eq = NaN, delta2 = NaN, a1 = NaN, rms = NaN, unc = NaN, n = used, Q = Q0)
    else
        fit = _fit_delta(deltas[1:used], Q0)
    end
    return (
        fit = fit,
        aborted = aborted,
        nturns = nturns,
        used = used,
        elapsed = time() - t0,
        n_edge = n_edge,
        Qs_note = Qs_note,
        Q_seed = Q0,
    )
end

function _judge_track(tag, formula, measured, unc)
    gap = abs(measured - formula)
    tol = max(TRACK_REL * abs(formula), TRACK_ABS)
    if !isfinite(measured)
        stop!("$tag fringe track has no finite delta_eq")
    elseif gap > tol
        stop!("$tag fringe delta_eq $measured versus formula $formula, gap $gap, tolerance $tol")
    else
        @printf("  agree  %-10s  |gap| %.3e   tol %.3e   unc %.3e\n", tag, gap, tol, unc)
        flush(stdout)
    end
    if isfinite(unc) && unc > tol
        println("  note   ", tag, " half-window split ", unc,
                " is above the tolerance. delta_eq is the fit on the full window.")
        flush(stdout)
    end
    return gap
end

# --- report ------------------------------------------------------------------

function _fmt(x)
    return x isa AbstractFloat ? (isfinite(x) ? @sprintf("%.8e", float(x)) : "nan") : string(x)
end

function _md_table(headers, rows)
    lines = String["| " * join(headers, " | ") * " |",
                   "| " * join(fill("---", length(headers)), " | ") * " |"]
    for r in rows
        push!(lines, "| " * join(r, " | ") * " |")
    end
    return join(lines, "\n")
end

function _push_num!(rows, model, quantity, value)
    push!(rows, (model = model, quantity = quantity, value = value))
    return rows
end

function _write_numbers(path, setup, body, ramp, dipole, readings, tracks, probe)
    rows = NamedTuple[]
    reg = setup.regime
    _push_num!(rows, "setup", "delta_m", reg.delta_m)
    _push_num!(rows, "setup", "Jx", reg.Jx)
    _push_num!(rows, "setup", "Jy", reg.Jy)
    _push_num!(rows, "setup", "Qs_file", reg.Qs)
    _push_num!(rows, "setup", "voltage_MV", reg.voltage_MV)
    _push_num!(rows, "setup", "ramp_width_m", ramp.width_m)
    _push_num!(rows, "setup", "ramp_half_width_m", ramp.half_width_m)
    _push_num!(rows, "setup", "ramp_fringe_integral_m", ramp.fringe_integral_m)
    _push_num!(rows, "setup", "bm_L_m", dipole.L_m)
    _push_num!(rows, "setup", "bm_gap_m", dipole.gap_m)
    _push_num!(rows, "setup", "bm_fint", dipole.fint)
    _push_num!(rows, "setup", "bm_e1", dipole.e1)
    _push_num!(rows, "setup", "body_Qx", body.along.Qx)
    _push_num!(rows, "setup", "body_Qy", body.along.Qy)
    _push_num!(rows, "setup", "k2_SXD", setup.k2["SXD"])
    _push_num!(rows, "setup", "k2_SXB", setup.k2["SXB"])
    _push_num!(rows, "setup", "k2_SYB", setup.k2["SYB"])
    for (rd, tr) in zip(readings, tracks)
        for (q, v) in (
            ("Qx", rd.Qx), ("Qy", rd.Qy),
            ("Qx_frac", rd.Qx_frac), ("Qy_frac", rd.Qy_frac),
            ("dQx", rd.dQx), ("dQy", rd.dQy),
            ("xi_x", rd.xi_x), ("xi_y", rd.xi_y),
            ("xi_orders_x", rd.xi_orders_x), ("xi_orders_y", rd.xi_orders_y),
            ("eta0", rd.eta0), ("eta1", rd.eta1),
            ("alpha0", rd.alpha0), ("alpha1", rd.alpha1),
            ("Lambda", tr.Lambda),
            ("formula_delta_eq", tr.formula_delta_eq),
            ("formula_gamma_eff", tr.formula_gamma),
            ("formula_dE_over_E", tr.formula_dE),
            ("track_delta_eq", tr.track.fit.delta_eq),
            ("track_delta2", tr.track.fit.delta2),
            ("track_unc_delta_eq", tr.track.fit.unc),
            ("track_gamma_eff", tr.track_gamma),
            ("track_dE_over_E", tr.track_dE),
            ("track_rms", tr.track.fit.rms),
            ("track_Q", tr.track.fit.Q),
            ("Qs_note", tr.track.Qs_note),
            ("matrix_gap", rd.matrix_gap),
            ("symplectic_error", rd.symplectic_error),
            ("rotation_residual", rd.gen.rotation_residual),
            ("beta_x", rd.po.x.beta),
            ("beta_y", rd.po.y.beta),
            ("D_x", rd.po.D4[1]),
            ("n_edge", float(tr.track.n_edge)),
            ("turns", float(tr.track.nturns)),
        )
            _push_num!(rows, rd.key, q, v)
        end
    end
    _push_num!(rows, "probe", "exact_track_accepted_EdgeMap", probe.accepted ? 1.0 : 0.0)
    open(path, "w") do io
        println(io, "model,quantity,value")
        for r in rows
            println(io, r.model, ",", r.quantity, ",", _fmt(r.value))
        end
    end
    return path
end

function _comparison_rows(rows)
    out = Vector{String}[]
    for r in rows
        push!(out, String[
            r.dipole,
            _fmt(r.Qx), _fmt(r.Qy),
            _fmt(r.xi_x), _fmt(r.xi_y),
            _fmt(r.eta0), _fmt(r.eta1),
            _fmt(r.Lambda),
            _fmt(r.formula_delta_eq),
            _fmt(r.track.fit.delta_eq),
            _fmt(r.formula_gamma),
            _fmt(r.track_gamma),
        ])
    end
    return out
end

function _write_report(path, setup, body, ramp, dipole, rows, probe)
    ref = setup.ref
    reg = setup.regime
    lines = String[]
    push!(lines, "# Dipole fringe")
    push!(lines, "")
    push!(lines, "Lattice file: `$(closed_lattice_path())`. Reference: $(ref.particle.name), kinetic energy $(ref.Ekin_MeV) MeV, harmonic $(setup.lattice.harmonic). Placement is `place_sextupole_families` on `optics_along` at ds = $(DS) m of the closed ring, then `set_family_strengths` with `k2_per_m3.own` from `output/05_sext_strengths.json`. Strengths are not solved again.")
    push!(lines, "")
    push!(lines, "Own strengths k2 [1/m³]: SXD = $(_fmt(setup.k2["SXD"])), SXB = $(_fmt(setup.k2["SXB"])), SYB = $(_fmt(setup.k2["SYB"])). Family indices match step 5.")
    push!(lines, "")
    push!(lines, "Each row is `ring_lie_factor(lat, model)` with the fringe model on dipoles. Quadrupoles and sextupoles stay `NoEdge`. `Qx` and `Qy` are the body tunes from `optics_along` plus the shortest fringe shift of μ/2π from `periodic_optics` of that factor. `ξ` is `generator_chromaticity` with the step-2b factor `ξ = −c / (2π)`. The normalising matrix is built from `periodic_optics` of the same factor. `η1` is `path_expansion` of `factor_to_map` of that factor.")
    push!(lines, "")
    push!(lines, "## Assumption")
    push!(lines, "")
    push!(lines, "The soft-edge width is an assumption of this step: `LinearRamp($(RAMP_WIDTH))`, width $(RAMP_WIDTH) m. The closed lattice does not carry this width. BELL's half-width of that profile is W = $(_fmt(ramp.half_width_m)) m (profile shift $(_fmt(ramp.profile_shift_m)) m). The fringe integral ∫ F(1−F) dz = $(_fmt(ramp.fringe_integral_m)) m. The default soft-edge method is `$(ramp.method)`.")
    push!(lines, "")
    push!(lines, "Loaded BM card: $(dipole.count) magnets, L = $(_fmt(dipole.L_m)) m, angle = $(_fmt(dipole.angle_rad)) rad, e1 = $(_fmt(dipole.e1)), e2 = $(_fmt(dipole.e2)), gap = $(_fmt(dipole.gap_m)) m, FINT = $(_fmt(dipole.fint)), k1 = $(_fmt(dipole.k1)). The ramp satisfies 2W < L.")
    push!(lines, "")
    push!(lines, "Body tunes, with no edge model: Qx = $(_fmt(body.along.Qx)), Qy = $(_fmt(body.along.Qy)).")
    push!(lines, "")
    push!(lines, "## Comparison")
    push!(lines, "")
    push!(lines, "Regime file `output/02_rf_regimes.json`, object `linear`. Voltage $(reg.voltage_MV) MV, harmonic $(reg.harmonic), δ_m = $(reg.delta_m), Q_s = $(reg.Qs), φ̂ = $(reg.phi_hat). The notes' stationary phase in that file is $(reg.notes_phi_s). The cavity in the track is at BELL phase π. Actions are the placeholder Jx = $(reg.Jx), Jy = $(reg.Jy), the same values as steps 4 and 5. Λ = −(2π/C) (Jx ξx + Jy ξy) uses the ξ of that fringe model. The linear-regime formula is δ_eq = −(η1 δ_m²/2 + Λ) / η0. η0 is the slip of that same factor. γ_eff = γ + β² γ δ_eq + β² ⟨δ²⟩ / (2γ), with ⟨δ²⟩ = δ_m²/2 in the formula and the harmonic mean square in the track.")
    push!(lines, "")
    push!(lines, _md_table(
        ["model", "Qx", "Qy", "ξx", "ξy", "η0", "η1", "Λ", "formula δ_eq", "track δ_eq", "formula γ_eff", "track γ_eff"],
        _comparison_rows(rows),
    ))
    push!(lines, "")
    push!(lines, "`SoftEdge` uses the assumed width $(RAMP_WIDTH) m. `HardEdge()` is `HardEdge(order=2)` on dipoles, with the magnet gap and FINT (both zero on this card).")
    push!(lines, "")
    push!(lines, "## Each model")
    push!(lines, "")
    for r in rows
        push!(lines, "### $(r.key)")
        push!(lines, "")
        push!(lines, "Dipole model: `$(r.dipole)`. Fractional tunes μ/2π: Qx = $(_fmt(r.Qx_frac)), Qy = $(_fmt(r.Qy_frac)). Fringe shift: ΔQx = $(_fmt(r.dQx)), ΔQy = $(_fmt(r.dQy)). |R − body matrix| = $(_fmt(r.matrix_gap)). Symplectic error of R = $(_fmt(r.symplectic_error)).")
        push!(lines, "")
        if r.optics_changed
            push!(lines, "The linear optics differ from the body, so the normalising matrix was built from this factor: βx = $(_fmt(r.po.x.beta)) m, βy = $(_fmt(r.po.y.beta)) m, Dx = $(_fmt(r.po.D4[1])) m.")
        else
            push!(lines, "The linear map matches the body within $MATRIX_TOL. The normalising matrix was still built from this factor's `periodic_optics`. βx = $(_fmt(r.po.x.beta)) m, βy = $(_fmt(r.po.y.beta)) m, Dx = $(_fmt(r.po.D4[1])) m.")
        end
        push!(lines, "")
        push!(lines, "Generator chromaticity ξx = $(_fmt(r.xi_x)), ξy = $(_fmt(r.xi_y)). `chromaticity_orders` reads ξx = $(_fmt(r.xi_orders_x)), ξy = $(_fmt(r.xi_orders_y)). Gaps: $(_fmt(r.xi_gap_x)), $(_fmt(r.xi_gap_y)). Rotation residual of the normalised linear map: $(_fmt(r.gen.rotation_residual)). Normalisation within the step-2b tolerance: $(r.norm_ok).")
        push!(lines, "")
        push!(lines, "Slip: η0 = $(_fmt(r.eta0)), η1 = $(_fmt(r.eta1)), α0 = $(_fmt(r.alpha0)), α1 = $(_fmt(r.alpha1)). Identity η1 = α1 − α0/γ² + (3γ² − 1)/(2γ⁴) gives $(_fmt(r.slip)).")
        push!(lines, "")
        tr = r.track
        push!(lines, "Track: $(tr.n_edge) edge items, $(tr.nturns) turns, ds = $(TRACK_DS) m, seed Q_s = $(tr.Q_seed), note Q_s from this η0 = $(tr.Qs_note). Fit Q_s = $(_fmt(tr.fit.Q)), rms = $(_fmt(tr.fit.rms)), half-window split = $(_fmt(tr.fit.unc)). Window left early: $(tr.aborted).")
        push!(lines, "")
    end
    push!(lines, "## Tracking")
    push!(lines, "")
    push!(lines, "A turn follows `fringe_sequence`. An `EdgeMap` is `apply_factor` of `dragt_finn_factor(edge_map(...))`, the same Lie factor the one-turn product uses, truncated at map order $(MAP_ORDER). A magnet is `exact_track_element` at ds = $(TRACK_DS) m. The thin cavity is `rf_kick` at phase π. `exact_track_element` is not called on an `EdgeMap`.")
    push!(lines, "")
    if probe.accepted
        push!(lines, "Probe: `exact_track_element` accepted `$(probe.edge_name)` ($(probe.parent), side $(probe.side)). The expected method gap was not observed.")
    else
        push!(lines, "Probe: `exact_track_element` on `$(probe.edge_name)` ($(probe.parent), side $(probe.side)) raises `$(probe.exception)`. The message is recorded in `output/06_bell_finding.md`.")
    end
    push!(lines, "")
    if isempty(STOPS)
        push!(lines, "No stop. Formula and fringe-track δ_eq agree within the step-4 tolerance, generator chromaticity agrees with `chromaticity_orders` within $XI_READ_TOL, and η1 agrees with the slip identity.")
    else
        push!(lines, "Stops ($(length(STOPS))):")
        push!(lines, "")
        for s in STOPS
            push!(lines, "- " * s)
        end
    end
    push!(lines, "")
    write(path, join(lines, "\n"))
    return path
end

function _write_finding(path, probe, dipole)
    lines = String[
        "# BELL finding: exact_track_element has no method for EdgeMap",
        "",
        "Observed while tracking the closed magnetic ring with a dipole fringe model (`studies/magnetic_ring/lattices/magnetic_2_closed.bell`, every BM angle π/8, sextupoles at the step-5 own strengths).",
        "",
        "## What happened",
        "",
        "`fringe_sequence` inserts an `EdgeMap` before and after every dipole whose model is not `NoEdge`. `EdgeMap` is not an `AbstractElement`. `exact_track_element` is defined for `Marker`, `Drift`, `RFCavity`, and `AbstractElement`.",
        "",
        "Calling it on the entrance edge of a BM with `HardEdge()` raises `$(probe.exception)`.",
        "",
        "```",
        probe.message,
        "```",
        "",
        "The edge item in that call is `$(probe.edge_name)`, parent `$(probe.parent)`, side `$(probe.side)`. The loaded BM card has L = $(_fmt(dipole.L_m)) m, e1 = $(_fmt(dipole.e1)), e2 = $(_fmt(dipole.e2)), gap = $(_fmt(dipole.gap_m)) m, FINT = $(_fmt(dipole.fint)).",
        "",
        "The one-turn Lie factor already includes these edges: `ring_lie_factor(lattice, model)` builds a Lie factor of `edge_map` for each `EdgeMap`. There is no matching exact, untruncated tracker.",
        "",
        "## Reproducer",
        "",
        "```julia",
        "using Bell",
        "lat = load_structure(\"studies/magnetic_ring/lattices/magnetic_2_closed.bell\")",
        "model = FringeModel(dipole = HardEdge())",
        "seq = fringe_sequence(lat, model)",
        "edge = first(it for it in seq if it isa EdgeMap)",
        "exact_track_element(zeros(6), edge, lat.ref, path_length(lat))",
        "```",
        "",
        "## What this study did",
        "",
        "The study does not add a method. A turn of `fringe_sequence` applies `dragt_finn_factor(edge_map(edge, ref))` on each `EdgeMap` (`apply_factor`), `exact_track_element` on each magnet, and `rf_kick` on the thin cavity. The edge piece is the Lie factor at map order $(MAP_ORDER). The magnet body is the exact flow.",
        "",
        "## Recommendation",
        "",
        "Give `EdgeMap` an exact tracker that integrates the edge Hamiltonian without a polynomial cut, and accept that item in `exact_track_element` (and in `exact_track_turn`). The Lie factor of `edge_map` would stay the map used by `ring_lie_factor`. Until that method exists, a mixed turn has to branch on `EdgeMap` itself.",
        "",
    ]
    write(path, join(lines, "\n"))
    return path
end

# --- main --------------------------------------------------------------------

function _equilibrium(reading, regime, ref)
    Lambda = betatron_lambda(reading.C, regime.Jx, regime.Jy, reading.xi_x, reading.xi_y)
    de = _delta_eq(reading.eta0, reading.eta1, regime.delta_m, Lambda)
    kin = _gamma_shift(ref.gamma, ref.beta, de, regime.delta_m^2 / 2)
    return (Lambda = Lambda, formula_delta_eq = de, formula_gamma = kin.gamma_eff, formula_dE = kin.dE_over_E)
end

function main()
    _self_check()
    println("Loading ", closed_lattice_path())
    flush(stdout)
    bare = load_closed_lattice()
    _require_bare(bare)
    own = _load_own(STRENGTH_JSON)
    regime = _load_linear(REGIME_JSON)
    regime.notes_phi_s == 0.0 ||
        error("linear regime phase is $(regime.notes_phi_s); expected the notes' 0")
    regime.harmonic == bare.harmonic ||
        error("regime harmonic $(regime.harmonic) differs from the lattice harmonic $(bare.harmonic)")
    placed = _place(bare)
    _check_placement(placed, own)
    lat = set_family_strengths(placed.lattice, placed.fam_idx, own.k2)
    for f in SEXT_FAMILIES
        for i in placed.fam_idx[f]
            lat.elements[i].k2 == own.k2[f] ||
                error("family $f was not set to the own strength")
        end
    end
    ref = lat.ref
    C = path_length(lat)
    R = ring_matrix(lat)
    po = periodic_optics(R, ref, C)
    po.x.stable && po.y.stable || error("periodic optics of the strengthened lattice is unstable")
    along = optics_along(lat, po; ds = DS)
    along.closure.betx < 1.0e-8 || error("strengthened optics_along beta_x did not close")
    body = (R = R, po = po, along = along)
    ramp = _ramp_record()
    dipole = _dipole_card(lat)
    2 * ramp.half_width_m < dipole.L_m ||
        error("assumed fringe width does not fit in BM: 2W = $(2 * ramp.half_width_m), L = $(dipole.L_m)")
    @printf("BM L %.6f m  e1 %.3e  gap %.3e  FINT %.3e  ramp W %.6f m  I %.6e m\n",
            dipole.L_m, dipole.e1, dipole.gap_m, dipole.fint, ramp.half_width_m, ramp.fringe_integral_m)
    flush(stdout)
    @printf("Own k2  SXD %.8e  SXB %.8e  SYB %.8e\n", own.k2["SXD"], own.k2["SXB"], own.k2["SYB"])
    flush(stdout)

    println("Probe exact_track_element on EdgeMap")
    flush(stdout)
    probe = _probe_exact_track(lat, FringeModel(dipole = HardEdge()))
    println("  accepted = ", probe.accepted, "  ", probe.exception)
    flush(stdout)
    if probe.accepted
        stop!("exact_track_element accepted an EdgeMap; the tracking split is unnecessary")
    elseif !probe.method_error
        stop!("exact_track_element on EdgeMap raised $(probe.exception), not MethodError")
    end

    readings = NamedTuple[]
    for spec in _models()
        push!(readings, _read_model(lat, spec, body, ramp))
    end
    # NoEdge is the plain ring. A linear mismatch here is a broken factor, not a fringe effect.
    noedge = readings[1]
    noedge.key == "NoEdge" || error("first model is $(noedge.key)")
    if noedge.matrix_gap > MATRIX_TOL
        stop!("NoEdge one-turn matrix differs from ring_matrix by $(noedge.matrix_gap)")
    end

    n_file = _turns(regime.Qs)
    println("Short track length from the regime file: ", n_file, " turns")
    flush(stdout)
    rows = NamedTuple[]
    for rd in readings
        nturns = n_file
        if !_turns_fixed()
            Qs_note = synchrotron_tune_note(ref, rd.eta0, regime.voltage_MV, regime.harmonic)
            nturns = max(n_file, max(64, round(Int, 1 / Qs_note)))
        end
        eq = _equilibrium(rd, regime, ref)
        @printf("  formula %-8s  delta_eq %.6e   Lambda %.6e   eta0 %.6e   eta1 %.6e\n",
                rd.key, eq.formula_delta_eq, eq.Lambda, rd.eta0, rd.eta1)
        flush(stdout)
        tracked = _short_track(lat, rd, regime, nturns)
        tracked.aborted && stop!("$(rd.key) fringe track left the momentum window")
        kin = _gamma_shift(ref.gamma, ref.beta, tracked.fit.delta_eq, tracked.fit.delta2)
        if !tracked.aborted && isfinite(tracked.fit.delta_eq)
            _judge_track(rd.key, eq.formula_delta_eq, tracked.fit.delta_eq, tracked.fit.unc)
        end
        # Below transition the BELL phase π is the stable stationary phase.
        if rd.eta0 * cos(BELL_PHI_S) <= 0
            stop!("$(rd.key) eta0 $(rd.eta0) is not stable at BELL phase π")
        end
        push!(rows, merge(rd, eq, (
            track = tracked,
            track_gamma = kin.gamma_eff,
            track_dE = kin.dE_over_E,
        )))
    end

    mkpath(OUT_DIR)
    paths = (
        report = joinpath(OUT_DIR, "06_fringe.md"),
        numbers = joinpath(OUT_DIR, "06_numbers.csv"),
        finding = joinpath(OUT_DIR, "06_bell_finding.md"),
    )
    setup = (lattice = lat, ref = ref, regime = regime, k2 = own.k2)
    _write_report(paths.report, setup, body, ramp, dipole, rows, probe)
    _write_numbers(paths.numbers, setup, body, ramp, dipole, readings, rows, probe)
    _write_finding(paths.finding, probe, dipole)
    println("Wrote ", paths.report)
    println("Wrote ", paths.numbers)
    println("Wrote ", paths.finding)
    println(isempty(STOPS) ? "Done." : "Done with $(length(STOPS)) stop(s).")
    flush(stdout)
    return nothing
end

main()
