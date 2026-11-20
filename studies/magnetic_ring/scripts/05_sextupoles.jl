# Step 5. Sextupole families on the closed magnetic ring.
#
#   julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/05_sextupoles.jl
#
# Place SXD, SXB, and SYB with place_sextupole_families. The in-study
# response is a central difference of generator chromaticity and of η1 from
# path_expansion of the magnetic factor. Its solution targets ξx = ξy = η1 = 0.
# solve_A is the neighbour: ξx = ξy = α1 = 0 from ring_metrics. Both settings,
# and k2 = 0, are read in the linear regime of output/02_rf_regimes.json.
# The cavity phase stored in BELL is π. --turns N replaces the short track.

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
const DS = 0.02
const DK = 1.0
const TRACK_DS = 0.1
const BELL_PHI_S = Float64(π)
# Step 2b: ξ1 = −c / (2π), c the coefficient of J δ.
const XI_FACTOR = -1 / (2π)
# Bare |ξ| and |η1| are O(1). Above this, the re-evaluated linear solve
# has not reached the target and NelderMead starts from that solve.
const RESIDUAL_ABS = 1.0e-3
const SLIP_TOL = 1.0e-8
const LINEAR_TOL = 1.0e-8
# Step 4 accepted formula versus exact track within half a percent.
# The floor sits above the residual of that step and below the natural δ_eq.
const TRACK_REL = 0.005
const TRACK_ABS = 1.0e-11
# Generator ξ and chromaticity_orders (default step 1e-4) already differ
# in the sixth digit on this ring. A gap above this is a new disagreement.
const XI_READ_TOL = 1.0e-3
# 0 selects one synchrotron period. A positive value is the default length.
# The k2 = 0 column still has to agree with the formula.
const SHORT_TURNS = 0

const STOPS = String[]

function stop!(msg)
    push!(STOPS, msg)
    println("STOP: ", msg)
    flush(stdout)
    return nothing
end

# --- regime file -------------------------------------------------------------

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

function _json_string(s)
    return "\"" * replace(replace(replace(string(s), "\\" => "\\\\"), "\"" => "\\\""), "\n" => "\\n") * "\""
end

function _json(io, x::Bool, ::Int)
    print(io, x ? "true" : "false")
end
function _json(io, x::Integer, ::Int)
    print(io, x)
end
function _json(io, x::AbstractFloat, ::Int)
    isfinite(x) || error("refusing a non-finite number: $x")
    @printf(io, "%.16g", float(x))
end
function _json(io, ::Nothing, ::Int)
    print(io, "null")
end
function _json(io, x::AbstractString, ::Int)
    print(io, _json_string(x))
end
function _json(io, x::Symbol, ::Int)
    print(io, _json_string(string(x)))
end
function _json(io, x::NamedTuple, level::Int)
    ks = collect(keys(x))
    if isempty(ks)
        print(io, "{}")
        return
    end
    println(io, "{")
    pad = "  "^(level + 1)
    for (i, k) in enumerate(ks)
        print(io, pad, _json_string(string(k)), ": ")
        _json(io, getfield(x, k), level + 1)
        i == length(ks) || print(io, ",")
        println(io)
    end
    print(io, "  "^level, "}")
end
function _json(io, x::AbstractVector, level::Int)
    if isempty(x) || !(first(x) isa NamedTuple)
        print(io, "[")
        for (i, v) in enumerate(x)
            i > 1 && print(io, ", ")
            _json(io, v, level)
        end
        print(io, "]")
        return
    end
    println(io, "[")
    pad = "  "^(level + 1)
    for (i, v) in enumerate(x)
        print(io, pad)
        _json(io, v, level + 1)
        i == length(x) || print(io, ",")
        println(io)
    end
    print(io, "  "^level, "]")
end

function _write_json(path, doc)
    open(path, "w") do io
        _json(io, doc, 0)
        println(io)
    end
    return path
end

# --- small physics -----------------------------------------------------------

function _self_check()
    gamma = 1.5
    alpha0 = 0.2
    target = _target_alpha1(alpha0, gamma)
    abs(eta1_formula(alpha0, target, gamma)) < 1.0e-14 ||
        error("target alpha1 does not null the slip identity")
    abs(XI_FACTOR + 1 / (2π)) < 1.0e-15 || error("chromaticity factor is not -1/(2π)")
    return nothing
end

"""α1 that makes η1 = α1 − α0/γ² + (3γ² − 1)/(2γ⁴) vanish."""
function _target_alpha1(alpha0, gamma)
    g2 = float(gamma)^2
    g2 > 0 || error("gamma^2 is not positive")
    return float(alpha0) / g2 - (3 * g2 - 1) / (2 * g2^2)
end

function _strengths(k)
    return Dict{String,Float64}(f => Float64(get(k, f, 0.0)) for f in SEXT_FAMILIES)
end

function _k_tuple(k)
    d = _strengths(k)
    return (SXD = d["SXD"], SXB = d["SXB"], SYB = d["SYB"])
end

function _residual(st)
    return (xi_x = st.xi_x, xi_y = st.xi_y, eta1 = st.eta1,
            maxabs = maximum(abs, (st.xi_x, st.xi_y, st.eta1)))
end

function _gamma_shift(gamma, beta, delta_eq, delta2)
    b2 = float(beta)^2
    g = float(gamma)
    g_eff = g + b2 * g * delta_eq + b2 * delta2 / (2 * g)
    dE = b2 * (delta_eq + delta2 / (2 * g^2))
    return (gamma_eff = g_eff, dE_over_E = dE)
end

function _delta_eq(eta0, eta1, delta_m, Lambda)
    eta0 == 0 && error("eta0 is zero")
    return -(eta1 * delta_m^2 / 2 + Lambda) / eta0
end

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

function _turns(Qs)
    i = findfirst(==("--turns"), ARGS)
    if i !== nothing
        i < length(ARGS) || error("--turns needs an integer")
        n = parse(Int, ARGS[i + 1])
        n >= 16 || error("--turns must be at least 16, got $n")
        return n
    end
    SHORT_TURNS > 0 && return SHORT_TURNS
    Qs > 0 || error("synchrotron tune is not positive")
    return max(64, round(Int, 1 / Qs))
end

# --- Lie readings ------------------------------------------------------------

function _own_state(lattice; label = "")
    ref = lattice.ref
    C = path_length(lattice)
    mag = magnetic_lattice(lattice)
    length(mag.elements) == length(lattice.elements) ||
        error("$label: magnetic_lattice dropped elements")
    any(e -> e isa RFCavity, mag.elements) &&
        error("$label: magnetic factor was asked for a lattice that still has a cavity")
    t0 = time()
    println("  ring_lie_factor  ", label)
    flush(stdout)
    F = ring_lie_factor(mag)
    M = factor_to_map(F)
    pe = path_expansion(M, ref, C)
    po = periodic_optics(F.R, ref, C)
    po.x.stable && po.y.stable || error("$label: periodic optics is unstable")
    gen = generator_chromaticity(F, po)
    xi_x = XI_FACTOR * gen.x1.c
    xi_y = XI_FACTOR * gen.y1.c
    slip = eta1_formula(pe.alpha0, pe.alpha1, ref.gamma)
    if abs(pe.eta1 - slip) > SLIP_TOL
        stop!(string(
            label, " path_expansion eta1 ", pe.eta1,
            " and alpha formula ", slip, " differ by ", pe.eta1 - slip,
        ))
    end
    @printf("  state %-10s %6.2f s   xi % .6e % .6e   eta1 % .6e   alpha1 % .6e\n",
            label, time() - t0, xi_x, xi_y, pe.eta1, pe.alpha1)
    flush(stdout)
    return (
        xi_x = xi_x, xi_y = xi_y,
        eta0 = pe.eta0, eta1 = pe.eta1,
        alpha0 = pe.alpha0, alpha1 = pe.alpha1,
        eta1_from_alpha = slip,
        C = C,
    )
end

function _own_response(lattice, fam_idx)
    families = [f for f in SEXT_FAMILIES if !isempty(fam_idx[f])]
    println("Own response, dk = ", DK)
    flush(stdout)
    base = _own_state(set_family_strengths(lattice, fam_idx, _strengths(Dict{String,Float64}()));
                      label = "k2=0")
    cols = Dict{String,NamedTuple}()
    for f in families
        plus = _own_state(set_family_strengths(lattice, fam_idx, _strengths(Dict(f => DK)));
                          label = f * "+")
        minus = _own_state(set_family_strengths(lattice, fam_idx, _strengths(Dict(f => -DK)));
                           label = f * "-")
        cols[f] = (
            xi_x = (plus.xi_x - minus.xi_x) / (2 * DK),
            xi_y = (plus.xi_y - minus.xi_y) / (2 * DK),
            eta1 = (plus.eta1 - minus.eta1) / (2 * DK),
        )
    end
    A = [getfield(cols[f], key) for key in (:xi_x, :xi_y, :eta1), f in families]
    return (base = base, families = families, cols = cols, A = A)
end

function _solve_own(base, families, A)
    b = -[base.xi_x, base.xi_y, base.eta1]
    values = A \ b
    return Dict{String,Float64}(f => values[i] for (i, f) in enumerate(families))
end

function _maybe_nelder(lattice, fam_idx, families, k_linear, achieved)
    res = _residual(achieved)
    missing = length(families) < length(SEXT_FAMILIES)
    if !missing && res.maxabs <= RESIDUAL_ABS
        println("NelderMead not used. Three families, re-evaluated residual ", res.maxabs)
        flush(stdout)
        return (used = false, converged = false, n_eval = 0, k = k_linear,
                state = achieved, objective = res.maxabs^2, message = "not needed")
    end
    println("NelderMead from the linear solve. Missing family = ", missing,
            ". Re-evaluated residual = ", res.maxabs)
    flush(stdout)
    x0 = [get(k_linear, f, 0.0) for f in families]
    function objective(x)
        st = _own_state(set_family_strengths(lattice, fam_idx, _strengths(Dict(families[i] => x[i] for i in eachindex(families))));
                        label = "NM")
        val = st.xi_x^2 + st.xi_y^2 + st.eta1^2
        isfinite(val) || error("NelderMead objective is not finite")
        return val
    end
    result = optimize(NelderMead(tol = 1.0e-12, max_iter = 25, initial_step = 0.05), objective, x0)
    k_nm = Dict{String,Float64}(families[i] => result.x[i] for i in eachindex(families))
    state = _own_state(set_family_strengths(lattice, fam_idx, _strengths(k_nm)); label = "NM best")
    use_nm = _residual(state).maxabs < res.maxabs
    @printf("  NelderMead f = %.3e  converged = %s  n_eval = %d  %s\n",
            result.f, result.converged, result.n_eval, result.message)
    println(use_nm ? "  keeping the NelderMead point" : "  keeping the linear solve")
    flush(stdout)
    return (
        used = true,
        converged = result.converged,
        n_eval = result.n_eval,
        k = use_nm ? k_nm : k_linear,
        state = use_nm ? state : achieved,
        objective = use_nm ? result.f : res.maxabs^2,
        message = result.message,
    )
end

# --- placement ---------------------------------------------------------------

function _nearest_quad(opt, s)
    best = nothing
    bd = Inf
    for (s0, s1, e) in opt.bounds
        e isa Quadrupole || continue
        mid = (s0 + s1) / 2
        dist = abs(mid - s)
        if dist < bd
            bd = dist
            best = (name = e.name, k1 = e.k1, s = mid, dist = dist)
        end
    end
    best === nothing && error("no quadrupole in the optics table")
    return best
end

function _global_peak(opt, key)
    arr = getfield(opt.data, key)
    k = argmax(arr)
    d = opt.data
    return (s = d.s[k], value = arr[k], betx = d.betx[k], bety = d.bety[k], D = d.D[k])
end

function _peak_kept(maxima, quantity, s)
    return any(m -> m.quantity == quantity && abs(m.s - s) <= 0.3, maxima)
end

function _picture(rows, notes, maxima, quads)
    lines = String[]
    push!(lines, "The placement picture is SXD at a dispersion maximum, SXB near a horizontal beta maximum, and SYB near a vertical beta maximum.")
    labels = (
        ("SXD", "D_x", "dispersion"),
        ("SXB", "beta_x", "horizontal beta"),
        ("SYB", "beta_y", "vertical beta"),
    )
    for (fam, qty, noun) in labels
        mine = filter(r -> r.family == fam, rows)
        peaks = filter(m -> m.quantity == qty, maxima)
        if isempty(mine)
            push!(lines, fam * " has no magnet, so it is not on a " * noun * " maximum.")
            continue
        end
        all(r -> r.quantity == qty, mine) ||
            push!(lines, fam * " has a row whose quantity is not " * qty * ".")
        peak = maximum(m.value for m in peaks)
        got = maximum(r.value for r in mine)
        shared = any(occursin(fam, n) && occursin("shares slot", n) for n in notes)
        # Notes name the family that already owns the slot, not the family that was skipped.
        skipped = any(occursin("shares slot", n) && occursin(qty == "D_x" ? "D_x" : qty, n) for n in notes)
        quad_names = join(unique(q.name for q in quads[findall(r -> r.family == fam, rows)]), ", ")
        focus = all(q.k1 > 0 for q in quads[findall(r -> r.family == fam, rows)]) ? "focusing" :
                all(q.k1 < 0 for q in quads[findall(r -> r.family == fam, rows)]) ? "defocusing" :
                "mixed focusing and defocusing"
        same = abs(got - peak) <= 1.0e-6 * max(peak, 1.0)
        where = same ? "at the kept $(noun) maximum" : "below the highest kept $(noun) maximum"
        extra = ""
        if fam == "SXB" && any(occursin("beta_x", n) && occursin("shares slot", n) for n in notes)
            extra = " The higher horizontal peaks share their slot with SXD, so this family takes the remaining horizontal peaks."
        end
        if fam == "SYB" && skipped
            extra = " A vertical peak shared a slot and was not given a magnet."
        end
        push!(lines, string(
            fam, " has ", length(mine), " magnets, ", where,
            " (placed ", @sprintf("%.6g", got), ", highest kept ", @sprintf("%.6g", peak),
            "). Nearest quadrupoles are ", focus, " (", quad_names, ").", extra,
        ))
        shared && fam == "SXD" &&
            push!(lines, "SXD keeps the slot that a horizontal beta peak also asked for.")
    end
    return lines
end

function _place(lattice, opt)
    println("place_sextupole_families")
    flush(stdout)
    placed, fam_idx, rows, notes = place_sextupole_families(lattice, opt)
    abs(path_length(placed) - path_length(lattice)) < 1.0e-12 ||
        error("placement changed the path length")
    maxima = find_optics_maxima(opt)
    quads = [_nearest_quad(opt, r.s_slot) for r in rows]
    for f in SEXT_FAMILIES
        n = count(r -> r.family == f, rows)
        n == length(fam_idx[f]) ||
            error("family $f has $n rows and $(length(fam_idx[f])) magnets")
        for i in fam_idx[f]
            e = placed.elements[i]
            (e isa Sextupole && e.name == f && e.k2 == 0) ||
                error("element $i is not a zero-strength $f")
        end
    end
    orphans = [e.name for e in placed.elements if e isa Sextupole && !(e.name in SEXT_FAMILIES)]
    isempty(orphans) || println("Sextupoles left outside the three families: ", join(orphans, ", "))
    peaks = (
        betx = _global_peak(opt, :betx),
        bety = _global_peak(opt, :bety),
        D = _global_peak(opt, :D),
    )
    return (
        lattice = placed,
        fam_idx = fam_idx,
        rows = rows,
        notes = notes,
        maxima = maxima,
        quads = quads,
        orphans = orphans,
        peaks = peaks,
        betx_peak_kept = _peak_kept(maxima, "beta_x", peaks.betx.s),
        bety_peak_kept = _peak_kept(maxima, "beta_y", peaks.bety.s),
        picture = _picture(rows, notes, maxima, quads),
    )
end

function _print_placement(placed)
    @printf("  %-4s %-7s %10s %-8s %10s %10s %10s %10s  nearest quad\n",
            "fam", "quantity", "s_max", "slot", "s_slot", "beta_x", "beta_y", "D")
    for (r, q) in zip(placed.rows, placed.quads)
        kind = q.k1 > 0 ? "F" : "D"
        @printf("  %-4s %-7s %10.4f %-8s %10.4f %10.4f %10.4f %10.4f  %s %s k1=% .4f  %.3f m\n",
                r.family, r.quantity, r.s_max, r.slot, r.s_slot, r.betx, r.bety, r.D,
                q.name, kind, q.k1, q.dist)
    end
    if isempty(placed.notes)
        println("Notes: none. Every kept maximum found a free slot.")
    else
        for note in placed.notes
            println("Note: ", note)
        end
    end
    for f in SEXT_FAMILIES
        n = length(placed.fam_idx[f])
        L = n == 0 ? 0.0 : placed.lattice.elements[placed.fam_idx[f][1]].L
        println("Family ", f, ": ", n, " magnets, L = ", L, " m")
    end
    for line in placed.picture
        println(line)
    end
    gx = placed.peaks.betx
    @printf("Largest beta_x sample: s = %.4f m, beta_x = %.4f m, D = %.4f m, kept by find_optics_maxima: %s\n",
            gx.s, gx.value, gx.D, placed.betx_peak_kept)
    flush(stdout)
    return nothing
end

# --- equilibrium track -------------------------------------------------------

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
    # A fixed Q leaks the fundamental into the constant when the window is
    # not an integer number of true periods. Q is refined; δ_eq is not fitted to a target.
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

function _short_track(lattice, po, regime, ref, C, nturns; label)
    println("Exact track ", label, "  turns ", nturns, "  ds ", TRACK_DS)
    flush(stdout)
    lat, cav = append_thin_cavity(lattice, regime.voltage_MV, regime.harmonic, BELL_PHI_S)
    cav.L == 0 || error("cavity is not thin")
    abs(path_length(lat) - C) < 1.0e-12 || error("the thin cavity changed the path length")
    z = _launch(po, regime.delta_m, regime.Jx, regime.Jy)
    deltas = Vector{Float64}(undef, nturns + 1)
    aborted = false
    used = 0
    t0 = time()
    for i in 0:nturns
        if i > 0
            z = track_turn_exact(z, lat, ref, C; ds = TRACK_DS)
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
        stop!("$label exact track kept $used samples")
        fit = (delta_eq = NaN, delta2 = NaN, a1 = NaN, rms = NaN, unc = NaN, n = used, Q = regime.Qs)
    else
        fit = _fit_delta(deltas[1:used], regime.Qs)
    end
    return (fit = fit, aborted = aborted, nturns = nturns, elapsed = time() - t0)
end

function _judge_track(tag, formula, measured, unc)
    gap = abs(measured - formula)
    tol = max(TRACK_REL * abs(formula), TRACK_ABS)
    if gap > tol
        stop!("$tag exact delta_eq $measured versus formula $formula, gap $gap, tolerance $tol")
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
    return x isa AbstractFloat ? @sprintf("%.8e", float(x)) : string(x)
end

function _md_table(headers, rows)
    lines = String["| " * join(headers, " | ") * " |",
                   "| " * join(fill("---", length(headers)), " | ") * " |"]
    for r in rows
        push!(lines, "| " * join(r, " | ") * " |")
    end
    return join(lines, "\n")
end

function _slot_rows(placed)
    out = NamedTuple[]
    for (r, q) in zip(placed.rows, placed.quads)
        push!(out, (
            family = r.family,
            quantity = r.quantity,
            s_max_m = r.s_max,
            slot = r.slot,
            s_slot_m = r.s_slot,
            beta_x_m = r.betx,
            beta_y_m = r.bety,
            D_m = r.D,
            nearest_quad = q.name,
            quad_k1 = q.k1,
        ))
    end
    return out
end

function _family_block(placed)
    blocks = NamedTuple[]
    for f in SEXT_FAMILIES
        idx = placed.fam_idx[f]
        L = isempty(idx) ? 0.0 : placed.lattice.elements[idx[1]].L
        slots = [r.slot for r in placed.rows if r.family == f]
        push!(blocks, (name = f, count = length(idx), length_m = L,
                       element_index = idx, slots = slots))
    end
    return blocks
end

function _push_num!(rows, setting, source, quantity, value)
    push!(rows, (setting = setting, source = source, quantity = quantity, value = value))
    return rows
end

function _write_numbers(path, own, neighbour, cases, target, condA)
    rows = NamedTuple[]
    _push_num!(rows, "matrix", "own", "condition_number", condA)
    _push_num!(rows, "matrix", "own", "target_alpha1", target)
    for (i, f) in enumerate(own.families)
        for (j, q) in enumerate(("dxi_x", "dxi_y", "deta1"))
            _push_num!(rows, "response", f, q, own.A[j, i])
        end
    end
    for f in SEXT_FAMILIES
        haskey(neighbour.cols, f) || continue
        c = neighbour.cols[f]
        for q in (:xi1x, :xi1y, :alpha1, :eta1)
            _push_num!(rows, "response_bell", f, string(q), c[q])
        end
    end
    for case in cases
        st = case.state
        for (q, v) in (("xi_x", st.xi_x), ("xi_y", st.xi_y), ("eta0", st.eta0),
                       ("eta1", st.eta1), ("alpha0", st.alpha0), ("alpha1", st.alpha1),
                       ("Lambda", case.Lambda),
                       ("formula_delta_eq", case.formula_delta_eq),
                       ("formula_gamma_eff", case.formula_gamma),
                       ("formula_dE_over_E", case.formula_dE),
                       ("exact_delta_eq", case.track.fit.delta_eq),
                       ("exact_delta2", case.track.fit.delta2),
                       ("exact_unc_delta_eq", case.track.fit.unc),
                       ("exact_gamma_eff", case.exact_gamma),
                       ("exact_dE_over_E", case.exact_dE),
                       ("exact_rms", case.track.fit.rms),
                       ("k2_SXD", get(_strengths(case.k), "SXD", 0.0)),
                       ("k2_SXB", get(_strengths(case.k), "SXB", 0.0)),
                       ("k2_SYB", get(_strengths(case.k), "SYB", 0.0)))
            _push_num!(rows, case.name, "study_reading", q, v)
        end
    end
    _push_num!(rows, "solve_A", "bell_metrics", "xi_x", neighbour.metrics.xi1x)
    _push_num!(rows, "solve_A", "bell_metrics", "xi_y", neighbour.metrics.xi1y)
    _push_num!(rows, "solve_A", "bell_metrics", "alpha1", neighbour.metrics.alpha1)
    _push_num!(rows, "solve_A", "bell_metrics", "eta1", neighbour.metrics.eta1)
    _push_num!(rows, "solve_A", "bell_metrics", "formula_delta_eq", neighbour.delta_eq)
    open(path, "w") do io
        println(io, "setting,source,quantity,value")
        for r in rows
            println(io, r.setting, ",", r.source, ",", r.quantity, ",", _fmt(r.value))
        end
    end
    return path
end

function _write_report(path, setup, placed, own, neighbour, cases, target, condA, rankA, nm)
    ref = setup.ref
    lines = String[]
    push!(lines, "# Sextupole families")
    push!(lines, "")
    push!(lines, "Lattice file: `$(closed_lattice_path())`. Reference: $(ref.particle.name), kinetic energy $(ref.Ekin_MeV) MeV, harmonic $(setup.lattice.harmonic). Placement uses `optics_along` at ds = $(DS) m and `place_sextupole_families`. Sextupoles in the closed file are at `k2 = 0`.")
    push!(lines, "")
    push!(lines, "Chromaticity in the in-study matrix is `generator_chromaticity` with the step-2b factor `ξ = −c / (2π)`. `η1` and `α1` are `path_expansion` of `factor_to_map` of `ring_lie_factor` on `magnetic_lattice`. The neighbour cancels `ξx`, `ξy`, and `α1` with `response` and `solve_A`. The finite-difference step is dk = $(DK) for both matrices.")
    push!(lines, "")
    push!(lines, "## Slots")
    push!(lines, "")
    push!(lines, "Each row is one magnet. `beta_x`, `beta_y`, and `D` are the values at the optics maximum that claimed the slot. The slot is the nearest sextupole or short drift.")
    push!(lines, "")
    slot_body = Vector{String}[]
    for (r, q) in zip(placed.rows, placed.quads)
        kind = q.k1 > 0 ? "F" : "D"
        push!(slot_body, [
            r.family, r.quantity, _fmt(r.s_max), r.slot, _fmt(r.s_slot),
            _fmt(r.betx), _fmt(r.bety), _fmt(r.D),
            q.name * " (" * kind * ")",
        ])
    end
    push!(lines, _md_table(
        ["family", "quantity", "s_max [m]", "slot", "s_slot [m]", "beta_x [m]", "beta_y [m]", "D [m]", "nearest quad"],
        slot_body,
    ))
    push!(lines, "")
    if isempty(placed.notes)
        push!(lines, "Notes: none. Every kept maximum found a free slot within 1 m.")
    else
        push!(lines, "Notes from `place_sextupole_families`:")
        push!(lines, "")
        for note in placed.notes
            push!(lines, "- " * note)
        end
    end
    push!(lines, "")
    for f in SEXT_FAMILIES
        n = length(placed.fam_idx[f])
        L = n == 0 ? 0.0 : placed.lattice.elements[placed.fam_idx[f][1]].L
        push!(lines, "Family $(f): $(n) magnets, L = $(L) m.")
    end
    push!(lines, "")
    push!(lines, "## Picture")
    push!(lines, "")
    for line in placed.picture
        push!(lines, line)
        push!(lines, "")
    end
    gx = placed.peaks.betx
    gy = placed.peaks.bety
    gD = placed.peaks.D
    push!(lines, "The largest `beta_x` sample is at s = $(_fmt(gx.s)) m, beta_x = $(_fmt(gx.value)) m, D = $(_fmt(gx.D)) m. `find_optics_maxima` keeps a horizontal peak within 0.3 m of it: $(placed.betx_peak_kept). The largest `beta_y` sample is at s = $(_fmt(gy.s)) m, D = $(_fmt(gy.D)) m, kept: $(placed.bety_peak_kept). The largest dispersion sample is at s = $(_fmt(gD.s)) m, D = $(_fmt(gD.value)) m.")
    push!(lines, "")
    push!(lines, "`find_optics_maxima` keeps a local maximum only where D > 0.5 m, including the beta curves. A horizontal peak that falls in a slot already given to SXD is omitted and printed as a shared-slot note.")
    push!(lines, "")
    push!(lines, "## Own response")
    push!(lines, "")
    push!(lines, "Columns are family strengths in 1/m³. Rows are dξx/dk2, dξy/dk2, and dη1/dk2. Bare readings at k2 = 0: ξx = $(_fmt(own.base.xi_x)), ξy = $(_fmt(own.base.xi_y)), η1 = $(_fmt(own.base.eta1)), α1 = $(_fmt(own.base.alpha1)), α0 = $(_fmt(own.base.alpha0)).")
    push!(lines, "")
    headers = vcat(["row"], own.families)
    body = Vector{String}[]
    for (j, name) in enumerate(("dξx/dk2", "dξy/dk2", "dη1/dk2"))
        push!(body, vcat([name], [_fmt(own.A[j, i]) for i in eachindex(own.families)]))
    end
    push!(lines, _md_table(headers, body))
    push!(lines, "")
    push!(lines, "Condition number cond(A) = $(_fmt(condA)). Rank = $(rankA). Square = $(size(own.A, 1) == size(own.A, 2)).")
    push!(lines, "")
    push!(lines, "The α1 that sets η1 to zero is α0/γ² − (3γ² − 1)/(2γ⁴) = $(_fmt(target)). γ = $(_fmt(ref.gamma)).")
    push!(lines, "")
    st = own.state
    push!(lines, "Achieved after applying the own solution: α1 = $(_fmt(st.alpha1)), η1 = $(_fmt(st.eta1)), ξx = $(_fmt(st.xi_x)), ξy = $(_fmt(st.xi_y)). α1 − target = $(_fmt(st.alpha1 - target)). Residual max(|ξx|, |ξy|, |η1|) = $(_fmt(_residual(st).maxabs)).")
    push!(lines, "")
    push!(lines, "Strengths k2 [1/m³]:")
    push!(lines, "")
    krows = Vector{String}[]
    for f in SEXT_FAMILIES
        push!(krows, [f, _fmt(get(_strengths(own.k_linear), f, 0.0)), _fmt(get(_strengths(own.k), f, 0.0)),
                      _fmt(get(_strengths(neighbour.k), f, 0.0))])
    end
    push!(lines, _md_table(["family", "own linear", "own used", "solve_A"], krows))
    push!(lines, "")
    if nm.used
        push!(lines, "NelderMead ran from the linear solve because a family was missing or the re-evaluated residual exceeded $(RESIDUAL_ABS). Converged = $(nm.converged), evaluations = $(nm.n_eval), $(nm.message). The kept point is the one with the smaller true residual.")
    else
        push!(lines, "NelderMead did not run. All three families are present and the re-evaluated residual is at most $(RESIDUAL_ABS).")
    end
    push!(lines, "")
    push!(lines, "## Neighbour, solve_A")
    push!(lines, "")
    push!(lines, "`solve_A` cancels ξx, ξy, and α1 using `ring_metrics` (`chromaticity_orders` and `path_expansion`). It does not target η1.")
    push!(lines, "")
    bheaders = ["family", "dξx/dk2", "dξy/dk2", "dα1/dk2", "dη1/dk2"]
    bbody = Vector{String}[]
    for f in SEXT_FAMILIES
        haskey(neighbour.cols, f) || continue
        c = neighbour.cols[f]
        push!(bbody, [f, _fmt(c[:xi1x]), _fmt(c[:xi1y]), _fmt(c[:alpha1]), _fmt(c[:eta1])])
    end
    push!(lines, _md_table(bheaders, bbody))
    push!(lines, "")
    m = neighbour.metrics
    push!(lines, "At the `solve_A` strengths, `ring_metrics` reads ξx = $(_fmt(m.xi1x)), ξy = $(_fmt(m.xi1y)), α1 = $(_fmt(m.alpha1)), η1 = $(_fmt(m.eta1)), η0 = $(_fmt(m.eta0)). Leftover η1 is the slip that remains when α1 is removed and the kinematic piece (3γ² − 1)/(2γ⁴) − α0/γ² is left. Formula δ_eq from these BELL readings, with Λ from the same ξ and the regime actions, is $(_fmt(neighbour.delta_eq)).")
    push!(lines, "")
    push!(lines, "The study reading of that same lattice is in the table below. |η1_study − η1_bell| = $(_fmt(abs(neighbour.study.eta1 - m.eta1))). |ξx_study − ξx_bell| = $(_fmt(abs(neighbour.study.xi_x - m.xi1x))).")
    push!(lines, "")
    push!(lines, "## Linear-regime equilibrium")
    push!(lines, "")
    reg = setup.regime
    push!(lines, "Regime file `output/02_rf_regimes.json`, object `linear`. Voltage $(reg.voltage_MV) MV, harmonic $(reg.harmonic), δ_m = $(reg.delta_m), Q_s = $(reg.Qs), φ̂ = $(reg.phi_hat). The notes' stationary phase in that file is $(reg.notes_phi_s). The cavity in the track is at BELL phase π. Actions are the placeholder Jx = $(reg.Jx), Jy = $(reg.Jy). Λ = −(2π/C) (Jx ξx + Jy ξy) uses the ξ of that setting. The formula is δ_eq = −(η1 δ_m²/2 + Λ) / η0. γ_eff = γ + β² γ δ_eq + β² ⟨δ²⟩ / (2γ), with ⟨δ²⟩ = δ_m²/2 in the formula and the harmonic mean square in the track.")
    push!(lines, "")
    push!(lines, "The exact column is `exact_track_element` on the magnets and `rf_kick` on one thin cavity, ds = $(TRACK_DS) m, $(cases[1].track.nturns) turns. δ_eq is the constant term of a harmonic fit. The tune starts at the regime Q_s and is refined on the residual, so a fraction of a turn of phase error does not leak into the constant.")
    push!(lines, "")
    erows = Vector{String}[]
    for case in cases
        push!(erows, [
            case.name,
            _fmt(case.state.xi_x), _fmt(case.state.xi_y),
            _fmt(case.state.eta1), _fmt(case.state.alpha1),
            _fmt(case.Lambda),
            _fmt(case.formula_delta_eq), _fmt(case.track.fit.delta_eq),
            _fmt(case.formula_gamma), _fmt(case.exact_gamma),
        ])
    end
    push!(lines, _md_table(
        ["setting", "ξx", "ξy", "η1", "α1", "Λ", "formula δ_eq", "exact δ_eq", "formula γ_eff", "exact γ_eff"],
        erows,
    ))
    push!(lines, "")
    if isempty(STOPS)
        push!(lines, "No stop. Formula and exact δ_eq agree within the step-4 tolerance, and η1 agrees with the slip identity.")
    else
        push!(lines, "Stops:")
        push!(lines, "")
        for s in STOPS
            push!(lines, "- " * s)
        end
    end
    push!(lines, "")
    push!(lines, "Strengths for the later steps are in `output/05_sext_strengths.json`. The primary vector is `own`: ξx = ξy = η1 = 0. Rebuild by calling `place_sextupole_families` on `optics_along` of the closed lattice at ds = $(DS) m, then `set_family_strengths`.")
    push!(lines, "")
    write(path, join(lines, "\n"))
    return path
end

function _write_strengths(path, placed, own, neighbour, condA, target, nm)
    names = [e.name for e in placed.lattice.elements if e isa Sextupole]
    doc = (
        lattice = "lattices/magnetic_2_closed.bell",
        step = "sextupoles",
        optics_ds_m = DS,
        finite_difference_dk = DK,
        primary = "own",
        primary_means = "k2 that sets xi_x = xi_y = eta1 = 0, from the in-study response. solve_A is the neighbour that sets xi_x = xi_y = alpha1 = 0.",
        rebuild = "Load the closed lattice, run optics_along at optics_ds_m from periodic_optics of ring_matrix, call place_sextupole_families, then set_family_strengths with k2_per_m3.own. element_index is 1-based in that placed lattice.",
        condition_number = condA,
        target_alpha1 = target,
        nelder_mead = nm.used,
        notes = collect(String, placed.notes),
        sextupole_names_in_order = names,
        families = _family_block(placed),
        k2_per_m3 = (
            zero = _k_tuple(Dict{String,Float64}()),
            own_linear = _k_tuple(own.k_linear),
            own = _k_tuple(own.k),
            solve_A = _k_tuple(neighbour.k),
        ),
    )
    return _write_json(path, doc)
end

function _write_finding(path, placed)
    gx = placed.peaks.betx
    lines = String[
        "# BELL finding: dispersion mask and shared sextupole slots",
        "",
        "Observed while placing SXD, SXB, and SYB on the closed magnetic ring (`studies/magnetic_ring/lattices/magnetic_2_closed.bell`, every BM angle π/8, sextupoles at k2 = 0).",
        "",
        "## What happened",
        "",
        "`find_optics_maxima` keeps a local maximum of D_x, beta_x, or beta_y only on samples with D > 0.5 m. The largest beta_x sample on this ring is at s = $(_fmt(gx.s)) m with D = $(_fmt(gx.D)) m, and it is not in that list (`kept = $(placed.betx_peak_kept)`).",
        "",
        "`place_sextupole_families` assigns each kept maximum to the nearest sextupole or short drift. The first family to claim a slot keeps it. Later maxima that fall in that slot are skipped.",
        "",
        "Notes returned on this lattice:",
        "",
    ]
    if isempty(placed.notes)
        push!(lines, "None.")
    else
        for note in placed.notes
            push!(lines, "- " * note)
        end
    end
    push!(lines, "")
    push!(lines, "The two horizontal peaks that sit on the dispersion maxima therefore do not become SXB magnets. SXB is placed on the remaining horizontal peaks. SXD is at the dispersion maxima. SYB is at the vertical maxima.")
    push!(lines, "")
    push!(lines, "## Reproducer")
    push!(lines, "")
    push!(lines, "```julia")
    push!(lines, "lat = load_structure(\"studies/magnetic_ring/lattices/magnetic_2_closed.bell\")")
    push!(lines, "C = path_length(lat)")
    push!(lines, "po = periodic_optics(ring_matrix(lat), lat.ref, C)")
    push!(lines, "opt = optics_along(lat, po; ds = 0.02)")
    push!(lines, "maxima = find_optics_maxima(opt)")
    push!(lines, "placed, fam, rows, notes = place_sextupole_families(lat, opt)")
    push!(lines, "```")
    push!(lines, "")
    push!(lines, "## Recommendation")
    push!(lines, "")
    push!(lines, "Document that the D > dmin mask applies to the beta curves as well as to D_x, and that a maximum whose nearest slot is already taken is omitted rather than moved to the next free slot. The returned note is the right report for a shared slot. No change to the response or the solve is requested from this run: the study uses the families and the notes the functions returned.")
    push!(lines, "")
    write(path, join(lines, "\n"))
    return path
end

# --- main --------------------------------------------------------------------

function _check_linear(tag, base, state)
    for (name, a, b) in (("alpha0", base.alpha0, state.alpha0), ("eta0", base.eta0, state.eta0))
        if abs(a - b) > LINEAR_TOL
            stop!("$tag changed $name from $a to $b")
        end
    end
    return nothing
end

function _matrix_gap(lattice, fam_idx, k)
    R0 = ring_matrix(set_family_strengths(lattice, fam_idx, _strengths(Dict{String,Float64}())))
    R1 = ring_matrix(set_family_strengths(lattice, fam_idx, _strengths(k)))
    return maximum(abs, R0 .- R1)
end

function _equilibrium(name, k, state, setup)
    reg = setup.regime
    Lambda = betatron_lambda(state.C, reg.Jx, reg.Jy, state.xi_x, state.xi_y)
    de = _delta_eq(state.eta0, state.eta1, reg.delta_m, Lambda)
    kin = _gamma_shift(setup.ref.gamma, setup.ref.beta, de, reg.delta_m^2 / 2)
    return (
        name = name, k = k, state = state, Lambda = Lambda,
        formula_delta_eq = de, formula_gamma = kin.gamma_eff, formula_dE = kin.dE_over_E,
    )
end

function main()
    _self_check()
    path = closed_lattice_path()
    println("Loading ", path)
    flush(stdout)
    lattice = load_closed_lattice()
    _require_bare(lattice)
    regime = _load_linear(REGIME_JSON)
    regime.notes_phi_s == 0.0 ||
        error("linear regime phase is $(regime.notes_phi_s); expected the notes' 0")
    regime.harmonic == lattice.harmonic ||
        error("regime harmonic $(regime.harmonic) differs from the lattice harmonic $(lattice.harmonic)")
    ref = lattice.ref
    C = path_length(lattice)
    println("Linear optics, ds = ", DS)
    flush(stdout)
    R = ring_matrix(lattice)
    po = periodic_optics(R, ref, C)
    po.x.stable && po.y.stable || error("periodic optics is unstable")
    opt = optics_along(lattice, po; ds = DS)
    opt.closure.betx < 1.0e-8 || error("optics_along beta_x did not close")
    opt.closure.D < 1.0e-8 || error("optics_along dispersion did not close")

    placed = _place(lattice, opt)
    _print_placement(placed)

    own_resp = _own_response(placed.lattice, placed.fam_idx)
    condA = cond(own_resp.A)
    rankA = rank(own_resp.A)
    @printf("Condition number %.8e   rank %d   size %s\n", condA, rankA, size(own_resp.A))
    flush(stdout)
    for (j, name) in enumerate(("dξx", "dξy", "dη1"))
        println(name, "  ", join([@sprintf("%.6e", own_resp.A[j, i]) for i in eachindex(own_resp.families)], "  "))
    end
    flush(stdout)

    k_linear = _solve_own(own_resp.base, own_resp.families, own_resp.A)
    achieved = _own_state(set_family_strengths(placed.lattice, placed.fam_idx, _strengths(k_linear));
                          label = "own")
    target = _target_alpha1(own_resp.base.alpha0, ref.gamma)
    @printf("Target alpha1 %.8e\n", target)
    @printf("Achieved alpha1 %.8e   eta1 %.8e   xi %.8e %.8e\n",
            achieved.alpha1, achieved.eta1, achieved.xi_x, achieved.xi_y)
    flush(stdout)
    nm = _maybe_nelder(placed.lattice, placed.fam_idx, own_resp.families, k_linear, achieved)
    own = (
        base = own_resp.base,
        families = own_resp.families,
        cols = own_resp.cols,
        A = own_resp.A,
        k_linear = k_linear,
        k = nm.k,
        state = nm.state,
    )
    _check_linear("own solution", own.base, own.state)
    gap_own = _matrix_gap(placed.lattice, placed.fam_idx, own.k)
    @printf("Linear-matrix change at the own solution: %.3e\n", gap_own)
    flush(stdout)
    gap_own > LINEAR_TOL && stop!("own solution changed the one-turn linear matrix by $gap_own")

    println("BELL response and solve_A")
    flush(stdout)
    t_b = time()
    base_b, cols_b = response(placed.lattice, placed.fam_idx; dk = DK)
    @printf("  response in %.2f s\n", time() - t_b)
    flush(stdout)
    if abs(base_b.eta1 - own.base.eta1) > SLIP_TOL
        stop!("k2=0 eta1 from ring_metrics $(base_b.eta1) and path_expansion $(own.base.eta1)")
    end
    if abs(base_b.xi1x - own.base.xi_x) > XI_READ_TOL || abs(base_b.xi1y - own.base.xi_y) > XI_READ_TOL
        stop!("k2=0 generator chromaticity and chromaticity_orders disagree")
    end
    @printf("  k2=0 xi gap generator-orders  %.3e  %.3e\n",
            own.base.xi_x - base_b.xi1x, own.base.xi_y - base_b.xi1y)
    flush(stdout)
    kA = solve_A(base_b, cols_b)
    metrics = ring_metrics(set_family_strengths(placed.lattice, placed.fam_idx, _strengths(kA)))
    study_A = _own_state(set_family_strengths(placed.lattice, placed.fam_idx, _strengths(kA));
                         label = "solve_A")
    if abs(study_A.eta1 - metrics.eta1) > SLIP_TOL
        stop!("solve_A eta1 study $(study_A.eta1) and ring_metrics $(metrics.eta1)")
    end
    if abs(study_A.xi_x - metrics.xi1x) > XI_READ_TOL || abs(study_A.xi_y - metrics.xi1y) > XI_READ_TOL
        stop!("solve_A generator chromaticity and chromaticity_orders disagree")
    end
    _check_linear("solve_A", own.base, study_A)
    gap_A = _matrix_gap(placed.lattice, placed.fam_idx, kA)
    gap_A > LINEAR_TOL && stop!("solve_A changed the one-turn linear matrix by $gap_A")
    Lambda_b = betatron_lambda(C, regime.Jx, regime.Jy, metrics.xi1x, metrics.xi1y)
    delta_b = _delta_eq(metrics.eta0, metrics.eta1, regime.delta_m, Lambda_b)
    @printf("solve_A leftover eta1 %.8e   alpha1 %.8e   xi %.8e %.8e   formula delta_eq %.8e\n",
            metrics.eta1, metrics.alpha1, metrics.xi1x, metrics.xi1y, delta_b)
    flush(stdout)
    neighbour = (
        cols = Dict(f => cols_b[f][1] for f in keys(cols_b)),
        k = kA,
        metrics = metrics,
        study = study_A,
        delta_eq = delta_b,
    )

    ref.gamma > 1 || error("gamma is not above 1")
    own.base.eta0 * cos(BELL_PHI_S) > 0 || error("phi_s = π is not stable for this eta0")
    setup = (lattice = lattice, ref = ref, regime = regime, po = po)
    nturns = _turns(regime.Qs)
    println("Short track length ", nturns, " turns")
    flush(stdout)
    mkpath(OUT_DIR)
    _write_strengths(joinpath(OUT_DIR, "05_sext_strengths.json"), placed, own, neighbour, condA, target, nm)
    println("Wrote strengths before the tracks")
    flush(stdout)
    bare_case = _equilibrium("k2=0", Dict{String,Float64}(), own.base, setup)
    own_case = _equilibrium("own", own.k, own.state, setup)
    bell_case = _equilibrium("solve_A", kA, study_A, setup)
    cases = NamedTuple[]
    for (spec, latk) in (
        (bare_case, set_family_strengths(placed.lattice, placed.fam_idx, _strengths(Dict{String,Float64}()))),
        (own_case, set_family_strengths(placed.lattice, placed.fam_idx, _strengths(own.k))),
        (bell_case, set_family_strengths(placed.lattice, placed.fam_idx, _strengths(kA))),
    )
        tracked = _short_track(latk, po, regime, ref, C, nturns; label = spec.name)
        tracked.aborted && stop!("$(spec.name) exact track left the momentum window")
        kin = _gamma_shift(ref.gamma, ref.beta, tracked.fit.delta_eq, tracked.fit.delta2)
        if !tracked.aborted && isfinite(tracked.fit.delta_eq)
            _judge_track(spec.name, spec.formula_delta_eq, tracked.fit.delta_eq, tracked.fit.unc)
        end
        push!(cases, merge(spec, (track = tracked, exact_gamma = kin.gamma_eff, exact_dE = kin.dE_over_E)))
    end

    mkpath(OUT_DIR)
    paths = (
        report = joinpath(OUT_DIR, "05_sextupoles.md"),
        numbers = joinpath(OUT_DIR, "05_numbers.csv"),
        strengths = joinpath(OUT_DIR, "05_sext_strengths.json"),
        finding = joinpath(OUT_DIR, "05_bell_finding.md"),
    )
    _write_report(paths.report, setup, placed, own, neighbour, cases, target, condA, rankA, nm)
    _write_numbers(paths.numbers, own, neighbour, cases, target, condA)
    _write_strengths(paths.strengths, placed, own, neighbour, condA, target, nm)
    _write_finding(paths.finding, placed)
    println("Wrote ", paths.report)
    println("Wrote ", paths.numbers)
    println("Wrote ", paths.strengths)
    println("Wrote ", paths.finding)
    println(isempty(STOPS) ? "Done." : "Done with $(length(STOPS)) stop(s).")
    flush(stdout)
    return nothing
end

main()
