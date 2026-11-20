# Step 7. Dynamic aperture on the own sextupole set, NoEdge.
#
#   julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/07_dynamic_aperture.jl
#
# Rebuild the step-5 placement and set k2 from output/05_sext_strengths.json,
# object k2_per_m3.own. The fringe model is FringeModel() (NoEdge).
#
# Lie readings of that one-turn factor: tunes, distance of n_x Q_x + n_y Q_y
# to an integer, max |f3|, and the optional 6D normal form of the linear-regime
# cavity (h4 / c_ij). The magnetic factor has a longitudinal shear, so the
# 6D normal form is not the transverse reading. Guignard and the track use the
# magnetic factor only.
#
# Guignard boundary (sequential linearisation, |Tr R| < 2):
#   approximate {Q_x} by p/q with q <= 100, minimising |q {Q_x} - p|;
#   linear orbit z_{n+1} = R z_n from (x0, 0, 0, 0, 0, δ);
#   at each of the q points, the one-turn Jacobian of the polynomial map
#   truncated at degree 2 (primary) and degree 3 (secondary);
#   Rprod is the product of the horizontal 2×2 blocks;
#   the boundary is the first |x0| at which |Tr Rprod| leaves (-2, 2).
#   δ = 0 and δ = ±δ_m from the linear regime.
#
# Tracking is exact_track_turn. A particle survives when x and y at the ring
# start stay finite and strictly inside the aperture cut after every turn.
# The default is 300 turns. LONG_DA=1 or --long uses 2000 turns and prints
# progress with flush. The figure overlays the |Tr| = 2 lines and the survival
# samples. The two readings are not averaged.

import Pkg

const STUDY_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(STUDY_DIR)

using LinearAlgebra
using Printf
using Bell
using Bell.NormalForm: dragt_finn_nf6d, one_turn_factor_rf
using MagneticRing

const OUT_DIR = joinpath(STUDY_DIR, "output")
const REGIME_JSON = joinpath(OUT_DIR, "02_rf_regimes.json")
const STRENGTH_JSON = joinpath(OUT_DIR, "05_sext_strengths.json")
const DS = 0.02
const TRACK_DS = 0.1
const APERTURE = 0.05
const Q_MAX = 100
const SEARCH_M = 2.0
const SHORT_TURNS = 300
const LONG_TURNS = 2000
const BELL_PHI_S = Float64(π)
const MATRIX_TOL = 1.0e-8
const TRACE_TOL = 1.0e-6
const TUNE_TOL = 1.0e-8

const RESONANCES = (
    (3, 0, "3Qx"),
    (4, 0, "4Qx"),
    (0, 3, "3Qy"),
    (0, 4, "4Qy"),
    (2, 1, "2Qx+Qy"),
    (2, -1, "2Qx-Qy"),
    (1, 2, "Qx+2Qy"),
    (1, -2, "Qx-2Qy"),
    (1, 1, "Qx+Qy"),
    (1, -1, "Qx-Qy"),
    (2, 2, "2Qx+2Qy"),
    (2, 0, "2Qx"),
    (0, 2, "2Qy"),
    (3, 1, "3Qx+Qy"),
    (3, -1, "3Qx-Qy"),
    (1, 3, "Qx+3Qy"),
    (1, -3, "Qx-3Qy"),
)

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
    return String(m.captures[1])
end

function _load_linear(path)
    isfile(path) || error("regime file is missing: $path")
    text = read(path, String)
    lin = _object_after(text, "linear")
    return (
        voltage_MV = _json_number(lin, "voltage_MV"),
        harmonic = Int(round(_json_number(lin, "harmonic"))),
        delta_m = _json_number(lin, "delta_m"),
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
    _json_string(text, "primary") == "own" || error("step 5 primary is not own")
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

function _csv_escape(s)
    t = string(s)
    if occursin(',', t) || occursin('"', t) || occursin('\n', t)
        return "\"" * replace(t, "\"" => "\"\"") * "\""
    end
    return t
end

_long_run() = get(ENV, "LONG_DA", "0") == "1" || any(==("--long"), ARGS)

_turns() = _long_run() ? LONG_TURNS : SHORT_TURNS

# --- arithmetic that does not need the lattice -------------------------------

"Best p/q for a fractional tune, q <= qmax. Minimises |q ν - p|, then q."
function _best_pq(nu; qmax::Int = Q_MAX)
    target = mod(float(nu), 1.0)
    best = (p = 0, q = 1, err = target)
    for q in 1:qmax
        p = round(Int, q * target)
        err = abs(q * target - p)
        if err < best.err - 1.0e-15 || (abs(err - best.err) <= 1.0e-15 && q < best.q)
            best = (p = p, q = q, err = err)
        end
    end
    return best
end

"Distance of a real number to the nearest integer, in [0, 1/2]."
_dist_int(x) = abs(float(x) - round(float(x)))

"Distance on the circle R/Z between two tunes."
function _circle_dist(a, b)
    d = abs(mod(float(a) - float(b), 1.0))
    return min(d, 1.0 - d)
end

"Distance of a fractional tune to the nearest multiple of 1/q."
function _family_dist(frac, q::Int)
    f = mod(float(frac), 1.0)
    k = round(Int, f * q)
    return abs(f - k / q)
end

function _mono_name(e)
    bits = String[]
    for (i, p) in enumerate(e)
        p == 0 && continue
        name = VAR_NAMES[i]
        push!(bits, p == 1 ? name : "$(name)^$(p)")
    end
    return isempty(bits) ? "1" : join(bits, " ")
end

function _leading(poly, degree::Int, n::Int)
    rows = Tuple{Float64,Vector{Int}}[]
    for e in exponents_of_order(degree)
        push!(rows, (Float64(coeff(poly, e)), e))
    end
    sort!(rows, by = r -> -abs(r[1]))
    return rows[1:min(n, length(rows))]
end

function _max_map_degree(M, degree::Int)
    m = 0.0
    for i in 1:6, e in exponents_of_order(degree)
        m = max(m, abs(Float64(coeff(M[i], e))))
    end
    return m
end

"""
Horizontal Guignard trace along the linear orbit.

`Rprod` is J_{q-1} … J_0 of the horizontal 2×2 blocks. `tr_block` is the
trace of the horizontal block of the product of the full 6×6 Jacobians.
The 1986 criterion |Tr R| < 2 is the 2×2 trace.
"""
function _guignard_trace(M, R, x0, delta, q)
    z = [float(x0), 0.0, 0.0, 0.0, 0.0, float(delta)]
    Px = Matrix{Float64}(I, 2, 2)
    P6 = Matrix{Float64}(I, 6, 6)
    for _ in 1:q
        J = jacobian_of_map(M, z)
        Px = J[1:2, 1:2] * Px
        P6 = J * P6
        z = R * z
        if !(all(isfinite, Px) && all(isfinite, P6))
            return (tr_x = Inf, tr_block = Inf, finite = false)
        end
    end
    return (tr_x = tr(Px), tr_block = tr(P6[1:2, 1:2]), finite = true)
end

_stable_trace(t) = t.finite && abs(t.tr_x) < 2

function _trace_disagrees(t)
    sx = t.finite && abs(t.tr_x) < 2
    sb = t.finite && abs(t.tr_block) < 2
    return sx != sb && abs(t.tr_x - t.tr_block) > TRACE_TOL
end

function _amplitudes()
    a = [i * 0.001 for i in 0:30]
    b = [0.03 + i * 0.005 for i in 1:34]
    c = [0.25 + (i - 1) * 0.05 for i in 1:35]
    c[end] = SEARCH_M
    amp = vcat(a, b, c)
    issorted(amp, lt = <) || error("Guignard amplitude grid is not strictly increasing")
    abs(amp[end] - SEARCH_M) < 1.0e-15 || error("Guignard search does not end at $(SEARCH_M)")
    return amp
end

"""
First |x0| at which |Tr| leaves (-2, 2), by bisection of the bracket that
contains that crossing. Samples stay in the result; a later stable sample
does not move the boundary.
"""
function _first_boundary(trace_at, amplitudes)
    samples = NamedTuple[]
    stable = Bool[]
    for a in amplitudes
        t = trace_at(a)
        push!(samples, (x = a, t = t))
        push!(stable, _stable_trace(t))
    end
    origin = samples[1]
    if !stable[1]
        return (
            x = origin.x,
            crossed = false,
            origin_stable = false,
            beyond = false,
            reentrant = false,
            x_unstable = origin.x,
            tr_stable = origin.t,
            tr_unstable = origin.t,
            x_largest_stable = origin.x,
            samples = samples,
        )
    end
    k = findfirst(!, stable)
    last = findlast(identity, stable)
    if k === nothing
        return (
            x = samples[end].x,
            crossed = false,
            origin_stable = true,
            beyond = true,
            reentrant = false,
            x_unstable = samples[end].x,
            tr_stable = samples[end].t,
            tr_unstable = samples[end].t,
            x_largest_stable = samples[end].x,
            samples = samples,
        )
    end
    lo = samples[k - 1].x
    hi = samples[k].x
    a = lo
    b = hi
    ta = samples[k - 1].t
    tb = samples[k].t
    for _ in 1:48
        m = (a + b) / 2
        tm = trace_at(m)
        push!(samples, (x = m, t = tm))
        if _stable_trace(tm)
            a = m
            ta = tm
        else
            b = m
            tb = tm
        end
    end
    return (
        x = a,
        crossed = true,
        origin_stable = true,
        beyond = false,
        reentrant = last > k,
        x_unstable = b,
        tr_stable = ta,
        tr_unstable = tb,
        x_largest_stable = samples[last].x,
        samples = samples,
    )
end

function _fate(z, aperture)
    all(isfinite, z) || return "nonfinite"
    abs(z[1]) < aperture || return "aperture_x"
    abs(z[3]) < aperture || return "aperture_y"
    return "inside"
end

function _track_exact(z0, lattice, nturns, aperture; ds = TRACK_DS, record = false, progress = false)
    z = collect(Float64, z0)
    xs = Float64[]
    pxs = Float64[]
    record && (push!(xs, z[1]); push!(pxs, z[2]))
    launch = _fate(z, aperture)
    if launch != "inside"
        return (turns = 0, survived = false, reason = "launch_$(launch)", x = xs, px = pxs)
    end
    reported = false
    for n in 1:nturns
        local w
        try
            w = exact_track_turn(z, lattice; ds = ds)
        catch err
            if !reported
                println("  exact_track_turn exception: ", sprint(showerror, err))
                flush(stdout)
                reported = true
            end
            return (turns = n - 1, survived = false, reason = "exception", x = xs, px = pxs)
        end
        z = w
        record && (push!(xs, z[1]); push!(pxs, z[2]))
        fate = _fate(z, aperture)
        if fate != "inside"
            return (turns = n, survived = false, reason = fate, x = xs, px = pxs)
        end
        if progress && (n % 100 == 0)
            @printf("  turn %d/%d  x %.6e  y %.6e\n", n, nturns, z[1], z[3])
            flush(stdout)
        end
    end
    return (turns = nturns, survived = true, reason = "survived", x = xs, px = pxs)
end

function _self_check()
    rat = _best_pq(0.35)
    (rat.p == 7 && rat.q == 20 && rat.err < 1.0e-14) ||
        error("best p/q of 0.35 is $(rat.p)/$(rat.q), expected 7/20")
    rat13 = _best_pq(6 / 13)
    (rat13.p == 6 && rat13.q == 13 && rat13.err < 1.0e-14) ||
        error("best p/q of 6/13 is $(rat13.p)/$(rat13.q)")
    abs(_dist_int(1.385) - 0.385) < 1.0e-12 || error("integer distance")
    abs(_circle_dist(0.46, 1 / 3) - abs(0.46 - 1 / 3)) < 1.0e-12 || error("circle distance")
    abs(_family_dist(1 / 3, 3)) < 1.0e-14 || error("family distance")
    M = [Z[1] + 0.5 * Z[1]^2, Z[2], Z[3], Z[4], Z[5], Z[6]]
    J = jacobian_of_map(M, [0.2, 0.0, 0.0, 0.0, 0.0, 0.0])
    abs(J[1, 1] - 1.2) < 1.0e-12 || error("jacobian_of_map returned $(J[1, 1])")
    maximum(abs, J[1, 2:6]) < 1.0e-12 || error("jacobian_of_map has a spurious row")
    # |1.5 + 100 x| reaches 2 at x = 0.005. The stable side is just below.
    toy(x) = (tr_x = 1.5 + 100 * x, tr_block = 1.5 + 100 * x, finite = true)
    b = _first_boundary(toy, [i * 0.001 for i in 0:20])
    b.crossed || error("toy boundary did not cross")
    abs(b.x - 0.005) < 1.0e-9 || error("toy boundary $(b.x) is not 0.005")
    b.x < 0.005 || error("toy boundary is not on the stable side")
    _stable_trace((tr_x = 2.0, tr_block = 2.0, finite = true)) && error("|Tr| = 2 is not inside the open interval")
    _trace_disagrees((tr_x = 1.0, tr_block = 3.0, finite = true)) || error("trace disagreement was missed")
    !_trace_disagrees((tr_x = 1.9999999, tr_block = 2.0000001, finite = true)) ||
        error("round-off across 2 was treated as a disagreement")
    return nothing
end

# --- lattice -----------------------------------------------------------------

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

function _place(lattice)
    ref = lattice.ref
    C = path_length(lattice)
    R = ring_matrix(lattice)
    po = periodic_optics(R, ref, C)
    po.x.stable && po.y.stable || error("periodic optics of the closed lattice is unstable")
    opt = optics_along(lattice, po; ds = DS)
    opt.closure.betx < 1.0e-8 || error("optics_along beta_x did not close")
    println("place_sextupole_families")
    flush(stdout)
    placed, fam_idx, rows, notes = place_sextupole_families(lattice, opt)
    return (lattice = placed, fam_idx = fam_idx, rows = rows, notes = notes)
end

function _check_placement(placed, own)
    for f in SEXT_FAMILIES
        got = placed.fam_idx[f]
        exp = own.index[f]
        got == exp || error("family $f indices $got differ from step 5 $exp")
        for i in got
            e = placed.lattice.elements[i]
            (e isa Sextupole && e.name == f) ||
                error("index $i is $(e.name), expected sextupole $f")
        end
    end
    return nothing
end

function _frac_diff(a, b)
    d = float(a) - float(b)
    return d - round(d)
end

# --- launches and survival edges ---------------------------------------------

function _is_portrait(x, y, delta)
    y == 0.0 || return false
    delta == 0.0 || return false
    for a in (0.005, 0.01, 0.02, 0.04)
        abs(x - a) < 1.0e-12 && return true
    end
    return false
end

function _launches(delta_m, long)
    rows = NamedTuple[]
    seen = Set{NTuple{3,Int}}()
    function add(x, y, d)
        key = (round(Int, x * 1.0e8), round(Int, y * 1.0e8), round(Int, d * 1.0e12))
        key in seen && return
        push!(seen, key)
        push!(rows, (x = float(x), y = float(y), delta = float(d), portrait = _is_portrait(x, y, d)))
    end
    step = long ? 0.0025 : 0.005
    n = round(Int, 0.045 / step)
    for d in (0.0, float(delta_m), -float(delta_m))
        for i in -n:n
            add(i * step, 0.0, d)
        end
    end
    for x in (0.008, 0.009, 0.012, -0.008, -0.009, -0.012)
        add(x, 0.0, 0.0)
    end
    for y in (0.01, -0.01, 0.02, -0.02), x in (0.0, 0.01, -0.01, 0.02, -0.02)
        add(x, y, 0.0)
    end
    sort!(rows, by = r -> (r.delta, r.y, r.x))
    return rows
end

function _survival_edge(tracked, y, delta, side::Int)
    pts = [r for r in tracked if r.y == y && r.delta == delta &&
           (r.x == 0.0 || (side > 0 ? r.x > 0.0 : r.x < 0.0))]
    sort!(pts, by = r -> abs(r.x))
    last_ok = 0.0
    first_bad = nothing
    reentered = Float64[]
    lost = false
    for r in pts
        if r.survived
            if !lost
                last_ok = r.x
            else
                push!(reentered, r.x)
            end
        else
            if !lost
                first_bad = r.x
                lost = true
            end
        end
    end
    return (last_surviving = last_ok, first_lost = first_bad, reentered = reentered)
end

# --- report ------------------------------------------------------------------

function _fmt(x)
    x isa Real || return string(x)
    isfinite(x) || return "Inf"
    return @sprintf("%.8e", float(x))
end

function _fmt6(x)
    isfinite(x) || return "Inf"
    return @sprintf("%.6e", float(x))
end

function _fmt_tr(x)
    isfinite(x) || return "Inf"
    return @sprintf("%.16e", float(x))
end

function _md_table(header, body)
    lines = String[]
    push!(lines, "| " * join(header, " | ") * " |")
    push!(lines, "| " * join(["---" for _ in header], " | ") * " |")
    for row in body
        push!(lines, "| " * join(string.(row), " | ") * " |")
    end
    return join(lines, "\n")
end

function _status_word(b)
    b.origin_stable || return "unstable at x0 = 0"
    b.beyond && return "no crossing inside $(SEARCH_M) m"
    b.crossed && return "first crossing"
    return "not classified"
end

function _write_numbers(path, lie, bounds, tracked, nturns)
    rows = Tuple{String,String,String}[]
    push!(rows, ("C_m", _fmt(lie.C), "path_length"))
    push!(rows, ("Qx", _fmt(lie.Qx), "optics_along"))
    push!(rows, ("Qy", _fmt(lie.Qy), "optics_along"))
    push!(rows, ("Qx_frac", _fmt(lie.Qx_frac), "periodic_optics mu/2pi"))
    push!(rows, ("Qy_frac", _fmt(lie.Qy_frac), "periodic_optics mu/2pi"))
    push!(rows, ("p_over_q_p", string(lie.rat.p), "best p/q of Qx_frac, q<=100"))
    push!(rows, ("p_over_q_q", string(lie.rat.q), "best p/q of Qx_frac, q<=100"))
    push!(rows, ("p_over_q_err", _fmt(lie.rat.err), "|q Qx_frac - p|"))
    push!(rows, ("Qy_same_q_err", _fmt(lie.qy_err), "|q Qy_frac - round| on the horizontal q"))
    push!(rows, ("max_abs_f3", _fmt(lie.max_f3), "max |coefficient| of f3"))
    push!(rows, ("max_abs_f4", _fmt(lie.max_f4), "max |coefficient| of f4"))
    push!(rows, ("linear_Tr", _fmt(lie.linear_tr), "2 cos(q mu_x)"))
    push!(rows, ("origin_Tr", _fmt(lie.origin_tr), "horizontal Tr of Rprod at x0 = 0, degree 3"))
    push!(rows, ("delta_m", _fmt(lie.delta_m), "regimes.linear"))
    push!(rows, ("aperture_m", _fmt(APERTURE), "survival cut at the ring start"))
    push!(rows, ("track_turns", string(nturns), "exact_track_turn"))
    push!(rows, ("track_ds_m", _fmt(TRACK_DS), "exact_track_turn step"))
    push!(rows, ("n_launches", string(length(tracked)), "survival grid"))
    n_ok = count(r -> r.survived, tracked)
    push!(rows, ("n_survived", string(n_ok), "exact_track_turn"))
    push!(rows, ("survival_fraction", _fmt(n_ok / length(tracked)), "n_survived / n_launches"))
    push!(rows, ("nf_magnetic_ok", lie.nf_magnetic_ok ? "1" : "0", "dragt_finn_nf6d of the magnetic factor"))
    push!(rows, ("nf_rf_ok", lie.nf_rf_ok ? "1" : "0", "dragt_finn_nf6d of one_turn_factor_rf"))
    if lie.nf_rf_ok
        nf = lie.nf
        for i in 1:3, j in 1:3
            push!(rows, ("c_$(i)$(j)", _fmt(nf.c[i, j]), "optional RF h4 = sum_{i<=j} c_ij J_i J_j"))
        end
        push!(rows, ("nf_f3_before", _fmt(nf.f3_before), "optional RF normal form"))
        push!(rows, ("nf_f3_after", _fmt(nf.f3_after), "optional RF cubic residual"))
        push!(rows, ("nf_f4_nonres_after", _fmt(nf.f4_nonres_after), "optional RF |f4 - h4|"))
        push!(rows, ("nf_Qs", _fmt(abs(nf.mu[3]) / (2π)), "optional RF |mu_s|/2pi"))
    end
    for b in bounds
        tag = "g$(b.degree)_d$(b.delta)_s$(b.side)"
        push!(rows, ("$(tag)_x_m", _fmt(b.side * b.result.x), "Guignard first |Tr|=2 boundary, stable side"))
        push!(rows, ("$(tag)_crossed", b.result.crossed ? "1" : "0", "1 if |Tr| leaves (-2, 2) inside the search"))
        push!(rows, ("$(tag)_reentrant", b.result.reentrant ? "1" : "0", "1 if a later grid sample is stable again"))
        push!(rows, ("$(tag)_x_largest_stable_m", _fmt(b.side * b.result.x_largest_stable),
                     "largest stable grid sample; equals the first boundary when the stable set is an interval"))
        push!(rows, ("$(tag)_Tr", _fmt(b.result.tr_stable.tr_x), "horizontal trace on the stable side of the first bracket"))
    end
    for (y, delta, side) in _edge_keys(tracked)
        e = _survival_edge(tracked, y, delta, side)
        tag = "surv_y$(y)_d$(delta)_s$(side)"
        push!(rows, ("$(tag)_last_x_m", _fmt(e.last_surviving), "outermost exact survivor connected to x = 0"))
        if e.first_lost !== nothing
            push!(rows, ("$(tag)_first_lost_x_m", _fmt(e.first_lost), "first exact loss on that ray"))
        end
    end
    open(path, "w") do io
        println(io, "quantity,value,source")
        for (q, v, src) in rows
            println(io, _csv_escape(q), ",", v, ",", _csv_escape(src))
        end
    end
    return path
end

function _edge_keys(tracked)
    keys = Tuple{Float64,Float64,Int}[]
    seen = Set{Tuple{Float64,Float64,Int}}()
    for r in tracked
        for side in (1, -1)
            key = (r.y, r.delta, side)
            key in seen && continue
            push!(seen, key)
            push!(keys, key)
        end
    end
    sort!(keys)
    return keys
end

function _write_guignard_csv(path, bounds)
    open(path, "w") do io
        println(io, "degree,delta,side,x0,tr_horizontal,tr_block,stable,kind")
        for b in bounds
            for s in b.result.samples
                x0 = b.side * s.x
                st = _stable_trace(s.t) ? 1 : 0
                @printf(io, "%d,%.8e,%d,%.8e,%.12e,%.12e,%d,sample\n",
                        b.degree, b.delta, b.side, x0, s.t.tr_x, s.t.tr_block, st)
            end
            @printf(io, "%d,%.8e,%d,%.8e,%.12e,%.12e,%d,boundary\n",
                    b.degree, b.delta, b.side, b.side * b.result.x,
                    b.result.tr_stable.tr_x, b.result.tr_stable.tr_block,
                    _stable_trace(b.result.tr_stable) ? 1 : 0)
        end
    end
    return path
end

function _write_survival_csv(path, tracked, nturns)
    open(path, "w") do io
        println(io, "x0,y0,delta,turns,survived,reason,aperture_m,nturns_requested,ds_m")
        for r in tracked
            @printf(io, "%.8e,%.8e,%.8e,%d,%d,%s,%.8e,%d,%.8e\n",
                    r.x, r.y, r.delta, r.turns, r.survived ? 1 : 0, r.reason,
                    APERTURE, nturns, TRACK_DS)
        end
    end
    return path
end

function _write_phase_csv(path, tracked)
    open(path, "w") do io
        println(io, "x0,turn,x,px")
        for r in tracked
            r.portrait || continue
            for (n, (x, px)) in enumerate(zip(r.phase_x, r.phase_px))
                @printf(io, "%.8e,%d,%.8e,%.8e\n", r.x, n - 1, x, px)
            end
        end
    end
    return path
end

function _sx(x, lo, hi, left, right)
    return left + (float(x) - lo) / (hi - lo) * (right - left)
end

function _sy(y, lo, hi, top, bottom)
    return bottom - (float(y) - lo) / (hi - lo) * (bottom - top)
end

function _poly_points(xs, ys, xlo, xhi, ylo, yhi, left, right, top, bottom)
    segs = Vector{String}[]
    cur = String[]
    for (x, y) in zip(xs, ys)
        if !(isfinite(x) && isfinite(y))
            isempty(cur) || push!(segs, cur)
            cur = String[]
            continue
        end
        push!(cur, @sprintf("%.2f,%.2f", _sx(x, xlo, xhi, left, right), _sy(y, ylo, yhi, top, bottom)))
    end
    isempty(cur) || push!(segs, cur)
    return segs
end

function _write_da_svg(path, bounds, tracked)
    width = 1100.0
    height = 540.0
    # Left: horizontal trace at δ = 0, |x| <= 25 mm.
    # Right: exact survival at y = 0.
    open(path, "w") do io
        @printf(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" viewBox=\"0 0 %.0f %.0f\">\n",
                width, height, width, height)
        println(io, "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>")
        println(io, "<text x=\"40\" y=\"28\" font-family=\"sans-serif\" font-size=\"16\">Guignard trace and exact survival, NoEdge, own sextupoles</text>")

        lL, lR, lT, lB = 70.0, 520.0, 56.0, 470.0
        xlo, xhi = -0.025, 0.025
        ys = Float64[]
        for b in bounds
            b.delta == 0.0 || continue
            for s in b.result.samples
                ax = b.side * s.x
                abs(ax) <= 0.025 || continue
                isfinite(s.t.tr_x) || continue
                push!(ys, s.t.tr_x)
            end
        end
        push!(ys, 2.0)
        ylo = minimum(ys) - 0.08 * (maximum(ys) - minimum(ys) + 1.0e-6)
        yhi = maximum(ys) + 0.08 * (maximum(ys) - minimum(ys) + 1.0e-6)
        if yhi - ylo < 1.0e-4
            ylo -= 5.0e-5
            yhi += 5.0e-5
        end
        println(io, "<rect x=\"$lL\" y=\"$lT\" width=\"$(lR - lL)\" height=\"$(lB - lT)\" fill=\"none\" stroke=\"#222\"/>")
        y2 = _sy(2.0, ylo, yhi, lT, lB)
        @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#666\" stroke-dasharray=\"4 3\"/>\n", lL, y2, lR, y2)
        @printf(io, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"12\" font-family=\"sans-serif\">|Tr|=2</text>\n", lR - 70, y2 - 4)
        colors = Dict(2 => "#1f4e79", 3 => "#c45911")
        for b in bounds
            b.delta == 0.0 || continue
            xs = Float64[]
            trs = Float64[]
            for s in sort(b.result.samples, by = s -> b.side * s.x)
                ax = b.side * s.x
                abs(ax) <= 0.025 || continue
                push!(xs, ax)
                push!(trs, s.t.tr_x)
            end
            for seg in _poly_points(xs, trs, xlo, xhi, ylo, yhi, lL, lR, lT, lB)
                length(seg) < 2 && continue
                @printf(io, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.6\" points=\"%s\"/>\n",
                        colors[b.degree], join(seg, " "))
            end
            if b.result.crossed && abs(b.result.x) <= 0.025
                xb = _sx(b.side * b.result.x, xlo, xhi, lL, lR)
                @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-dasharray=\"2 2\"/>\n",
                        xb, lT, xb, lB, colors[b.degree])
            end
        end
        println(io, "<text x=\"70\" y=\"500\" font-size=\"12\" font-family=\"sans-serif\" fill=\"#1f4e79\">degree 2</text>")
        println(io, "<text x=\"160\" y=\"500\" font-size=\"12\" font-family=\"sans-serif\" fill=\"#c45911\">degree 3</text>")
        println(io, "<text x=\"250\" y=\"500\" font-size=\"12\" font-family=\"sans-serif\">δ = 0, raw samples</text>")
        @printf(io, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"12\" font-family=\"sans-serif\" text-anchor=\"middle\">x0 [m]</text>\n",
                (lL + lR) / 2, lB + 36)
        @printf(io, "<text x=\"24\" y=\"%.2f\" font-size=\"12\" font-family=\"sans-serif\" transform=\"rotate(-90 24 %.2f)\">Tr of horizontal Rprod</text>\n",
                (lT + lB) / 2, (lT + lB) / 2)
        for x in (-0.02, -0.01, 0.0, 0.01, 0.02)
            px = _sx(x, xlo, xhi, lL, lR)
            @printf(io, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"11\" text-anchor=\"middle\" font-family=\"sans-serif\">%.0f mm</text>\n",
                    px, lB + 16, 1000 * x)
        end

        rL, rR, rT, rB = 620.0, 1060.0, 56.0, 470.0
        println(io, "<rect x=\"$rL\" y=\"$rT\" width=\"$(rR - rL)\" height=\"$(rB - rT)\" fill=\"none\" stroke=\"#222\"/>")
        deltas = unique(sort!(collect(r.delta for r in tracked if r.y == 0.0)))
        nrow = max(length(deltas), 1)
        for (i, d) in enumerate(deltas)
            cy = rT + (i - 0.5) * (rB - rT) / nrow
            @printf(io, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"11\" font-family=\"sans-serif\" text-anchor=\"end\">δ=%.1e</text>\n",
                    rL - 8, cy + 4, d)
            for r in tracked
                (r.y == 0.0 && r.delta == d) || continue
                px = _sx(r.x, -APERTURE, APERTURE, rL, rR)
                col = r.survived ? "#2e7d32" : "#c62828"
                @printf(io, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"3.2\" fill=\"%s\"/>\n", px, cy, col)
            end
        end
        for b in bounds
            b.delta == 0.0 || continue
            b.result.crossed || continue
            abs(b.result.x) > APERTURE && continue
            xb = _sx(b.side * b.result.x, -APERTURE, APERTURE, rL, rR)
            dash = b.degree == 2 ? "5 3" : "1 3"
            @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-dasharray=\"%s\"/>\n",
                    xb, rT, xb, rB, colors[b.degree], dash)
        end
        @printf(io, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"12\" font-family=\"sans-serif\" text-anchor=\"middle\">x0 [m], y = 0</text>\n",
                (rL + rR) / 2, rB + 36)
        println(io, "<text x=\"620\" y=\"500\" font-size=\"12\" font-family=\"sans-serif\" fill=\"#2e7d32\">survived</text>")
        println(io, "<text x=\"710\" y=\"500\" font-size=\"12\" font-family=\"sans-serif\" fill=\"#c62828\">lost</text>")
        println(io, "<text x=\"770\" y=\"500\" font-size=\"12\" font-family=\"sans-serif\">dashed lines: |Tr|=2 at δ=0</text>")
        for x in (-0.04, -0.02, 0.0, 0.02, 0.04)
            px = _sx(x, -APERTURE, APERTURE, rL, rR)
            @printf(io, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"11\" text-anchor=\"middle\" font-family=\"sans-serif\">%.0f</text>\n",
                    px, rB + 16, 1000 * x)
        end
        println(io, "</svg>")
    end
    return path
end

function _write_phase_svg(path, tracked)
    portraits = [r for r in tracked if r.portrait && !isempty(r.phase_x)]
    n = length(portraits)
    n == 0 && (write(path, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"200\" height=\"40\"><text y=\"20\">no portraits</text></svg>"); return path)
    width = 280.0 * n
    height = 320.0
    open(path, "w") do io
        @printf(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" viewBox=\"0 0 %.0f %.0f\">\n",
                width, height, width, height)
        println(io, "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>")
        println(io, "<text x=\"16\" y=\"24\" font-family=\"sans-serif\" font-size=\"16\">Phase portraits at the ring start, exact track, y = 0, δ = 0</text>")
        for (i, r) in enumerate(portraits)
            left = 40.0 + (i - 1) * 280.0
            right = left + 210.0
            top = 50.0
            bottom = 270.0
            xs = r.phase_x
            pxs = r.phase_px
            finite = [(x, p) for (x, p) in zip(xs, pxs) if isfinite(x) && isfinite(p)]
            println(io, "<rect x=\"$left\" y=\"$top\" width=\"$(right - left)\" height=\"$(bottom - top)\" fill=\"none\" stroke=\"#222\"/>")
            @printf(io, "<text x=\"%.2f\" y=\"42\" font-size=\"13\" font-family=\"sans-serif\">x0 = %.1f mm, %s</text>\n",
                    left, 1000 * r.x, r.reason)
            if length(finite) >= 2
                xmax = min(max(maximum(abs(p[1]) for p in finite), 1.0e-4), 5 * APERTURE) * 1.15
                pmax = min(max(maximum(abs(p[2]) for p in finite), 1.0e-6), 0.05) * 1.15
                segs = _poly_points([p[1] for p in finite], [p[2] for p in finite],
                                    -xmax, xmax, -pmax, pmax, left, right, top, bottom)
                for seg in segs
                    length(seg) < 2 && continue
                    @printf(io, "<polyline fill=\"none\" stroke=\"#1f4e79\" stroke-width=\"1.2\" points=\"%s\"/>\n",
                            join(seg, " "))
                end
            end
            @printf(io, "<text x=\"%.2f\" y=\"292\" font-size=\"11\" text-anchor=\"middle\" font-family=\"sans-serif\">x [m]</text>\n",
                    (left + right) / 2)
        end
        println(io, "</svg>")
    end
    return path
end

function _write_md(path, lie, bounds, tracked, nturns, seconds)
    lines = String[]
    push!(lines, "# Dynamic aperture")
    push!(lines, "")
    push!(lines, "NoEdge (`FringeModel()`) on the step-5 own sextupoles. The one-turn object is `ring_lie_factor`. Guignard uses `factor_to_map` truncated at degree 2 and degree 3, then `jacobian_of_map` on the linear orbit of `F.R`. Survival is `exact_track_turn` at the ring start. The two are overlaid and not averaged.")
    push!(lines, "")
    push!(lines, "Search window for the trace: 0 to $(SEARCH_M) m. Aperture cut: $(APERTURE) m. Turns: $nturns. Step: $(TRACK_DS) m. q ≤ $(Q_MAX). Runtime $(@sprintf("%.1f", seconds)) s.")
    push!(lines, "")
    if !isempty(STOPS)
        push!(lines, "## Stops")
        push!(lines, "")
        for s in STOPS
            push!(lines, "- $s")
        end
        push!(lines, "")
    end
    push!(lines, "## Lie readings")
    push!(lines, "")
    push!(lines, "Tunes are `optics_along` (integer part included) and μ/2π from `periodic_optics` of the factor. `{Q_x}` is approximated by $(_best_label(lie.rat)), the fraction with q ≤ $(Q_MAX) that minimises |q {Q_x} − p|.")
    push!(lines, "")
    push!(lines, _md_table(
        ["quantity", "value"],
        [
            ["Qx", _fmt6(lie.Qx)],
            ["Qy", _fmt6(lie.Qy)],
            ["Qx fractional", _fmt6(lie.Qx_frac)],
            ["Qy fractional", _fmt6(lie.Qy_frac)],
            ["p/q", _best_label(lie.rat)],
            ["|q Qx_frac − p|", _fmt6(lie.rat.err)],
            ["|q Qy_frac − nearest integer|", _fmt6(lie.qy_err)],
            ["max |f3|", _fmt6(lie.max_f3)],
            ["max |f4|", _fmt6(lie.max_f4)],
            ["2 cos(q μx)", _fmt6(lie.linear_tr)],
            ["Tr at x0 = 0, degree 3", _fmt6(lie.origin_tr)],
            ["αx", _fmt6(lie.alpha_x)],
            ["βx [m]", _fmt6(lie.beta_x)],
            ["coupling", _fmt6(lie.coupling)],
            ["δ_m", _fmt6(lie.delta_m)],
        ],
    ))
    push!(lines, "")
    push!(lines, "Distance of the fractional tune to 1/3 and to 1/4, and to the nearest multiple of those families:")
    push!(lines, "")
    push!(lines, _md_table(
        ["plane", "to 1/3", "to 1/4", "nearest third", "nearest quarter"],
        [
            ["x", _fmt6(_circle_dist(lie.Qx_frac, 1 / 3)), _fmt6(_circle_dist(lie.Qx_frac, 1 / 4)),
             _fmt6(_family_dist(lie.Qx_frac, 3)), _fmt6(_family_dist(lie.Qx_frac, 4))],
            ["y", _fmt6(_circle_dist(lie.Qy_frac, 1 / 3)), _fmt6(_circle_dist(lie.Qy_frac, 1 / 4)),
             _fmt6(_family_dist(lie.Qy_frac, 3)), _fmt6(_family_dist(lie.Qy_frac, 4))],
        ],
    ))
    push!(lines, "")
    push!(lines, "Distance of n_x Q_x + n_y Q_y to the nearest integer. The tunes are the full `optics_along` values.")
    push!(lines, "")
    body = Any[]
    for (nx, ny, name) in RESONANCES
        push!(body, [name, _fmt6(_dist_int(nx * lie.Qx + ny * lie.Qy))])
    end
    push!(lines, _md_table(["resonance", "|n·Q − integer|"], body))
    push!(lines, "")
    push!(lines, "Largest |coefficients| of f3:")
    push!(lines, "")
    mono = Any[]
    for (c, e) in lie.leading_f3
        push!(mono, [_mono_name(e), _fmt6(c)])
    end
    push!(lines, _md_table(["monomial", "coefficient"], mono))
    push!(lines, "")
    if lie.nf_magnetic_ok
        push!(lines, "`dragt_finn_nf6d` accepted the magnetic factor.")
    else
        push!(lines, "The magnetic factor was not accepted by `dragt_finn_nf6d`: $(_one_line(lie.nf_magnetic_error)). The longitudinal block is a shear, so this is the documented rejection (fewer than three elliptic pairs), not a second tune reading.")
    end
    push!(lines, "")
    if lie.nf_rf_ok
        nf = lie.nf
        push!(lines, "Optional normal form of `one_turn_factor_rf` at the linear-regime voltage $(lie.voltage_MV) MV, harmonic $(lie.harmonic), BELL phase π. This h4 is not an input of the Guignard scan. `Q_s = |μ_s| / 2π = $(_fmt6(abs(nf.mu[3]) / (2π)))`. max |f3| before the homological step is $(_fmt6(nf.f3_before)); the cubic residual is $(_fmt6(nf.f3_after)); |f4 − h4| is $(_fmt6(nf.f4_nonres_after)).")
        push!(lines, "")
        crow = Any[]
        for i in 1:3
            push!(crow, [_fmt6(nf.c[i, j]) for j in 1:3])
        end
        push!(lines, _md_table(["c_i1", "c_i2", "c_i3"], crow))
    else
        push!(lines, "Optional RF normal form did not return: $(_one_line(lie.nf_rf_error))")
    end
    push!(lines, "")
    push!(lines, "## Guignard boundary")
    push!(lines, "")
    push!(lines, "The linear orbit is R applied to (x0, 0, 0, 0, 0, δ). Rprod is the product of the horizontal 2×2 blocks of the one-turn Jacobian. The boundary below is the stable side of the first crossing of |Tr| = 2. A sample that is stable again after that crossing is listed as re-entrant and does not move the boundary.")
    push!(lines, "")
    grow = Any[]
    for b in bounds
        push!(grow, [
            string(b.degree),
            _fmt6(b.delta),
            b.side > 0 ? "+" : "−",
            _status_word(b.result),
            _fmt6(b.side * b.result.x),
            _fmt_tr(b.result.tr_stable.tr_x),
            b.result.reentrant ? _fmt6(b.side * b.result.x_largest_stable) : "same",
        ])
    end
    push!(lines, _md_table(
        ["degree", "δ", "sign", "status", "first x0 [m]", "Tr stable side", "largest stable sample [m]"],
        grow,
    ))
    push!(lines, "")
    push!(lines, "The reported boundary is that first crossing, the edge of the stable interval that contains x0 = 0. Where the largest stable grid sample differs, the trace has become unstable and then stable again. That later sample is not used as the boundary.")
    push!(lines, "")
    push!(lines, "## Tracking survival")
    push!(lines, "")
    n_ok = count(r -> r.survived, tracked)
    push!(lines, "Launches: $(length(tracked)). Survived all $nturns turns: $n_ok. Fraction: $(_fmt6(n_ok / length(tracked))). The cut is |x| < $(APERTURE) m and |y| < $(APERTURE) m at the ring start, and a finite state vector. The initial condition is the Guignard launch (x0, 0, y0, 0, 0, δ), not a matched Courant–Snyder momentum.")
    if n_ok == length(tracked)
        push!(lines, "")
        push!(lines, "Every launch survived. The outermost survivor is the edge of this grid. No dynamic-aperture loss was observed inside the aperture cut.")
    end
    push!(lines, "")
    srow = Any[]
    for (y, delta, side) in _edge_keys(tracked)
        e = _survival_edge(tracked, y, delta, side)
        lost = e.first_lost === nothing ? "none" : _fmt6(e.first_lost)
        re = isempty(e.reentered) ? "no" : join(_fmt6.(e.reentered), ", ")
        push!(srow, [_fmt6(y), _fmt6(delta), side > 0 ? "+" : "−", _fmt6(e.last_surviving), lost, re])
    end
    push!(lines, _md_table(
        ["y [m]", "δ", "sign", "last survivor connected to 0 [m]", "first loss [m]", "re-entrant survivors"],
        srow,
    ))
    push!(lines, "")
    off = [r for r in tracked if r.y != 0.0]
    if !isempty(off)
        n_off = count(r -> r.survived, off)
        push!(lines, "Off-axis launches (y ≠ 0, δ = 0): $(length(off)), survived $n_off.")
        push!(lines, "")
    end
    push!(lines, "Phase portraits are the raw turn samples of the exact track at x0 = 5, 10, 20 and 40 mm, y = 0, δ = 0.")
    push!(lines, "")
    push!(lines, "## Figures")
    push!(lines, "")
    push!(lines, "- `07_da.svg` — degree-2 and degree-3 traces at δ = 0, and the exact survival samples at y = 0. Vertical lines are the |Tr| = 2 boundaries that fall inside the aperture window. Nothing is interpolated between a trace and a survival point.")
    push!(lines, "- `07_phase.svg` — phase portraits.")
    push!(lines, "")
    push!(lines, "Tables: `07_numbers.csv`, `07_guignard.csv`, `07_survival.csv`, `07_phase.csv`.")
    push!(lines, "")
    write(path, join(lines, "\n"))
    return path
end

_best_label(rat) = "$(rat.p)/$(rat.q)"

_one_line(s) = replace(string(s), "\n" => " ")

# --- main --------------------------------------------------------------------

function _scan_all(maps, R, q, deltas)
    amplitudes = _amplitudes()
    amplitudes[1] == 0.0 || error("amplitude grid does not start at 0")
    abs(amplitudes[end] - SEARCH_M) < 1.0e-15 ||
        error("amplitude grid does not end at the search limit")
    bounds = NamedTuple[]
    for (degree, M) in maps
        for delta in deltas
            for side in (1, -1)
                println("Guignard degree $degree  δ $delta  side $side")
                flush(stdout)
                trace_at(a) = _guignard_trace(M, R, side * a, delta, q)
                result = _first_boundary(trace_at, amplitudes)
                for s in result.samples
                    if _trace_disagrees(s.t)
                        stop!(string(
                            "degree $degree, δ $delta, x0 $(side * s.x): ",
                            "horizontal-block product Tr $(s.t.tr_x) and the horizontal block of the 6×6 product Tr $(s.t.tr_block) disagree on |Tr| < 2",
                        ))
                        break
                    end
                end
                @printf("  %s  x = %.6e  Tr = %.8f\n",
                        _status_word(result), side * result.x, result.tr_stable.tr_x)
                flush(stdout)
                push!(bounds, (degree = degree, delta = delta, side = side, result = result))
            end
        end
    end
    return bounds
end

function _try_magnetic_nf(F)
    try
        nf = dragt_finn_nf6d(F)
        return (ok = true, nf = nf, error = "")
    catch err
        return (ok = false, nf = nothing, error = sprint(showerror, err))
    end
end

function _try_rf_nf(lat, regime, C)
    cav = RFCavity("RF", regime.voltage_MV, regime.harmonic, BELL_PHI_S)
    try
        Frf = one_turn_factor_rf(lat, cav, lat.ref, C)
        nf = dragt_finn_nf6d(Frf)
        return (ok = true, nf = nf, error = "")
    catch err
        return (ok = false, nf = nothing, error = sprint(showerror, err))
    end
end

function main()
    t_all = time()
    _self_check()
    long = _long_run()
    nturns = _turns()
    println(long ? "LONG_DA track, $nturns turns" : "Short track, $nturns turns")
    flush(stdout)

    println("Loading ", closed_lattice_path())
    flush(stdout)
    bare = load_closed_lattice()
    _require_bare(bare)
    own = _load_own(STRENGTH_JSON)
    regime = _load_linear(REGIME_JSON)
    regime.delta_m > 0 || error("delta_m is not positive")
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
    @printf("Own k2  SXD %.8e  SXB %.8e  SYB %.8e\n",
            own.k2["SXD"], own.k2["SXB"], own.k2["SYB"])
    flush(stdout)

    println("ring_lie_factor NoEdge")
    flush(stdout)
    F = ring_lie_factor(lat, FringeModel())
    ref = lat.ref
    C = path_length(lat)
    Rlin = ring_matrix(lat)
    gap = maximum(abs, F.R - Rlin)
    if gap > MATRIX_TOL
        stop!("NoEdge Lie matrix differs from ring_matrix by $gap")
    end
    po = periodic_optics(F.R, ref, C)
    po.x.stable && po.y.stable || error("periodic optics of the sextupole lattice is unstable")
    along = optics_along(lat, po; ds = DS)
    if abs(_frac_diff(along.Qx, po.x.mu / (2π))) > TUNE_TOL ||
       abs(_frac_diff(along.Qy, po.y.mu / (2π))) > TUNE_TOL
        stop!("optics_along and μ/2π differ by $(_frac_diff(along.Qx, po.x.mu / (2π))), $(_frac_diff(along.Qy, po.y.mu / (2π)))")
    end
    ch = chromaticity_orders(factor_to_map(F))
    if abs(_frac_diff(ch.x[1], po.x.mu / (2π))) > TUNE_TOL ||
       abs(_frac_diff(ch.y[1], po.y.mu / (2π))) > TUNE_TOL
        stop!("chromaticity_orders tune and μ/2π differ by $(_frac_diff(ch.x[1], po.x.mu / (2π))), $(_frac_diff(ch.y[1], po.y.mu / (2π)))")
    end
    if maximum(abs, F.R[6, 1:5]) > MATRIX_TOL || abs(F.R[6, 6] - 1) > MATRIX_TOL
        stop!("the linear map does not keep δ fixed: row 6 = $(F.R[6, :])")
    end

    println("factor_to_map")
    flush(stdout)
    M2 = factor_to_map(F; order = 2)
    M3 = factor_to_map(F; order = 3)
    cubic2 = _max_map_degree(M2, 3)
    cubic2 == 0.0 || stop!("degree-2 map still has a cubic coefficient $cubic2")
    if maximum(abs, map_linear_part(M3) - F.R) > MATRIX_TOL
        stop!("degree-3 map linear part differs from F.R")
    end
    J0 = jacobian_of_map(M3, zeros(6))
    if maximum(abs, J0 - F.R) > MATRIX_TOL
        stop!("jacobian_of_map at the origin differs from F.R by $(maximum(abs, J0 - F.R))")
    end

    rat = _best_pq(po.x.mu / (2π))
    q = rat.q
    linear_tr = 2 * cos(q * po.x.mu)
    origin = _guignard_trace(M3, F.R, 0.0, 0.0, q)
    if abs(origin.tr_x - linear_tr) > TRACE_TOL
        stop!("horizontal Tr at x0 = 0 is $(origin.tr_x); 2 cos(q μx) is $linear_tr")
    end
    qy_err = abs(q * (po.y.mu / (2π)) - round(q * (po.y.mu / (2π))))
    @printf("Qx %.10f  Qy %.10f  p/q = %d/%d  |q Qx - p| = %.3e  origin Tr = %.10f\n",
            along.Qx, along.Qy, rat.p, rat.q, rat.err, origin.tr_x)
    flush(stdout)

    println("normal form")
    flush(stdout)
    nf_m = _try_magnetic_nf(F)
    nf_r = _try_rf_nf(lat, regime, C)
    println("  magnetic nf ok = ", nf_m.ok, nf_m.ok ? "" : "  " * nf_m.error)
    println("  RF nf ok = ", nf_r.ok)
    flush(stdout)

    lie = (
        C = C,
        Qx = along.Qx,
        Qy = along.Qy,
        Qx_frac = po.x.mu / (2π),
        Qy_frac = po.y.mu / (2π),
        rat = rat,
        qy_err = qy_err,
        max_f3 = maxabs(F.f3),
        max_f4 = maxabs(F.f4),
        leading_f3 = _leading(F.f3, 3, 8),
        linear_tr = linear_tr,
        origin_tr = origin.tr_x,
        delta_m = regime.delta_m,
        alpha_x = po.x.alpha,
        beta_x = po.x.beta,
        coupling = po.coupling,
        voltage_MV = regime.voltage_MV,
        harmonic = regime.harmonic,
        nf_magnetic_ok = nf_m.ok,
        nf_magnetic_error = nf_m.error,
        nf_rf_ok = nf_r.ok,
        nf_rf_error = nf_r.error,
        nf = nf_r.nf,
    )

    bounds = _scan_all(((2, M2), (3, M3)), F.R, q, (0.0, regime.delta_m, -regime.delta_m))

    launches = _launches(regime.delta_m, long)
    println("Exact track: $(length(launches)) launches, $nturns turns, aperture $APERTURE m")
    flush(stdout)
    exact_track_turn([1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0], lat; ds = TRACK_DS)
    tracked = NamedTuple[]
    for (i, L) in enumerate(launches)
        @printf("track %d/%d  x % .6e  y % .6e  δ % .6e\n", i, length(launches), L.x, L.y, L.delta)
        flush(stdout)
        trk = _track_exact([L.x, 0.0, L.y, 0.0, 0.0, L.delta], lat, nturns, APERTURE;
                           record = L.portrait, progress = long)
        @printf("  %s  turns %d\n", trk.reason, trk.turns)
        flush(stdout)
        push!(tracked, (
            x = L.x, y = L.y, delta = L.delta, portrait = L.portrait,
            turns = trk.turns, survived = trk.survived, reason = trk.reason,
            phase_x = trk.x, phase_px = trk.px,
        ))
    end

    mkpath(OUT_DIR)
    _write_numbers(joinpath(OUT_DIR, "07_numbers.csv"), lie, bounds, tracked, nturns)
    _write_guignard_csv(joinpath(OUT_DIR, "07_guignard.csv"), bounds)
    _write_survival_csv(joinpath(OUT_DIR, "07_survival.csv"), tracked, nturns)
    _write_phase_csv(joinpath(OUT_DIR, "07_phase.csv"), tracked)
    _write_da_svg(joinpath(OUT_DIR, "07_da.svg"), bounds, tracked)
    _write_phase_svg(joinpath(OUT_DIR, "07_phase.svg"), tracked)
    seconds = time() - t_all
    _write_md(joinpath(OUT_DIR, "07_dynamic_aperture.md"), lie, bounds, tracked, nturns, seconds)
    n_ok = count(r -> r.survived, tracked)
    @printf("Done in %.1f s. Survived %d / %d. Stops: %d\n", seconds, n_ok, length(tracked), length(STOPS))
    flush(stdout)
    return nothing
end

main()
