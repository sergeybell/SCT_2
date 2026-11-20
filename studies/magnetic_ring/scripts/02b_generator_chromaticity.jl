# Step 2b. Chromaticity of the closed magnetic ring from the Lie generators.
#
#   julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/02b_generator_chromaticity.jl
#
# The one-turn factor is ring_lie_factor with no RF. The normalising matrix
# is built in src/generator_chromaticity.jl. symplectic_eigen_A is called
# once, only to record what it does with the magnetic shear. Its output is
# not a chromaticity.

import Pkg

const STUDY_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(STUDY_DIR)

using Bell
using Bell.NormalForm: symplectic_eigen_A
using MagneticRing
using LinearAlgebra
using Printf

include(joinpath(STUDY_DIR, "src", "generator_chromaticity.jl"))
using .GeneratorChromaticity

const OUT_DIR = joinpath(STUDY_DIR, "output")
const SERIES_ORDER = 4
const ORDER_STEPS = (1.0e-4, 1.0e-5, 1.0e-6)
const MATRIX_TOL = 1.0e-8
const NORM_TOL = 1.0e-8
# Settled finite difference versus chromaticity_series.
const AGREE_XI1 = 1.0e-4
const AGREE_XI2 = 1.0e-3
# Relative error at which a candidate factor is accepted for ξ1.
const FACTOR_REL = 1.0e-3
# Relative error at which an extra ξ2 factor (1, 2, 1/2, …) is called clear.
const XI2_REL = 1.0e-2

function main()
    _self_check()
    path = closed_lattice_path()
    println("Loading ", path)
    flush(stdout)
    lattice = load_closed_lattice()
    _require_bare_ring(lattice)

    println("Linear optics for the sextupole scan")
    flush(stdout)
    C = path_length(lattice)
    R0 = ring_matrix(lattice)
    po0 = periodic_optics(R0, lattice.ref, C)
    opt = optics_along(lattice, po0; ds = 0.05)
    rows = sextupole_rows(lattice, opt)
    isempty(rows) && error("the closed lattice has no sextupole slot")
    slot = rows[argmax([abs(r.beta_x_D) for r in rows])]
    k2 = _select_k2(slot)
    println("Sextupole slot ", slot.index, " ", slot.name,
            "  s0 = ", slot.s0, " m  βx D = ", slot.beta_x_D, "  k2 = ", k2)
    flush(stdout)

    println("Case k2 = 0")
    flush(stdout)
    bare = _read_case(lattice, "k2=0")
    println("Case one sextupole")
    flush(stdout)
    excited = _read_case(with_one_sextupole(lattice, slot.index, k2), "one_sextupole")

    factor = _match_factor(bare.gen.x1.c, bare.gen.y1.c, _xi1(bare, :x), _xi1(bare, :y))
    stops = _stops(bare, excited, factor)
    xi2 = _xi2_verdict(bare, excited, factor, stops)

    mkpath(OUT_DIR)
    paths = (
        report = joinpath(OUT_DIR, "02b_chromaticity.md"),
        table = joinpath(OUT_DIR, "02b_chromaticity.csv"),
        setup = joinpath(OUT_DIR, "02b_normalisation.csv"),
        draft = joinpath(OUT_DIR, "bell_issue_chromaticity_from_generators.md"),
    )
    _write_table(paths.table, bare, excited, factor, slot, k2)
    _write_setup(paths.setup, bare, excited, factor, slot, k2, xi2, stops)
    _write_report(paths.report, path, lattice, rows, slot, k2, bare, excited, factor, xi2, stops)
    _write_draft(paths.draft, bare, factor, xi2, stops)

    _print_summary(bare, excited, factor, slot, k2, xi2, stops, paths)
    return nothing
end

function _self_check()
    A = dispersion_shear([1.1, -0.3, 0.2, 0.05])
    symplectic_error(A) < 1.0e-12 || error("dispersion shear is not symplectic")
    B = twiss_normalisation(3.0, 0.4, 2.0, -0.2)
    symplectic_error(B) < 1.0e-12 || error("Twiss matrix is not symplectic")
    symplectic_error(A * B) < 1.0e-12 || error("full normalisation is not symplectic")
    R = transverse_rotation(1.2, 2.3)
    Jx = (Z[1]^2 + Z[2]^2) / 2
    maxabs(hom(Jx * Z[6] - pull_generator(Jx * Z[6], R), 3)) < 1.0e-12 ||
        error("J δ is not invariant under the transverse rotation")
    println("Normalisation self-check passed")
    flush(stdout)
    return nothing
end

function _require_bare_ring(lattice)
    ref = lattice.ref
    ref.particle.name == "deuteron" || error("reference particle is $(ref.particle.name)")
    ref.Ekin_MeV == 270.0 || error("kinetic energy is $(ref.Ekin_MeV) MeV")
    abs(total_bend_angle(lattice) - 2π) < 1.0e-9 ||
        error("total bend is not 2π: $(total_bend_angle(lattice))")
    any(e -> e isa RFCavity, lattice.elements) && error("the lattice contains an RF cavity")
    for e in lattice.elements
        e isa Sextupole && e.k2 != 0 &&
            error("sextupole $(e.name) has k2 = $(e.k2); the bare case needs k2 = 0")
    end
    n_mag = length(magnetic_lattice(lattice).elements)
    n_mag == length(lattice.elements) ||
        error("magnetic_lattice dropped $(length(lattice.elements) - n_mag) elements")
    return nothing
end

function _select_k2(slot)
    denom = slot.L * slot.beta_x_D
    abs(denom) < 1.0e-4 &&
        error("sextupole $(slot.name) at index $(slot.index) has |βx D| = $(abs(slot.beta_x_D)); no feed-down lever")
    for m in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0)
        if abs(m * denom) / (4π) >= 0.25
            return m
        end
    end
    return 50.0
end

function _read_case(lattice, label)
    mag = magnetic_lattice(lattice)
    length(mag.elements) == length(lattice.elements) ||
        error("$label: magnetic_lattice dropped elements")
    println("  ring_matrix")
    flush(stdout)
    R = ring_matrix(lattice)
    C = path_length(lattice)
    println("  ring_lie_factor")
    flush(stdout)
    F = ring_lie_factor(mag)
    println("  factor_to_map, dispersion_orders, chromaticity_orders")
    flush(stdout)
    M = factor_to_map(F)
    orders = Dict{Float64,Any}(dd => chromaticity_orders(M; dd = dd) for dd in ORDER_STEPS)
    dord = dispersion_orders(M)
    po = periodic_optics(F.R, lattice.ref, C)
    po_R = periodic_optics(R, lattice.ref, C)
    println("  generators")
    flush(stdout)
    gen = generator_chromaticity(F, po)
    println("  chromaticity_series order ", SERIES_ORDER)
    flush(stdout)
    series = chromaticity_series(mag, SERIES_ORDER)
    eigen_text = label == "k2=0" ? _probe_eigen(F.R) : ""
    println("  done ", label)
    flush(stdout)
    return (
        label = label,
        R = R,
        F = F,
        orders = orders,
        series = series,
        dord = dord,
        po = po,
        po_R = po_R,
        gen = gen,
        eigen_text = eigen_text,
        matrix_gap = maximum(abs, F.R .- R),
        dispersion_gap = maximum(abs, Float64.(po.D4) .- Float64.(dord.D1)),
        optics_gap = _optics_gap(po, po_R),
    )
end

function _optics_gap(a, b)
    return max(
        abs(a.x.beta - b.x.beta),
        abs(a.x.alpha - b.x.alpha),
        abs(a.x.mu - b.x.mu),
        abs(a.y.beta - b.y.beta),
        abs(a.y.alpha - b.y.alpha),
        abs(a.y.mu - b.y.mu),
        maximum(abs, Float64.(a.D4) .- Float64.(b.D4)),
    )
end

function _probe_eigen(R)
    ev = eigen(R)
    n_ell = count(λ -> imag(λ) > 1.0e-8, ev.values)
    values = join([@sprintf("%+.6e%+.6eim", real(λ), imag(λ)) for λ in ev.values], ", ")
    thrown = ""
    returned = ""
    try
        A, mu = symplectic_eigen_A(R)
        returned = @sprintf("returned without throwing, mu = (%.6e, %.6e, %.6e), |A|max = %.6e",
                            mu[1], mu[2], mu[3], maximum(abs, A))
    catch err
        thrown = sprint(showerror, err)
    end
    return """
    Elliptic pairs with Im λ > 1e-8: $n_ell.
    Eigenvalues: $values.
    symplectic_eigen_A: $(isempty(thrown) ? returned : "threw: " * thrown)
    """
end

_xi1(case, plane) = plane == :x ? case.orders[1.0e-6].x[2] : case.orders[1.0e-6].y[2]
_xi2(case, plane) = plane == :x ? case.orders[1.0e-6].x[3] : case.orders[1.0e-6].y[3]

function _orders_xi(case, plane, which, dd)
    ch = case.orders[dd]
    block = plane == :x ? ch.x : ch.y
    return which == :xi1 ? block[2] : block[3]
end

function _series_xi(case, plane, which)
    v = plane == :x ? case.series.xi_x : case.series.xi_y
    return Float64(which == :xi1 ? v[1] : v[2])
end

function _series_Q(case, plane)
    v = plane == :x ? case.series.x : case.series.y
    return Float64(v[1])
end

# Candidates are the clear multiples of 1 and 1/π. The run picks one.
# It does not interpolate.
const FACTOR_CANDIDATES = (
    ("-1/(2π)", -1 / (2π)),
    ("+1/(2π)", +1 / (2π)),
    ("-1/π", -1 / π),
    ("+1/π", +1 / π),
    ("-1/(4π)", -1 / (4π)),
    ("+1/(4π)", +1 / (4π)),
    ("-2/π", -2 / π),
    ("+2/π", +2 / π),
    ("-1", -1.0),
    ("+1", +1.0),
    ("-2", -2.0),
    ("+2", +2.0),
    ("-1/2", -0.5),
    ("+1/2", +0.5),
)

function _rel(pred, obs)
    return abs(pred - obs) / max(abs(obs), 1.0e-12)
end

function _match_factor(cx, cy, xi_x, xi_y)
    best = nothing
    for (name, fac) in FACTOR_CANDIDATES
        rx = _rel(fac * cx, xi_x)
        ry = _rel(fac * cy, xi_y)
        score = max(rx, ry)
        if best === nothing || score < best.score
            best = (name = name, factor = fac, score = score, rel_x = rx, rel_y = ry)
        end
    end
    accepted = best.score <= FACTOR_REL
    ratio_x = xi_x / cx
    ratio_y = xi_y / cy
    return merge(best, (accepted = accepted, ratio_x = ratio_x, ratio_y = ratio_y, cx = cx, cy = cy))
end

function _predict(factor, c1, c2)
    # μ_extra = c / factor_to_mu, with ξ1 = factor * c1 = (1/2π) dμ/dδ.
    # A term c2 J δ² then contributes factor * c2 * δ² to Q, and the
    # second derivative is twice that coefficient.
    return (xi1 = factor * c1, xi2 = 2 * factor * c2)
end

const XI2_CLEAR = (
    ("1", 1.0),
    ("-1", -1.0),
    ("2", 2.0),
    ("-2", -2.0),
    ("1/2", 0.5),
    ("-1/2", -0.5),
    ("4", 4.0),
    ("-4", -4.0),
    ("1/4", 0.25),
    ("-1/4", -0.25),
)

function _scaled_gap(pred, obs)
    return abs(pred - obs) / max(abs(obs), abs(pred), 1.0e-4)
end

function _clear_ratio(pred_x, pred_y, obs_x, obs_y)
    abs(pred_x) < 1.0e-8 && abs(pred_y) < 1.0e-8 && abs(obs_x) < 1.0e-6 && abs(obs_y) < 1.0e-6 &&
        return (
            matched = true, name = "1", factor = 1.0, score = 0.0,
            reason = "both sides are at round-off",
        )
    best = nothing
    for (name, fac) in XI2_CLEAR
        sx = _scaled_gap(fac * pred_x, obs_x)
        sy = _scaled_gap(fac * pred_y, obs_y)
        score = max(sx, sy)
        if best === nothing || score < best.score
            best = (name = name, factor = fac, score = score)
        end
    end
    matched = best.score <= XI2_REL
    return merge(best, (matched = matched, reason = matched ? "" : "no clear factor"))
end

function _xi2_verdict(bare, excited, factor, stops)
    if !factor.accepted
        return (
            worked = false,
            extra = "not applied",
            detail = "ξ1 factor was not accepted, so ξ2 was not converted",
        )
    end
    bare_x = _predict(factor.factor, bare.gen.x1.c, bare.gen.x2.c)
    bare_y = _predict(factor.factor, bare.gen.y1.c, bare.gen.y2.c)
    versus_orders = _clear_ratio(bare_x.xi2, bare_y.xi2, _xi2(bare, :x), _xi2(bare, :y))
    versus_series = _clear_ratio(bare_x.xi2, bare_y.xi2, _series_xi(bare, :x, :xi2), _series_xi(bare, :y, :xi2))
    ex_x = _predict(factor.factor, excited.gen.x1.c, excited.gen.x2.c)
    ex_y = _predict(factor.factor, excited.gen.y1.c, excited.gen.y2.c)
    ex_orders = _clear_ratio(ex_x.xi2, ex_y.xi2, _xi2(excited, :x), _xi2(excited, :y))
    bell_split = stops.bare_xi2 || stops.excited_xi2
    same = versus_orders.matched && versus_orders.name == "1" &&
           ex_orders.matched && ex_orders.name == "1"
    worked = same && !bell_split
    detail = if bell_split
        "chromaticity_orders and chromaticity_series disagree on ξ2; no column was selected. The generator column is f4 − ½[f3, F3]"
    elseif same
        "f4 − ½[f3, F3] reproduces d²Q/dδ² with the ξ1 factor, on the bare ring and with one sextupole. f4 + ½[f3, F3] is stored as xi2_plus and is not the reported value"
    elseif versus_orders.matched
        "bare-ring J δ² matches d²Q/dδ² after an extra factor $(versus_orders.name); the sextupole case extra factor is $(ex_orders.matched ? ex_orders.name : "none")"
    else
        "f4 − ½[f3, F3] does not match chromaticity_orders within a clear factor"
    end
    return (
        worked = worked,
        extra = versus_orders.matched ? versus_orders.name : "none",
        detail = detail,
        versus_orders = versus_orders,
        versus_series = versus_series,
        excited_orders = ex_orders,
        pred_bare = (x = bare_x, y = bare_y),
        pred_excited = (x = ex_x, y = ex_y),
    )
end

function _reading_gap(case, which, atol)
    gaps = Float64[]
    flags = String[]
    for plane in (:x, :y)
        settled = which == :xi1 ? _xi1(case, plane) : _xi2(case, plane)
        series = _series_xi(case, plane, which)
        step_spread = maximum(abs, [(_orders_xi(case, plane, which, dd) - settled) for dd in ORDER_STEPS])
        gap = abs(settled - series)
        push!(gaps, gap)
        if gap > atol && gap > 5 * max(step_spread, 1.0e-12)
            push!(flags, "$plane $which |orders − series| = $gap, step spread = $step_spread")
        end
    end
    return (gap = maximum(gaps), disagree = !isempty(flags), flags = flags)
end

function _stops(bare, excited, factor)
    flags = String[]
    bare.matrix_gap <= MATRIX_TOL || push!(flags, "k2=0 |Lie R − ring_matrix| = $(bare.matrix_gap)")
    excited.matrix_gap <= MATRIX_TOL || push!(flags, "sextupole |Lie R − ring_matrix| = $(excited.matrix_gap)")
    bare.dispersion_gap <= MATRIX_TOL ||
        push!(flags, "k2=0 |periodic D − dispersion_orders| = $(bare.dispersion_gap)")
    excited.dispersion_gap <= MATRIX_TOL ||
        push!(flags, "sextupole |periodic D − dispersion_orders| = $(excited.dispersion_gap)")
    bare.optics_gap <= MATRIX_TOL ||
        push!(flags, "k2=0 periodic optics of Lie R and ring_matrix differ by $(bare.optics_gap)")
    dR = maximum(abs, bare.R .- excited.R)
    dR <= 1.0e-10 || push!(flags, "turning on one sextupole changed the one-turn matrix by $dR")
    for (case, tag) in ((bare, "k2=0"), (excited, "sextupole"))
        for plane in (:x, :y)
            q_orders = plane == :x ? case.orders[1.0e-4].x[1] : case.orders[1.0e-4].y[1]
            q_mu = plane == :x ? case.po.x.mu / 2π : case.po.y.mu / 2π
            q_series = _series_Q(case, plane)
            if abs(q_orders - q_mu) > MATRIX_TOL
                push!(flags, "$tag fractional tune, orders versus μ/2π, $plane, |diff| = $(abs(q_orders - q_mu))")
            end
            if abs(q_series - q_mu) > MATRIX_TOL
                push!(flags, "$tag fractional tune, series versus μ/2π, $plane, |diff| = $(abs(q_series - q_mu))")
            end
        end
    end
    b1 = _reading_gap(bare, :xi1, AGREE_XI1)
    b2 = _reading_gap(bare, :xi2, AGREE_XI2)
    e1 = _reading_gap(excited, :xi1, AGREE_XI1)
    e2 = _reading_gap(excited, :xi2, AGREE_XI2)
    for (tag, g) in (("k2=0", b1), ("k2=0", b2), ("sextupole", e1), ("sextupole", e2))
        g.disagree || continue
        for line in g.flags
            push!(flags, "$tag $line")
        end
    end
    norm_ok = _normalisation_ok(bare.gen) && _normalisation_ok(excited.gen)
    return (
        flags = flags,
        stopped = !isempty(flags),
        bare_xi1 = b1.disagree,
        bare_xi2 = b2.disagree,
        excited_xi1 = e1.disagree,
        excited_xi2 = e2.disagree,
        norm_ok = norm_ok,
        factor_accepted = factor.accepted,
    )
end

function _normalisation_ok(g)
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

function _fmt(x)
    return isfinite(x) ? @sprintf("%.8e", float(x)) : "nan"
end

function _fmt_rel(x)
    return isfinite(x) ? @sprintf("%.3e", float(x)) : "nan"
end

function _csv_escape(s)
    t = string(s)
    if occursin(',', t) || occursin('"', t)
        return "\"" * replace(t, '"' => "\"\"") * "\""
    end
    return t
end

function _gen_xi1(gen, plane, factor)
    c = plane == :x ? gen.x1.c : gen.y1.c
    return factor.accepted ? factor.factor * c : NaN
end

function _write_table(path, bare, excited, factor, slot, k2)
    open(path, "w") do io
        println(io, "case,plane,quantity,generator,orders_dd_1e-4,orders_dd_1e-5,orders_dd_1e-6,series,formula_centre,formula_integral")
        for (case, case_name) in ((bare, "k2=0"), (excited, "one_sextupole"))
            for (plane, sym) in ((:x, "x"), (:y, "y"))
                for which in (:xi1, :xi2)
                    # ξ2 is the direct Lie prediction of d²Q/dδ², before any
                    # extra factor. The report states whether one was required.
                    g = if which == :xi1
                        _gen_xi1(case.gen, plane, factor)
                    else
                        c = plane == :x ? case.gen.x2.c : case.gen.y2.c
                        factor.accepted ? 2 * factor.factor * c : NaN
                    end
                    @printf(io, "%s,%s,%s,%s,%s,%s,%s,%s,,\n",
                            case_name, sym, string(which), _fmt(g),
                            _fmt(_orders_xi(case, plane, which, 1.0e-4)),
                            _fmt(_orders_xi(case, plane, which, 1.0e-5)),
                            _fmt(_orders_xi(case, plane, which, 1.0e-6)),
                            _fmt(_series_xi(case, plane, which)))
                    which == :xi2 || continue
                    cplus = plane == :x ? case.gen.x2_plus.c : case.gen.y2_plus.c
                    gplus = factor.accepted ? 2 * factor.factor * cplus : NaN
                    @printf(io, "%s,%s,xi2_plus,%s,%s,%s,%s,%s,,\n",
                            case_name, sym, _fmt(gplus),
                            _fmt(_orders_xi(case, plane, :xi2, 1.0e-4)),
                            _fmt(_orders_xi(case, plane, :xi2, 1.0e-5)),
                            _fmt(_orders_xi(case, plane, :xi2, 1.0e-6)),
                            _fmt(_series_xi(case, plane, :xi2)))
                end
            end
        end
        fx = feeddown_xi(k2, slot.L, slot.beta_x, slot.D)
        fy = -feeddown_xi(k2, slot.L, slot.beta_y, slot.D)
        ix = feeddown_xi_integral(k2, slot.integral_beta_x_D)
        iy = -feeddown_xi_integral(k2, slot.integral_beta_y_D)
        for (sym, plane, formula, integ) in (("x", :x, fx, ix), ("y", :y, fy, iy))
            dg = _gen_xi1(excited.gen, plane, factor) - _gen_xi1(bare.gen, plane, factor)
            @printf(io, "one_sextupole_minus_k2=0,%s,delta_xi1,%s,%s,%s,%s,%s,%s,%s\n",
                    sym, _fmt(dg),
                    _fmt(_orders_xi(excited, plane, :xi1, 1.0e-4) - _orders_xi(bare, plane, :xi1, 1.0e-4)),
                    _fmt(_orders_xi(excited, plane, :xi1, 1.0e-5) - _orders_xi(bare, plane, :xi1, 1.0e-5)),
                    _fmt(_orders_xi(excited, plane, :xi1, 1.0e-6) - _orders_xi(bare, plane, :xi1, 1.0e-6)),
                    _fmt(_series_xi(excited, plane, :xi1) - _series_xi(bare, plane, :xi1)),
                    _fmt(formula), _fmt(integ))
        end
    end
    return path
end

function _write_setup(path, bare, excited, factor, slot, k2, xi2, stops)
    rows = Tuple{String,String}[]
    push!(rows, ("series_order", string(SERIES_ORDER)))
    push!(rows, ("factor_name", factor.name))
    push!(rows, ("factor_value", _fmt(factor.factor)))
    push!(rows, ("factor_accepted", string(factor.accepted)))
    push!(rows, ("factor_max_rel", _fmt_rel(factor.score)))
    push!(rows, ("ratio_xi1_over_c_x", _fmt(factor.ratio_x)))
    push!(rows, ("ratio_xi1_over_c_y", _fmt(factor.ratio_y)))
    push!(rows, ("xi2_extra_factor", xi2.extra))
    push!(rows, ("xi2_worked", string(xi2.worked)))
    push!(rows, ("stopped", string(stops.stopped)))
    push!(rows, ("sextupole_index", string(slot.index)))
    push!(rows, ("sextupole_name", slot.name))
    push!(rows, ("sextupole_k2", _fmt(k2)))
    push!(rows, ("sextupole_L_m", _fmt(slot.L)))
    push!(rows, ("sextupole_s0_m", _fmt(slot.s0)))
    push!(rows, ("sextupole_beta_x_m", _fmt(slot.beta_x)))
    push!(rows, ("sextupole_beta_y_m", _fmt(slot.beta_y)))
    push!(rows, ("sextupole_Dx_m", _fmt(slot.D)))
    for (tag, case) in (("k2=0", bare), ("sextupole", excited))
        g = case.gen
        push!(rows, ("$tag symplectic_error_A", _fmt(g.symplectic_error_A)))
        push!(rows, ("$tag symplectic_error_shear", _fmt(g.symplectic_error_shear)))
        push!(rows, ("$tag dispersion_column", _fmt(g.dispersion_column)))
        push!(rows, ("$tag shear_dispersion_column", _fmt(g.shear_dispersion_column)))
        push!(rows, ("$tag longitudinal_coupling", _fmt(g.longitudinal_coupling)))
        push!(rows, ("$tag rotation_residual", _fmt(g.rotation_residual)))
        push!(rows, ("$tag transverse_coupling", _fmt(g.transverse_coupling)))
        push!(rows, ("$tag delta_row", _fmt(g.delta_row)))
        push!(rows, ("$tag l_on_delta", _fmt(g.l_delta)))
        push!(rows, ("$tag max_l_f3", _fmt(g.max_l_f3)))
        push!(rows, ("$tag max_l_f4", _fmt(g.max_l_f4)))
        push!(rows, ("$tag c_Jx_delta", _fmt(g.x1.c)))
        push!(rows, ("$tag c_Jy_delta", _fmt(g.y1.c)))
        push!(rows, ("$tag q2_minus_p2_x", _fmt(g.x1.q2 - g.x1.p2)))
        push!(rows, ("$tag q2_minus_p2_y", _fmt(g.y1.q2 - g.y1.p2)))
        push!(rows, ("$tag cross_x", _fmt(g.x1.cross)))
        push!(rows, ("$tag cross_y", _fmt(g.y1.cross)))
        push!(rows, ("$tag c_Jx_delta2", _fmt(g.x2.c)))
        push!(rows, ("$tag c_Jy_delta2", _fmt(g.y2.c)))
        push!(rows, ("$tag c_Jx_delta2_plus", _fmt(g.x2_plus.c)))
        push!(rows, ("$tag c_Jy_delta2_plus", _fmt(g.y2_plus.c)))
        push!(rows, ("$tag c_Jx_delta2_bare_f4", _fmt(g.x2_bare.c)))
        push!(rows, ("$tag c_Jy_delta2_bare_f4", _fmt(g.y2_bare.c)))
        push!(rows, ("$tag homological_residual", _fmt(g.homological_residual)))
        push!(rows, ("$tag F3_max", _fmt(g.F3_max)))
        push!(rows, ("$tag kernel_dimension", string(g.kernel_dimension)))
        push!(rows, ("$tag smallest_kept_singular", _fmt(g.smallest_kept)))
        push!(rows, ("$tag matrix_gap", _fmt(case.matrix_gap)))
        push!(rows, ("$tag dispersion_gap", _fmt(case.dispersion_gap)))
    end
    open(path, "w") do io
        println(io, "quantity,value")
        for (q, v) in rows
            println(io, _csv_escape(q), ",", _csv_escape(v))
        end
    end
    return path
end

function _md_table(headers, rows)
    lines = String["| " * join(headers, " | ") * " |",
                   "| " * join(fill("---", length(headers)), " | ") * " |"]
    for r in rows
        push!(lines, "| " * join(r, " | ") * " |")
    end
    return join(lines, "\n")
end

function _plane_rows(case, factor)
    out = Vector{Vector{String}}()
    for (plane, sym) in ((:x, "x"), (:y, "y"))
        c1 = plane == :x ? case.gen.x1.c : case.gen.y1.c
        c2 = plane == :x ? case.gen.x2.c : case.gen.y2.c
        g1 = factor.accepted ? factor.factor * c1 : NaN
        g2 = factor.accepted ? 2 * factor.factor * c2 : NaN
        push!(out, [sym, "ξ1", _fmt(g1),
                    _fmt(_orders_xi(case, plane, :xi1, 1.0e-4)),
                    _fmt(_orders_xi(case, plane, :xi1, 1.0e-6)),
                    _fmt(_series_xi(case, plane, :xi1))])
        push!(out, [sym, "ξ2", _fmt(g2),
                    _fmt(_orders_xi(case, plane, :xi2, 1.0e-4)),
                    _fmt(_orders_xi(case, plane, :xi2, 1.0e-6)),
                    _fmt(_series_xi(case, plane, :xi2))])
    end
    return out
end

function _write_report(path, latpath, lattice, rows, slot, k2, bare, excited, factor, xi2, stops)
    g = bare.gen
    lines = String[]
    push!(lines, "# Chromaticity from the magnetic Lie generators")
    push!(lines, "")
    push!(lines, "Lattice file: `$latpath`. Reference: $(lattice.ref.particle.name), kinetic energy $(lattice.ref.Ekin_MeV) MeV, harmonic $(lattice.harmonic). No RF. The bare case has every sextupole at `k2 = 0`.")
    push!(lines, "")
    push!(lines, "The one-turn object is `ring_lie_factor` of `magnetic_lattice`. `chromaticity_orders` reads `factor_to_map` of that factor. `chromaticity_series` builds `one_turn_map_da` at order $SERIES_ORDER. The normalising matrix is not `symplectic_eigen_A`.")
    push!(lines, "")
    push!(lines, "## Normalisation")
    push!(lines, "")
    push!(lines, "`A` maps normal coordinates to canonical coordinates: dispersion shear, then the Courant–Snyder matrix of each plane at the ring start. The shear is `x → x + D δ`, `p_x → p_x + D' δ`, and `l → l − D' x + D p_x`, with the same pattern in `y`. Actions are `J = (q² + p²) / 2`.")
    push!(lines, "")
    push!(lines, _md_table(["quantity", "k2 = 0", "one sextupole"], [
        ["symplectic error of A", _fmt(g.symplectic_error_A), _fmt(excited.gen.symplectic_error_A)],
        ["symplectic error of the shear", _fmt(g.symplectic_error_shear), _fmt(excited.gen.symplectic_error_shear)],
        ["max |R_n[1:4, δ]|", _fmt(g.dispersion_column), _fmt(excited.gen.dispersion_column)],
        ["max |R_shear[1:4, δ]|", _fmt(g.shear_dispersion_column), _fmt(excited.gen.shear_dispersion_column)],
        ["max |R_n[l, 1:4]|", _fmt(g.longitudinal_coupling), _fmt(excited.gen.longitudinal_coupling)],
        ["rotation residual", _fmt(g.rotation_residual), _fmt(excited.gen.rotation_residual)],
        ["transverse coupling of R_n", _fmt(g.transverse_coupling), _fmt(excited.gen.transverse_coupling)],
        ["max |R_n[δ, 1:5]|", _fmt(g.delta_row), _fmt(excited.gen.delta_row)],
        ["(l|δ) in normal coordinates", _fmt(g.l_delta), _fmt(excited.gen.l_delta)],
        ["max |coef| of f3 monomials with l", _fmt(g.max_l_f3), _fmt(excited.gen.max_l_f3)],
        ["max |coef| of f4 monomials with l", _fmt(g.max_l_f4), _fmt(excited.gen.max_l_f4)],
        ["|q² − p²| coefficient, x, Jδ", _fmt(abs(g.x1.q2 - g.x1.p2)), _fmt(abs(excited.gen.x1.q2 - excited.gen.x1.p2))],
        ["|q p| coefficient, x, Jδ", _fmt(abs(g.x1.cross)), _fmt(abs(excited.gen.x1.cross))],
        ["|q² − p²| coefficient, y, Jδ", _fmt(abs(g.y1.q2 - g.y1.p2)), _fmt(abs(excited.gen.y1.q2 - excited.gen.y1.p2))],
        ["homological residual", _fmt(g.homological_residual), _fmt(excited.gen.homological_residual)],
        ["max |coef| of F3", _fmt(g.F3_max), _fmt(excited.gen.F3_max)],
        ["kernel dimension", string(g.kernel_dimension), string(excited.gen.kernel_dimension)],
        ["smallest kept singular value", _fmt(g.smallest_kept), _fmt(excited.gen.smallest_kept)],
        ["|Lie R − ring_matrix|", _fmt(bare.matrix_gap), _fmt(excited.matrix_gap)],
        ["|periodic D − dispersion_orders|", _fmt(bare.dispersion_gap), _fmt(excited.dispersion_gap)],
    ]))
    push!(lines, "")
    push!(lines, "Normalisation checks are within $NORM_TOL: $(stops.norm_ok).")
    push!(lines, "")
    push!(lines, "## First chromaticity")
    push!(lines, "")
    push!(lines, "The cubic generator in normal coordinates contains `c_x J_x δ + c_y J_y δ`. The map is `exp(:f:)`, and the angle that matches `rotation` by `+μ` is the generator `−μ J`. The overall factor was required to be one entry of a fixed list, the same number on both planes, and was compared with `chromaticity_orders` at step `1e-6` on the `k2 = 0` ring.")
    push!(lines, "")
    push!(lines, _md_table(["quantity", "value"], [
        ["c_x of J_x δ", _fmt(g.x1.c)],
        ["c_y of J_y δ", _fmt(g.y1.c)],
        ["ξ1_x / c_x", _fmt(factor.ratio_x)],
        ["ξ1_y / c_y", _fmt(factor.ratio_y)],
        ["candidate", factor.name],
        ["factor", _fmt(factor.factor)],
        ["max relative residual", _fmt_rel(factor.score)],
        ["accepted", factor.accepted ? "yes" : "no"],
    ]))
    push!(lines, "")
    if factor.accepted
        push!(lines, "Accepted factor: ξ1 = $(factor.name) × (coefficient of J δ). On this factor, `d²Q/dδ²` predicted from a term `c₂ J δ²` is `2 × $(factor.name) × c₂`.")
    else
        push!(lines, "No candidate matched both planes within relative error $FACTOR_REL. The generator column below is left empty. The ratios above are the raw readings.")
    end
    push!(lines, "")
    push!(lines, "`chromaticity_orders` is a central difference. Steps `1e-4` (the default), `1e-5`, and `1e-6` are all stored. The comparison with `chromaticity_series` uses the step `1e-6` once the three steps have settled. ξ1 is `dQ/dδ`. ξ2 is `d²Q/dδ²`, the third component of `chromaticity_orders` and the second derivative from `chromaticity_series`.")
    push!(lines, "")
    push!(lines, "## Bare ring, k2 = 0")
    push!(lines, "")
    push!(lines, _md_table(["plane", "quantity", "generator", "orders 1e-4", "orders 1e-6", "series"],
                           _plane_rows(bare, factor)))
    push!(lines, "")
    push!(lines, "Generator ξ2 in that table uses `c₂` from `f4 − ½[f3, F3]` and the direct prediction `2 × factor × c₂`. It is not rescaled. `F3` is the solution of `F3 − F3 ∘ r = −f3` on the non-resonant cubic. The row `ξ2 plus` below is the same prediction with the bracket added instead of subtracted.")
    push!(lines, "")
    push!(lines, _md_table(["plane", "ξ2 from f4 − ½[f3,F3]", "ξ2 from f4 + ½[f3,F3]", "orders 1e-6", "series"], [
        ["x", _fmt(factor.accepted ? 2 * factor.factor * g.x2.c : NaN),
         _fmt(factor.accepted ? 2 * factor.factor * g.x2_plus.c : NaN),
         _fmt(_xi2(bare, :x)), _fmt(_series_xi(bare, :x, :xi2))],
        ["y", _fmt(factor.accepted ? 2 * factor.factor * g.y2.c : NaN),
         _fmt(factor.accepted ? 2 * factor.factor * g.y2_plus.c : NaN),
         _fmt(_xi2(bare, :y)), _fmt(_series_xi(bare, :y, :xi2))],
    ]))
    push!(lines, "")
    push!(lines, _md_table(["plane", "c₂ of f4 − ½[f3,F3]", "c₂ of f4 + ½[f3,F3]", "c₂ in f4 alone", "|q² − p²| of the minus quartic"], [
        ["x", _fmt(g.x2.c), _fmt(g.x2_plus.c), _fmt(g.x2_bare.c), _fmt(abs(g.x2.q2 - g.x2.p2))],
        ["y", _fmt(g.y2.c), _fmt(g.y2_plus.c), _fmt(g.y2_bare.c), _fmt(abs(g.y2.q2 - g.y2.p2))],
    ]))
    push!(lines, "")
    push!(lines, "Second chromaticity: $(xi2.detail). Extra factor on top of the ξ1 conversion: $(xi2.extra).")
    push!(lines, "")
    push!(lines, "## One sextupole")
    push!(lines, "")
    push!(lines, "The slot is the existing sextupole with the largest `|β_x D_x|` at its centre. Centre values are the drift transport of the entrance Twiss across half the element. The horizontal formula is `+k2 L β_x D_x / (4π)`. The vertical formula is `−k2 L β_y D_x / (4π)`. `D_x` is the horizontal dispersion in both planes. The integral column uses `k2 ∫ β D_x ds / (4π)` with the same signs.")
    push!(lines, "")
    push!(lines, _md_table(["quantity", "value"], [
        ["element index", string(slot.index)],
        ["name", slot.name],
        ["s at entrance [m]", _fmt(slot.s0)],
        ["L [m]", _fmt(slot.L)],
        ["k2 [1/m³]", _fmt(k2)],
        ["βx at centre [m]", _fmt(slot.beta_x)],
        ["βy at centre [m]", _fmt(slot.beta_y)],
        ["Dx at centre [m]", _fmt(slot.D)],
        ["formula Δξx centre", _fmt(feeddown_xi(k2, slot.L, slot.beta_x, slot.D))],
        ["formula Δξy centre", _fmt(-feeddown_xi(k2, slot.L, slot.beta_y, slot.D))],
        ["formula Δξx integral", _fmt(feeddown_xi_integral(k2, slot.integral_beta_x_D))],
        ["formula Δξy integral", _fmt(-feeddown_xi_integral(k2, slot.integral_beta_y_D))],
    ]))
    push!(lines, "")
    push!(lines, "Other sextupoles stay at `k2 = 0`. Linear matrices of the bare and excited rings differ by $( _fmt(maximum(abs, bare.R .- excited.R)) ).")
    push!(lines, "")
    push!(lines, _md_table(["plane", "quantity", "generator", "orders 1e-4", "orders 1e-6", "series"],
                           _plane_rows(excited, factor)))
    push!(lines, "")
    eg = excited.gen
    push!(lines, "The same bracket test on the sextupole. The plus combination changes sign relative to `chromaticity_series`. The minus combination stays with it.")
    push!(lines, "")
    push!(lines, _md_table(["plane", "ξ2 from f4 − ½[f3,F3]", "ξ2 from f4 + ½[f3,F3]", "orders 1e-6", "series"], [
        ["x", _fmt(factor.accepted ? 2 * factor.factor * eg.x2.c : NaN),
         _fmt(factor.accepted ? 2 * factor.factor * eg.x2_plus.c : NaN),
         _fmt(_xi2(excited, :x)), _fmt(_series_xi(excited, :x, :xi2))],
        ["y", _fmt(factor.accepted ? 2 * factor.factor * eg.y2.c : NaN),
         _fmt(factor.accepted ? 2 * factor.factor * eg.y2_plus.c : NaN),
         _fmt(_xi2(excited, :y)), _fmt(_series_xi(excited, :y, :xi2))],
    ]))
    push!(lines, "")
    delta_rows = Vector{Vector{String}}()
    for (sym, plane, beta, integ_bD, sign) in (
            ("x", :x, slot.beta_x, slot.integral_beta_x_D, 1.0),
            ("y", :y, slot.beta_y, slot.integral_beta_y_D, -1.0),
        )
        formula = sign * feeddown_xi(k2, slot.L, beta, slot.D)
        integ = sign * feeddown_xi_integral(k2, integ_bD)
        dg = _gen_xi1(excited.gen, plane, factor) - _gen_xi1(bare.gen, plane, factor)
        dord = _orders_xi(excited, plane, :xi1, 1.0e-6) - _orders_xi(bare, plane, :xi1, 1.0e-6)
        dser = _series_xi(excited, plane, :xi1) - _series_xi(bare, plane, :xi1)
        push!(delta_rows, [sym, _fmt(dg), _fmt(dord), _fmt(dser), _fmt(formula), _fmt(integ)])
    end
    push!(lines, _md_table(["plane", "Δξ1 generator", "Δξ1 orders 1e-6", "Δξ1 series", "formula centre", "formula integral"],
                           delta_rows))
    push!(lines, "")
    push!(lines, "Sextupole survey, centre values, bare strengths:")
    push!(lines, "")
    survey = Vector{Vector{String}}()
    for r in rows
        push!(survey, [string(r.index), r.name, _fmt(r.s0), _fmt(r.L), _fmt(r.beta_x), _fmt(r.beta_y), _fmt(r.D),
                       r.index == slot.index ? "chosen" : ""])
    end
    push!(lines, _md_table(["index", "name", "s0 [m]", "L [m]", "βx [m]", "βy [m]", "Dx [m]", "role"], survey))
    push!(lines, "")
    push!(lines, "## Agreement")
    push!(lines, "")
    if stops.stopped
        push!(lines, "STOP. Two BELL readings of the same quantity disagree. No column was selected as the value.")
        push!(lines, "")
        for line in stops.flags
            push!(lines, "- $line")
        end
    else
        push!(lines, "On each case, `chromaticity_orders` at step `1e-6` and `chromaticity_series` agree within $AGREE_XI1 on ξ1 and within $AGREE_XI2 on ξ2. `ring_matrix` and the Lie factor `R` agree within $MATRIX_TOL, and so do `periodic_optics` dispersion and `dispersion_orders`. The finite-difference steps are printed so the truncation of the default step `1e-4` can be seen.")
    end
    push!(lines, "")
    push!(lines, "A generator value is this study's reading of the Lie factor. Where it disagrees with a BELL column, both numbers stay in the table.")
    push!(lines, "")
    push!(lines, "## What symplectic_eigen_A does on this factor")
    push!(lines, "")
    push!(lines, "Called on `F.R` of the bare ring and not used afterwards.")
    push!(lines, "")
    push!(lines, "```")
    push!(lines, strip(bare.eigen_text))
    push!(lines, "```")
    push!(lines, "")
    write(path, join(lines, "\n"))
    return path
end

function _write_draft(path, bare, factor, xi2, stops)
    g = bare.gen
    lines = String[]
    push!(lines, "# Magnetic chromaticity is an invariant of the RF-off Lie factor")
    push!(lines, "")
    push!(lines, "## Summary")
    push!(lines, "")
    push!(lines, "The magnetic one-turn factor has two elliptic planes and a longitudinal shear. `symplectic_eigen_A` asks for three elliptic pairs. `resonant_cubic` treats a cubic generator as having no pure invariant; that statement fits three oscillators. On the magnetic factor the cubic invariants include `J_x δ`, `J_y δ`, and `δ³`. Those `J δ` coefficients are the first chromaticity.")
    push!(lines, "")
    push!(lines, "## What the magnetic ring shows")
    push!(lines, "")
    push!(lines, "Lattice: `magnetic_2` with every sector bend named BM set to π/8, sextupoles at `k2 = 0`, no cavity. The closed file used here is `studies/magnetic_ring/lattices/magnetic_2_closed.bell`. The factor is `ring_lie_factor(magnetic_lattice(lat))`.")
    push!(lines, "")
    push!(lines, strip(bare.eigen_text))
    push!(lines, "")
    push!(lines, "A normalising matrix built as a symplectic dispersion shear (with the `l` row `l → l − D' x + D p_x`) times the Courant–Snyder matrix at the ring start has symplectic error $( _fmt(g.symplectic_error_A) ). After the pull, `max |R_n[1:4, δ]|` is $( _fmt(g.dispersion_column) ) and the transverse blocks match the rotations to $( _fmt(g.rotation_residual) ).")
    push!(lines, "")
    push!(lines, "Monomials that contain `l` have max |coefficient| $( _fmt(g.max_l_f3) ) in `f3` and $( _fmt(g.max_l_f4) ) in `f4`.")
    push!(lines, "")
    push!(lines, "Coefficients of `J δ`: horizontal $( _fmt(g.x1.c) ), vertical $( _fmt(g.y1.c) ). Against `chromaticity_orders` on this ring the conversion is ξ1 = $(factor.name) × c, relative residual $( _fmt_rel(factor.score) ), accepted: $(factor.accepted).")
    push!(lines, "")
    push!(lines, "The transverse homological operator, with `δ` held fixed, has kernel dimension $(g.kernel_dimension) at the cutoff used (analytic count is 8: `δ³`, `J δ`, and the same families times powers of `l`). `F3` solves `F3 − F3 ∘ r = −f3` on the non-resonant cubic, residual $( _fmt(g.homological_residual) ).")
    push!(lines, "")
    push!(lines, "Coefficients of `J δ²` in `f4 − ½[f3, F3]`: horizontal $( _fmt(g.x2.c) ), vertical $( _fmt(g.y2.c) ). In `f4 + ½[f3, F3]`: horizontal $( _fmt(g.x2_plus.c) ), vertical $( _fmt(g.y2_plus.c) ). With the same factor −1/(2π), the minus combination reproduces `d²Q/dδ²`. The plus combination does not: once a sextupole makes the bracket the leading term, the plus combination flips the sign of ξ2.")
    push!(lines, "")
    push!(lines, "Second chromaticity: $(xi2.detail).")
    push!(lines, "")
    if stops.stopped
        push!(lines, "On this run, two BELL readings disagreed. The study did not choose between them:")
        push!(lines, "")
        for line in stops.flags
            push!(lines, "- $line")
        end
        push!(lines, "")
    end
    push!(lines, "## Recommendation")
    push!(lines, "")
    push!(lines, "Do not pass the RF-off factor to `symplectic_eigen_A`. Document that the magnetic `(l, δ)` block is a shear. For that factor, read the first chromaticity from the invariant `J δ` in the cubic generator after a dispersion shear and a transverse Twiss normalisation. Read the second from `J δ²` in `f4 − ½[f3, F3]`, where `F3` solves `F3 − F3 ∘ r = −f3` on the transverse rotation and `δ` is a parameter. The plus sign in front of that bracket, as written for the three-oscillator normal form, flips the sextupole chromaticity on this ring. State that the kernel of the transverse homological operator includes `δ³` and `J δ`. `resonant_cubic`, as written, is the three-oscillator statement and drops this invariant.")
    push!(lines, "")
    push!(lines, "## Reproducer")
    push!(lines, "")
    push!(lines, "```julia")
    push!(lines, "using Bell")
    push!(lines, "using Bell.NormalForm: symplectic_eigen_A")
    push!(lines, "lat = load_structure(raw\"C:/GIT_REPOS/SCT_2/studies/magnetic_ring/lattices/magnetic_2_closed.bell\")")
    push!(lines, "F = ring_lie_factor(magnetic_lattice(lat))")
    push!(lines, "symplectic_eigen_A(F.R)   # three elliptic pairs are required; this R has two")
    push!(lines, "```")
    push!(lines, "")
    push!(lines, "The same lattice is `examples/lattices/magnetic_2.bell` with the BM angle set to `π/8` and with no cavity. SCT_2 study step `studies/magnetic_ring/scripts/02b_generator_chromaticity.jl` records the generator coefficients. This note recommends an API reading. It does not change BELL.")
    push!(lines, "")
    write(path, join(lines, "\n"))
    return path
end

function _print_summary(bare, excited, factor, slot, k2, xi2, stops, paths)
    g = bare.gen
    println()
    println("symplectic_error(A) = ", g.symplectic_error_A)
    println("dispersion column = ", g.dispersion_column)
    println("max |l| in f3, f4 = ", g.max_l_f3, "  ", g.max_l_f4)
    @printf("c(Jx δ) = %.8e    c(Jy δ) = %.8e\n", g.x1.c, g.y1.c)
    println("factor ", factor.name, " = ", factor.factor, "  accepted = ", factor.accepted,
            "  max rel = ", factor.score)
    @printf("ξ1 generator x, y = %.8e  %.8e\n",
            factor.accepted ? factor.factor * g.x1.c : NaN,
            factor.accepted ? factor.factor * g.y1.c : NaN)
    @printf("ξ1 orders   x, y = %.8e  %.8e\n", _xi1(bare, :x), _xi1(bare, :y))
    @printf("ξ1 series   x, y = %.8e  %.8e\n", _series_xi(bare, :x, :xi1), _series_xi(bare, :y, :xi1))
    @printf("ξ2 generator x, y = %.8e  %.8e\n",
            factor.accepted ? 2 * factor.factor * g.x2.c : NaN,
            factor.accepted ? 2 * factor.factor * g.y2.c : NaN)
    @printf("ξ2 orders   x, y = %.8e  %.8e\n", _xi2(bare, :x), _xi2(bare, :y))
    @printf("ξ2 series   x, y = %.8e  %.8e\n", _series_xi(bare, :x, :xi2), _series_xi(bare, :y, :xi2))
    println("ξ2: ", xi2.detail)
    fx = feeddown_xi(k2, slot.L, slot.beta_x, slot.D)
    fy = -feeddown_xi(k2, slot.L, slot.beta_y, slot.D)
    @printf("Δξ1 orders x, y = %.8e  %.8e    formula = %.8e  %.8e\n",
            _xi1(excited, :x) - _xi1(bare, :x), _xi1(excited, :y) - _xi1(bare, :y), fx, fy)
    println("sextupole ", slot.name, " index ", slot.index, " k2 = ", k2)
    if stops.stopped
        println("STOP")
        for line in stops.flags
            println("  ", line)
        end
    else
        println("BELL readings agree")
    end
    println("normalisation ok = ", stops.norm_ok)
    println("outputs:")
    for p in paths
        println("  ", p)
    end
    flush(stdout)
    return nothing
end

main()
