# Step 3. Second-order slip and chromaticity of the magnetic Lie factor,
# then the 6D normal form of the same ring with one thin cavity.
#
#   julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/03_lie_invariants.jl
#
# The default track covers two synchrotron periods and stops at 2000 turns.
# A longer track is the explicit flag --turns N.

import Pkg

const STUDY_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(STUDY_DIR)

using Bell
using Bell.NormalForm: dragt_finn_nf6d, rf_kick_taylor, symplectic_eigen_A
using MagneticRing
using Printf

include(joinpath(STUDY_DIR, "src", "generator_chromaticity.jl"))
include(joinpath(STUDY_DIR, "src", "readings.jl"))
using .GeneratorChromaticity
using .LieReadings

const OUT_DIR = joinpath(STUDY_DIR, "output")
const REGIME_JSON = joinpath(OUT_DIR, "02_rf_regimes.json")
const TRACK_DS = 0.1
const TRACK_TURNS_CAP = 2000
const TRACK_PERIODS = 2.0
# Below transition the notes use a stationary phase of 0. BELL's smooth
# tune is stable when eta0 * cos(phi_s) > 0, so the cavity phase is π.
const BELL_PHI_S = Float64(π)
# Step 2b: the generator of a rotation by +μ is −μ J, and ξ1 = −c / (2π).
const XI_FACTOR = -1 / (2π)
const SLIP_TOL = 1.0e-8
const XI1_TOL = 1.0e-4
const XI2_TOL = 1.0e-3

function main()
    _self_check()
    path = closed_lattice_path()
    println("Loading ", path)
    flush(stdout)
    lattice = load_closed_lattice()
    _require_bare_ring(lattice)
    ref = lattice.ref
    C = path_length(lattice)

    println("Magnetic one-turn factor")
    flush(stdout)
    mag = _magnetic_readings(lattice, ref, C)

    println("Linear regime from ", REGIME_JSON)
    flush(stdout)
    regime = read_linear_regime(REGIME_JSON)
    rf = _rf_setup(lattice, ref, C, mag, regime)

    println("One-turn factor of the ring with the thin cavity")
    flush(stdout)
    F = ring_lie_factor(rf.lattice)
    kicks = _kick_comparison(mag.F, F, rf.cav, ref, C)

    println("Normal form")
    flush(stdout)
    nf_result = _normal_form(F)
    if !nf_result.ok
        println("Normal form failed. Magnetic tables will still be written.")
        println(nf_result.error)
        flush(stdout)
    end

    mkpath(OUT_DIR)
    paths = (
        report = joinpath(OUT_DIR, "03_lie_invariants.md"),
        numbers = joinpath(OUT_DIR, "03_numbers.csv"),
        actions = joinpath(OUT_DIR, "03_actions.csv"),
        figure = joinpath(OUT_DIR, "03_actions.svg"),
        draft = joinpath(OUT_DIR, "03_bell_issue_one_turn_factor_rf.md"),
    )
    # Magnetic readings and the normal form are on disk before the long track.
    _write_numbers(paths.numbers, mag, rf, kicks, nf_result, nothing)
    _write_report(paths.report, path, mag, rf, kicks, nf_result, nothing, "")
    _write_draft(paths.draft, kicks, nf_result, rf)

    track = nothing
    track_error = ""
    if nf_result.ok
        nturns = _turn_count(nf_result.Qs)
        ax, ay = _amplitudes(mag.po, regime.delta_m)
        z0 = _launch(mag.po, regime.delta_m, ax, ay)
        println("Tracking ", nturns, " turns, ds = ", TRACK_DS, " m")
        flush(stdout)
        try
            track = _track(rf.lattice, ref, C, z0, nf_result.nf, mag.po, nturns, ax, ay)
        catch err
            track_error = sprint(showerror, err)
            println("Tracking failed. Tables already written will be updated with this error.")
            println(track_error)
            flush(stdout)
        end
    end
    _write_numbers(paths.numbers, mag, rf, kicks, nf_result, track)
    if track !== nothing
        _write_actions(paths.actions, track)
        _write_svg(paths.figure, track)
    end
    _write_report(paths.report, path, mag, rf, kicks, nf_result, track, track_error)
    _write_draft(paths.draft, kicks, nf_result, rf)

    _print_summary(mag, rf, kicks, nf_result, track, paths, track_error)
    return nothing
end

function _self_check()
    g = 1.5
    g2 = g^2
    expected = 0.3 - 0.2 / g2 + (3 * g2 - 1) / (2 * g2^2)
    abs(eta1_formula(0.2, 0.3, g) - expected) < 1e-15 ||
        error("eta1 formula does not match its expression")
    beta, alpha, a = 4.0, 0.5, 1.0e-3
    J = a^2 / (2 * beta)
    px = -alpha / beta * a
    abs(courant_snyder_action(a, px, beta, alpha) - J) < 1e-18 ||
        error("Courant-Snyder action at the crest is not a^2 / (2 beta)")
    println("Reading self-check passed")
    flush(stdout)
    return nothing
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

function _magnetic_readings(lattice, ref, C)
    maglat = magnetic_lattice(lattice)
    length(maglat.elements) == length(lattice.elements) ||
        error("magnetic_lattice dropped elements from the bare ring")
    println("  ring_lie_factor")
    flush(stdout)
    F = ring_lie_factor(maglat)
    println("  factor_to_map, path_expansion, chromaticity_orders")
    flush(stdout)
    M = factor_to_map(F)
    pe = path_expansion(M, ref, C)
    orders = (
        fine = chromaticity_orders(M; dd = 1.0e-6),
        coarse = chromaticity_orders(M; dd = 1.0e-4),
    )
    R = ring_matrix(lattice)
    po = periodic_optics(F.R, ref, C)
    po.x.stable && po.y.stable || error("magnetic periodic optics is unstable")
    println("  generator_chromaticity")
    flush(stdout)
    gen = generator_chromaticity(F, po)
    note = eta1_formula(pe.alpha0, pe.alpha1, ref.gamma)
    c1, c2, _ = velocity_expansion(ref)
    bell_piece = -c2 + c1^2
    note_piece = (3 * ref.gamma^2 - 1) / (2 * ref.gamma^4)
    return (
        F = F,
        pe = pe,
        orders = orders,
        po = po,
        gen = gen,
        matrix_gap = maximum(abs, F.R .- R),
        eta1_formula = note,
        eta1_minus_formula = pe.eta1 - note,
        bell_eta1_formula = pe.eta1_formula,
        eta1_minus_bell = pe.eta1 - pe.eta1_formula,
        formula_minus_bell = note - pe.eta1_formula,
        velocity_piece_minus_note = bell_piece - note_piece,
        C = C,
        ref = ref,
    )
end

function _rf_setup(lattice, ref, C, mag, regime)
    mag.pe.eta0 < 0 ||
        error("eta0 = $(mag.pe.eta0) is not below transition; phi_s = π would be unstable in BELL")
    mag.pe.eta0 * cos(BELL_PHI_S) > 0 ||
        error("BELL stationary phase π is not stable for this eta0")
    lat, cav = append_thin_cavity(lattice, regime.voltage_MV, regime.harmonic, BELL_PHI_S)
    abs(path_length(lat) - C) <= 1e-12 || error("the thin cavity changed the path length")
    count(e -> e isa RFCavity, lat.elements) == 1 || error("expected one cavity")
    Qs_note = synchrotron_tune_note(ref, mag.pe.eta0, regime.voltage_MV, regime.harmonic)
    Qs_smooth = synchrotron_tune(ref, mag.pe.eta0, regime.voltage_MV, regime.harmonic, BELL_PHI_S)
    isfinite(Qs_note) && isfinite(Qs_smooth) || error("synchrotron tune is not finite")
    return (
        lattice = lat,
        cav = cav,
        regime = regime,
        Qs_note = Qs_note,
        Qs_smooth = Qs_smooth,
        Qs_json = regime.Qs,
    )
end

function _kick_comparison(Fmag, F, cav, ref, C)
    # ring_lie_factor(lat_with_cavity) already contains the cavity.
    # one_turn_factor_rf would call ring_lie_factor on that same lattice and
    # then concatenate rf_kick_taylor. That second concatenate is measured
    # here. one_turn_factor_rf itself is not called.
    println("  cavity factor from rf_kick_taylor, for the double-count measurement")
    flush(stdout)
    Frf = dragt_finn_factor(rf_kick_taylor(cav, ref, C); name = cav.name)
    one = concatenate(Fmag, Frf; name = "magnetic+one_kick")
    two = concatenate(F, Frf; name = "cavity_lattice+second_kick")
    mu_one = _safe_mu(F.R)
    mu_two = _safe_mu(two.R)
    qs_one = mu_one.ok ? abs(mu_one.mu[3]) / (2π) : NaN
    qs_two = mu_two.ok ? abs(mu_two.mu[3]) / (2π) : NaN
    return (
        R_gap = maximum(abs, F.R .- one.R),
        f3_gap = maxabs(F.f3 - one.f3),
        f4_gap = maxabs(F.f4 - one.f4),
        mu_one = mu_one,
        mu_two = mu_two,
        Qs_one = qs_one,
        Qs_two = qs_two,
        Qs_ratio = (isfinite(qs_one) && qs_one > 0) ? qs_two / qs_one : NaN,
    )
end

function _safe_mu(R)
    try
        _, mu = symplectic_eigen_A(R)
        return (ok = true, mu = mu, error = "")
    catch err
        return (ok = false, mu = fill(NaN, 3), error = sprint(showerror, err))
    end
end

function _normal_form(F)
    local nf
    try
        nf = dragt_finn_nf6d(F)
    catch err
        return (
            ok = false,
            nf = nothing,
            error = sprint(showerror, err),
            Qs = NaN,
        )
    end
    return (
        ok = true,
        nf = nf,
        error = "",
        Qs = abs(nf.mu[3]) / (2π),
    )
end

function _turn_count(Qs)
    requested = _requested_turns()
    if requested !== nothing
        return requested
    end
    Qs > 0 || return TRACK_TURNS_CAP
    return min(TRACK_TURNS_CAP, max(ceil(Int, TRACK_PERIODS / Qs), 50))
end

function _requested_turns()
    i = findfirst(==("--turns"), ARGS)
    i === nothing && return nothing
    i < length(ARGS) || error("--turns needs an integer")
    n = parse(Int, ARGS[i + 1])
    n >= 1 || error("--turns must be at least 1, got $n")
    return n
end

function _amplitudes(po, delta_m)
    disp = abs(po.D4[1]) * delta_m
    ax = disp > 0 ? 2 * disp : 1.0e-5
    return ax, ax
end

function _launch(po, delta, ax, ay)
    Dx, Dpx, Dy, Dpy = po.D4
    x = ax + Dx * delta
    px = -po.x.alpha / po.x.beta * ax + Dpx * delta
    y = ay + Dy * delta
    py = -po.y.alpha / po.y.beta * ay + Dpy * delta
    return [x, px, y, py, 0.0, delta]
end

function _record_point!(store, i, turn, z, nf, po)
    J = normal_form_actions(nf, z)
    xb, pxb = dispersion_free_transverse(z[1], z[2], po.D4[1], po.D4[2], z[6])
    yb, pyb = dispersion_free_transverse(z[3], z[4], po.D4[3], po.D4[4], z[6])
    store.turn[i] = turn
    store.z[:, i] = z
    store.Jx_nf[i] = J[1]
    store.Jy_nf[i] = J[2]
    store.Js_nf[i] = J[3]
    store.Jx_cs[i] = courant_snyder_action(z[1], z[2], po.x.beta, po.x.alpha)
    store.Jy_cs[i] = courant_snyder_action(z[3], z[4], po.y.beta, po.y.alpha)
    store.Jx_beta[i] = courant_snyder_action(xb, pxb, po.x.beta, po.x.alpha)
    store.Jy_beta[i] = courant_snyder_action(yb, pyb, po.y.beta, po.y.alpha)
    return nothing
end

function _track(lattice, ref, C, z0, nf, po, nturns, ax, ay)
    n = nturns + 1
    store = (
        turn = Vector{Int}(undef, n),
        z = Matrix{Float64}(undef, 6, n),
        Jx_nf = Vector{Float64}(undef, n),
        Jy_nf = Vector{Float64}(undef, n),
        Js_nf = Vector{Float64}(undef, n),
        Jx_cs = Vector{Float64}(undef, n),
        Jy_cs = Vector{Float64}(undef, n),
        Jx_beta = Vector{Float64}(undef, n),
        Jy_beta = Vector{Float64}(undef, n),
    )
    z = collect(Float64, z0)
    t0 = time()
    for i in 1:n
        if i > 1
            z = track_turn_exact(z, lattice, ref, C; ds = TRACK_DS)
        end
        _record_point!(store, i, i - 1, z, nf, po)
        if i == 2 || i == n || (i - 1) % 100 == 0
            @printf("  turn %5d / %d    elapsed %7.1f s    Jx_nf %.6e    Jx_cs %.6e    delta %.6e\n",
                    i - 1, nturns, time() - t0, store.Jx_nf[i], store.Jx_cs[i], z[6])
            flush(stdout)
        end
    end
    stats = (
        Jx_nf = _stat(store.Jx_nf),
        Jy_nf = _stat(store.Jy_nf),
        Js_nf = _stat(store.Js_nf),
        Jx_cs = _stat(store.Jx_cs),
        Jy_cs = _stat(store.Jy_cs),
        Jx_beta = _stat(store.Jx_beta),
        Jy_beta = _stat(store.Jy_beta),
        delta = _stat(store.z[6, :]),
    )
    return merge(store, (
        nturns = nturns,
        ds = TRACK_DS,
        ax = ax,
        ay = ay,
        elapsed = time() - t0,
        stats = stats,
        z0 = z0,
    ))
end

function _stat(v)
    n = length(v)
    m = sum(v) / n
    lo, hi = extrema(v)
    rms = sqrt(sum(x -> (x - m)^2, v) / n)
    span = hi - lo
    rel = abs(m) > 0 ? span / abs(m) : NaN
    return (mean = m, minimum = lo, maximum = hi, span = span, rms = rms, span_over_abs_mean = rel)
end

function _xi(gen, plane, which)
    if which == :xi1
        c = plane == :x ? gen.x1.c : gen.y1.c
        return XI_FACTOR * c
    end
    c = plane == :x ? gen.x2.c : gen.y2.c
    return 2 * XI_FACTOR * c
end

function _orders_xi(orders, plane, which)
    block = plane == :x ? orders.x : orders.y
    return which == :xi1 ? block[2] : block[3]
end

function _fmt(x)
    return isfinite(x) ? @sprintf("%.8e", float(x)) : "nan"
end

function _fmt6(x)
    return isfinite(x) ? @sprintf("%.6e", float(x)) : "nan"
end

function _csv_escape(s)
    t = string(s)
    if occursin(',', t) || occursin('"', t)
        return "\"" * replace(t, '"' => "\"\"") * "\""
    end
    return t
end

function _push_num!(rows, q, v, src)
    push!(rows, (q, _fmt(v), src))
    return nothing
end

function _detuning(c, J)
    return ntuple(3) do k
        acc = 2 * c[k, k] * J[k]
        for j in 1:3
            j == k && continue
            acc += c[min(k, j), max(k, j)] * J[j]
        end
        acc
    end
end

function _frac_distance(a, b)
    d = mod(float(a) - float(b), 1.0)
    return min(d, 1.0 - d)
end

function _write_numbers(path, mag, rf, kicks, nf_result, track)
    rows = Tuple{String,String,String}[]
    pe = mag.pe
    _push_num!(rows, "C_m", mag.C, "path_length")
    _push_num!(rows, "gamma", mag.ref.gamma, "reference")
    _push_num!(rows, "alpha0", pe.alpha0, "path_expansion")
    _push_num!(rows, "alpha1", pe.alpha1, "path_expansion")
    _push_num!(rows, "eta0", pe.eta0, "path_expansion")
    _push_num!(rows, "eta1", pe.eta1, "path_expansion")
    _push_num!(rows, "eta1_formula", mag.eta1_formula, "alpha1 - alpha0/gamma^2 + (3*gamma^2-1)/(2*gamma^4)")
    _push_num!(rows, "eta1_minus_eta1_formula", mag.eta1_minus_formula, "path_expansion minus the note formula")
    _push_num!(rows, "bell_eta1_formula", mag.bell_eta1_formula, "path_expansion eta1_formula")
    _push_num!(rows, "eta1_minus_bell_eta1_formula", mag.eta1_minus_bell, "path_expansion identity")
    _push_num!(rows, "note_formula_minus_bell_formula", mag.formula_minus_bell, "same expression when the velocity piece matches")
    _push_num!(rows, "velocity_piece_minus_note_piece", mag.velocity_piece_minus_note, "-c2 + c1^2 minus (3*gamma^2-1)/(2*gamma^4)")
    _push_num!(rows, "matrix_gap", mag.matrix_gap, "max |Lie R - ring_matrix|")
    for (plane, sym) in ((:x, "x"), (:y, "y"))
        for which in (:xi1, :xi2)
            g = _xi(mag.gen, plane, which)
            fine = _orders_xi(mag.orders.fine, plane, which)
            coarse = _orders_xi(mag.orders.coarse, plane, which)
            tag = string(which)
            _push_num!(rows, "$(tag)_$(sym)_generator", g, "generator_chromaticity, factor -1/(2*pi)")
            _push_num!(rows, "$(tag)_$(sym)_orders_1e-6", fine, "chromaticity_orders")
            _push_num!(rows, "$(tag)_$(sym)_orders_1e-4", coarse, "chromaticity_orders")
            _push_num!(rows, "$(tag)_$(sym)_generator_minus_orders_1e-6", g - fine, "generator minus orders")
        end
        c1 = plane == :x ? mag.gen.x1.c : mag.gen.y1.c
        c2 = plane == :x ? mag.gen.x2.c : mag.gen.y2.c
        _push_num!(rows, "c_Jdelta_$(sym)", c1, "coefficient of J delta")
        _push_num!(rows, "c_Jdelta2_$(sym)", c2, "coefficient of J delta^2 in f4 - 1/2[f3, F3]")
    end
    _push_num!(rows, "voltage_MV", rf.cav.voltage_MV, "regimes.linear")
    _push_num!(rows, "harmonic", rf.cav.harmonic, "regimes.linear")
    _push_num!(rows, "phi_s_bell_rad", rf.cav.phi_s, "RFCavity")
    _push_num!(rows, "phi_s_notes_rad", rf.regime.phi_s_rad, "regimes.linear in the regime file")
    _push_num!(rows, "Qs_note", rf.Qs_note, "synchrotron_tune_note")
    _push_num!(rows, "Qs_smooth", rf.Qs_smooth, "synchrotron_tune at phi_s = pi")
    _push_num!(rows, "Qs_json", rf.Qs_json, "regimes.linear.Qs")
    _push_num!(rows, "kick_R_gap", kicks.R_gap, "ring_lie_factor(with cavity) versus concatenate(magnetic, one kick)")
    _push_num!(rows, "kick_f3_gap", kicks.f3_gap, "same comparison, f3")
    _push_num!(rows, "kick_f4_gap", kicks.f4_gap, "same comparison, f4")
    _push_num!(rows, "Qs_one_kick", kicks.Qs_one, "|mu_s|/(2*pi) of the cavity lattice")
    _push_num!(rows, "Qs_second_kick", kicks.Qs_two, "|mu_s|/(2*pi) after one more rf_kick_taylor")
    _push_num!(rows, "Qs_second_over_one", kicks.Qs_ratio, "second kick over one kick")
    _push_num!(rows, "nf_ok", nf_result.ok ? 1.0 : 0.0, "1 if dragt_finn_nf6d returned")
    if nf_result.ok
        nf = nf_result.nf
        _push_num!(rows, "mu_x_rad", nf.mu[1], "dragt_finn_nf6d")
        _push_num!(rows, "mu_y_rad", nf.mu[2], "dragt_finn_nf6d")
        _push_num!(rows, "mu_s_rad", nf.mu[3], "dragt_finn_nf6d")
        _push_num!(rows, "Qs_nf", nf_result.Qs, "|mu_s|/(2*pi)")
        _push_num!(rows, "Qs_nf_minus_Qs_note", nf_result.Qs - rf.Qs_note, "normal form minus the note tune")
        _push_num!(rows, "Qx_nf_frac_distance", _frac_distance(nf.mu[1] / (2π), mag.po.x.mu / (2π)),
                   "wrapped distance of mu_x/2pi to the magnetic fractional tune")
        _push_num!(rows, "Qy_nf_frac_distance", _frac_distance(nf.mu[2] / (2π), mag.po.y.mu / (2π)),
                   "wrapped distance of mu_y/2pi to the magnetic fractional tune")
        _push_num!(rows, "Rn_residual", nf.Rn_residual, "|A^-1 R A - rotation|")
        _push_num!(rows, "symp_residual", nf.symp_residual, "|A^T S A - S|")
        _push_num!(rows, "f3_before", nf.f3_before, "max |coef| of f3")
        _push_num!(rows, "f3_after", nf.f3_after, "cubic homological residual")
        _push_num!(rows, "f4_nonres_after", nf.f4_nonres_after, "max |f4 - h4| after F4")
        for i in 1:3, j in 1:3
            _push_num!(rows, "c_$(i)$(j)", nf.c[i, j], "h4 = sum_{i<=j} c_ij J_i J_j; lower triangle stored as 0")
        end
        for k in 1:3
            _push_num!(rows, "d_$(k)", nf.d[k], "delta_eq = zstar_6 + d · J")
        end
        _push_num!(rows, "sigma_delta_s", nf.sigma_delta_s, "delta amplitude per sqrt(2 Js)")
        _push_num!(rows, "sigma_l_s", nf.sigma_l_s, "l amplitude per sqrt(2 Js)")
        if track !== nothing
            J0 = (track.Jx_nf[1], track.Jy_nf[1], track.Js_nf[1])
            g = _detuning(nf.c, J0)
            for k in 1:3
                _push_num!(rows, "dh4_dJ_$(k)", g[k], "partial h4 / partial J at the initial actions")
                _push_num!(rows, "dQ_$(k)", g[k] / (2π), "tune shift -dQ from that gradient, in tune units")
            end
        end
    end
    if track !== nothing
        _push_num!(rows, "track_turns", track.nturns, "exact_track_element and rf_kick")
        _push_num!(rows, "track_ds_m", track.ds, "step passed to exact_track_element")
        _push_num!(rows, "ax_m", track.ax, "horizontal betatron offset")
        _push_num!(rows, "ay_m", track.ay, "vertical betatron offset")
        _push_num!(rows, "delta_m", rf.regime.delta_m, "regimes.linear")
        for (name, st) in (
                ("Jx_nf", track.stats.Jx_nf), ("Jy_nf", track.stats.Jy_nf), ("Js_nf", track.stats.Js_nf),
                ("Jx_cs", track.stats.Jx_cs), ("Jy_cs", track.stats.Jy_cs),
                ("Jx_cs_dispersion_free", track.stats.Jx_beta), ("Jy_cs_dispersion_free", track.stats.Jy_beta),
            )
            _push_num!(rows, "$(name)_mean", st.mean, "turn 0 through the last turn")
            _push_num!(rows, "$(name)_span", st.span, "max minus min")
            _push_num!(rows, "$(name)_span_over_abs_mean", st.span_over_abs_mean, "(max-min)/|mean|")
            _push_num!(rows, "$(name)_rms", st.rms, "root mean square about the mean")
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

function _write_actions(path, track)
    open(path, "w") do io
        println(io, "turn,x,px,y,py,l,delta,Jx_nf,Jy_nf,Js_nf,Jx_cs,Jy_cs,Jx_cs_dispersion_free,Jy_cs_dispersion_free")
        for i in eachindex(track.turn)
            z = track.z[:, i]
            @printf(io, "%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                    track.turn[i], z[1], z[2], z[3], z[4], z[5], z[6],
                    track.Jx_nf[i], track.Jy_nf[i], track.Js_nf[i],
                    track.Jx_cs[i], track.Jy_cs[i], track.Jx_beta[i], track.Jy_beta[i])
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

function _limits(v)
    lo = minimum(v)
    hi = maximum(v)
    span = hi - lo
    if span < 1.0e-30
        span = max(1.0e-18, abs(lo) * 0.2)
        mid = (lo + hi) / 2
        return mid - span / 2, mid + span / 2
    end
    pad = 0.08 * span
    return lo - pad, hi + pad
end

function _tick_label(x, span)
    ax = abs(x)
    if span == 0
        return @sprintf("%.1e", x)
    elseif span >= 20
        return @sprintf("%.0f", x)
    elseif ax >= 0.01 && span >= 0.01
        return @sprintf("%.4g", x)
    else
        return @sprintf("%.2e", x)
    end
end

function _polyline(xs, ys, xp, yp)
    parts = String[]
    for i in eachindex(xs)
        push!(parts, @sprintf("%.2f,%.2f", xp(xs[i]), yp(ys[i])))
    end
    return join(parts, " ")
end

function _xml(s)
    return replace(replace(replace(string(s), "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")
end

function _write_svg(path, track)
    W = 980.0
    left, right, top = 88.0, 24.0, 36.0
    pw = W - left - right
    ph, gap = 150.0, 28.0
    panels = (
        (title = "Jx", ylabel = "Jx",
         curves = (("NF", track.Jx_nf, "#1d4e89", ""),
                   ("CS", track.Jx_cs, "#c0392b", "6 4"),
                   ("CS, D removed", track.Jx_beta, "#1e7f4f", "2 3"))),
        (title = "Jy", ylabel = "Jy",
         curves = (("NF", track.Jy_nf, "#1d4e89", ""),
                   ("CS", track.Jy_cs, "#c0392b", "6 4"),
                   ("CS, D removed", track.Jy_beta, "#1e7f4f", "2 3"))),
        (title = "Js", ylabel = "Js",
         curves = (("NF", track.Js_nf, "#1d4e89", ""),)),
        (title = "delta", ylabel = "delta",
         curves = (("delta", track.z[6, :], "#7a4e12", ""),)),
    )
    H = top + length(panels) * (ph + gap) + 16
    turns = Float64.(track.turn)
    tmax = float(track.nturns)
    open(path, "w") do io
        println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
        @printf(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" viewBox=\"0 0 %.0f %.0f\">\n",
                W, H, W, H)
        println(io, "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>")
        println(io, "<text x=\"88\" y=\"22\" font-family=\"sans-serif\" font-size=\"16\" fill=\"#222\">normal-form actions and Courant-Snyder actions</text>")
        for (p, panel) in enumerate(panels)
            y0 = top + (p - 1) * (ph + gap)
            ys = Float64[]
            for c in panel.curves
                append!(ys, c[2])
            end
            ylo, yhi = _limits(ys)
            span = yhi - ylo
            xp = x -> left + (tmax == 0 ? 0.0 : (x / tmax) * pw)
            yp = y -> y0 + (yhi - y) / (yhi - ylo) * ph
            @printf(io, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"#fcfcfc\" stroke=\"#dddddd\"/>\n",
                    left, y0, pw, ph)
            println(io, "<text x=\"16\" y=\"", y0 + ph / 2,
                    "\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#333\" transform=\"rotate(-90 16,",
                    y0 + ph / 2, ")\">", _xml(panel.ylabel), "</text>")
            for tick in 0:4
                xv = tick * tmax / 4
                px = xp(xv)
                @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#eeeeee\"/>\n",
                        px, y0, px, y0 + ph)
                @printf(io, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#444\">%s</text>\n",
                        px, y0 + ph + 14, _tick_label(xv, max(tmax, 1.0)))
            end
            for tick in 0:3
                yv = ylo + tick * span / 3
                py = yp(yv)
                @printf(io, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#444\">%s</text>\n",
                        left - 6, py + 4, _xml(_tick_label(yv, span)))
            end
            for (name, y, color, dash) in panel.curves
                dashattr = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
                println(io, "<polyline fill=\"none\" stroke=\"", color, "\" stroke-width=\"1.4\"",
                        dashattr, " points=\"", _polyline(turns, y, xp, yp), "\"/>")
            end
            for (j, (name, _, color, dash)) in enumerate(panel.curves)
                lx = left + 12
                ly = y0 + 16 + (j - 1) * 14
                dashattr = isempty(dash) ? "" : " stroke-dasharray=\"$dash\""
                @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-width=\"2\"%s/>\n",
                        lx, ly - 3, lx + 18, ly - 3, color, dashattr)
                @printf(io, "<text x=\"%.2f\" y=\"%.2f\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#222\">%s</text>\n",
                        lx + 22, ly, _xml(name))
            end
        end
        println(io, "</svg>")
    end
    return path
end

function _chrom_rows(mag)
    out = Vector{Vector{String}}()
    for (plane, sym) in ((:x, "x"), (:y, "y"))
        for which in (:xi1, :xi2)
            g = _xi(mag.gen, plane, which)
            fine = _orders_xi(mag.orders.fine, plane, which)
            coarse = _orders_xi(mag.orders.coarse, plane, which)
            push!(out, [sym, string(which), _fmt(g), _fmt(coarse), _fmt(fine), _fmt(g - fine)])
        end
    end
    return out
end

function _action_rows(track)
    rows = Vector{Vector{String}}()
    for (name, st) in (
            ("Jx NF", track.stats.Jx_nf),
            ("Jy NF", track.stats.Jy_nf),
            ("Js NF", track.stats.Js_nf),
            ("Jx Courant-Snyder", track.stats.Jx_cs),
            ("Jy Courant-Snyder", track.stats.Jy_cs),
            ("Jx Courant-Snyder, dispersion removed", track.stats.Jx_beta),
            ("Jy Courant-Snyder, dispersion removed", track.stats.Jy_beta),
        )
        push!(rows, [name, _fmt(st.mean), _fmt(st.minimum), _fmt(st.maximum),
                     _fmt(st.span), _fmt6(st.span_over_abs_mean)])
    end
    return rows
end

function _xi_agree(mag)
    flags = String[]
    for (plane, sym) in ((:x, "x"), (:y, "y"))
        for (which, tol) in ((:xi1, XI1_TOL), (:xi2, XI2_TOL))
            d = abs(_xi(mag.gen, plane, which) - _orders_xi(mag.orders.fine, plane, which))
            if d > tol
                push!(flags, "$sym $which |generator − orders| = $d, tolerance $tol")
            end
        end
    end
    return (agree = isempty(flags), flags = flags)
end

function _write_report(path, latpath, mag, rf, kicks, nf_result, track, track_error)
    pe = mag.pe
    ref = mag.ref
    xi = _xi_agree(mag)
    lines = String[]
    push!(lines, "# Lie invariants of the closed magnetic ring")
    push!(lines, "")
    push!(lines, "Lattice file: `$latpath`. Sextupoles are at `k2 = 0`. Reference: $(ref.particle.name), kinetic energy $(ref.Ekin_MeV) MeV, gamma $( _fmt(ref.gamma) ), harmonic $(rf.lattice.harmonic).")
    push!(lines, "")
    push!(lines, "The magnetic factor is `ring_lie_factor` of `magnetic_lattice`. Slip coefficients are `path_expansion` of `factor_to_map`. Chromaticity in the generator column is `generator_chromaticity` with the factor fixed in step 2b, `ξ1 = −c / (2π)` and `ξ2 = 2 ξ1_factor c₂`, placed next to `chromaticity_orders`.")
    push!(lines, "")
    push!(lines, "## Slip")
    push!(lines, "")
    push!(lines, "The note formula is `η1 = α1 − α0/γ² + (3γ² − 1)/(2γ⁴)`. `eta1 − eta1_formula` uses that expression. `path_expansion` also returns its own `eta1_formula`.")
    push!(lines, "")
    push!(lines, _md_table(["quantity", "value"], [
        ["C [m]", _fmt(mag.C)],
        ["gamma", _fmt(ref.gamma)],
        ["1/gamma^2", _fmt(1 / ref.gamma^2)],
        ["alpha0", _fmt(pe.alpha0)],
        ["alpha1", _fmt(pe.alpha1)],
        ["eta0", _fmt(pe.eta0)],
        ["eta1", _fmt(pe.eta1)],
        ["eta1_formula", _fmt(mag.eta1_formula)],
        ["eta1 − eta1_formula", _fmt6(mag.eta1_minus_formula)],
        ["path_expansion eta1_formula", _fmt(mag.bell_eta1_formula)],
        ["eta1 − path_expansion eta1_formula", _fmt6(mag.eta1_minus_bell)],
        ["note formula − path_expansion formula", _fmt6(mag.formula_minus_bell)],
        ["(-c2 + c1^2) − (3γ²−1)/(2γ⁴)", _fmt6(mag.velocity_piece_minus_note)],
        ["|Lie R − ring_matrix|", _fmt6(mag.matrix_gap)],
    ]))
    push!(lines, "")
    if abs(mag.eta1_minus_formula) <= SLIP_TOL && abs(mag.formula_minus_bell) <= SLIP_TOL
        push!(lines, "η1 and the note formula agree within $SLIP_TOL. The note formula and `path_expansion`'s `eta1_formula` agree within the same tolerance. The velocity identity `−c2 + c1² = (3γ² − 1)/(2γ⁴)` is the algebraic reason those two formulas meet.")
    else
        push!(lines, "The slip readings differ by more than $SLIP_TOL. Both values are listed above. This step does not choose one of them.")
    end
    push!(lines, "")
    push!(lines, "## Chromaticity")
    push!(lines, "")
    push!(lines, "Generator ξ1 is `−1/(2π)` times the coefficient of `J δ`. Generator ξ2 is `−1/π` times the coefficient of `J δ²` in `f4 − ½[f3, F3]`. `chromaticity_orders` is the central difference of the same magnetic map. The step `1e-6` is the settled difference. The step `1e-4` is the default of that function.")
    push!(lines, "")
    push!(lines, _md_table(["plane", "quantity", "generator", "orders 1e-4", "orders 1e-6", "generator − orders 1e-6"],
                           _chrom_rows(mag)))
    push!(lines, "")
    push!(lines, _md_table(["plane", "c of J δ", "c of J δ²"], [
        ["x", _fmt(mag.gen.x1.c), _fmt(mag.gen.x2.c)],
        ["y", _fmt(mag.gen.y1.c), _fmt(mag.gen.y2.c)],
    ]))
    push!(lines, "")
    if xi.agree
        push!(lines, "Generator and `chromaticity_orders` at step `1e-6` agree within $XI1_TOL on ξ1 and within $XI2_TOL on ξ2.")
    else
        push!(lines, "Generator and `chromaticity_orders` differ. Both columns stay in the table.")
        push!(lines, "")
        for line in xi.flags
            push!(lines, "- $line")
        end
    end
    push!(lines, "")
    push!(lines, "## Thin cavity")
    push!(lines, "")
    push!(lines, "One thin `RFCavity` is appended to the element list. Voltage and harmonic are `regimes.linear` in `02_rf_regimes.json`. The cavity phase passed to BELL is `$( _fmt(BELL_PHI_S) )` rad. The regime file records the notes' stationary phase `$( _fmt(rf.regime.phi_s_rad) )` rad. `eta0` is negative, and `eta0 * cos(phi_s)` is positive at the BELL phase, so the smooth tune is real.")
    push!(lines, "")
    push!(lines, "The one-turn factor used below is `ring_lie_factor` of this lattice. The cavity is already an element, so its kick is inside that product. `one_turn_factor_rf` is not called.")
    push!(lines, "")
    push!(lines, _md_table(["quantity", "value"], [
        ["h", string(rf.cav.harmonic)],
        ["V [MV]", _fmt(rf.cav.voltage_MV)],
        ["phi_s on the cavity [rad]", _fmt(rf.cav.phi_s)],
        ["phi_s in the regime file [rad]", _fmt(rf.regime.phi_s_rad)],
        ["Q_s note", _fmt(rf.Qs_note)],
        ["Q_s smooth model at phi_s = π", _fmt(rf.Qs_smooth)],
        ["Q_s stored in the regime file", _fmt(rf.Qs_json)],
        ["Q_s note − regime file", _fmt6(rf.Qs_note - rf.Qs_json)],
        ["Q_s note − smooth model", _fmt6(rf.Qs_note - rf.Qs_smooth)],
    ]))
    push!(lines, "")
    push!(lines, "## Two ways to attach one kick, and a second kick")
    push!(lines, "")
    push!(lines, "`ring_lie_factor` of the lattice with the cavity is compared with `concatenate` of the magnetic factor and `dragt_finn_factor(rf_kick_taylor(...))`. A further `concatenate` of the cavity-lattice factor with that same kick is the map `one_turn_factor_rf` builds when the lattice already contains the cavity: `ring_lie_factor` includes the cavity, and the function appends `rf_kick_taylor` again. Synchrotron tunes in this section are `|μ_s| / 2π` from `symplectic_eigen_A` of those linear matrices.")
    push!(lines, "")
    push!(lines, _md_table(["quantity", "value"], [
        ["|R(cavity lattice) − R(magnetic + one kick)|", _fmt6(kicks.R_gap)],
        ["|f3| gap of that comparison", _fmt6(kicks.f3_gap)],
        ["|f4| gap of that comparison", _fmt6(kicks.f4_gap)],
        ["Q_s of the cavity lattice", _fmt(kicks.Qs_one)],
        ["Q_s after a second kick", _fmt(kicks.Qs_two)],
        ["Q_s second / Q_s one", _fmt(kicks.Qs_ratio)],
    ]))
    push!(lines, "")
    if !kicks.mu_one.ok
        push!(lines, "Eigenvalues of the cavity-lattice matrix: $(kicks.mu_one.error)")
        push!(lines, "")
    end
    if !kicks.mu_two.ok
        push!(lines, "Eigenvalues after the second kick: $(kicks.mu_two.error)")
        push!(lines, "")
    end
    push!(lines, "## Normal form")
    push!(lines, "")
    if !nf_result.ok
        push!(lines, "`dragt_finn_nf6d` threw. The magnetic tables above are the delivered reading. Tracking was not run.")
        push!(lines, "")
        push!(lines, "```")
        push!(lines, nf_result.error)
        push!(lines, "```")
        push!(lines, "")
    else
        nf = nf_result.nf
        push!(lines, "`dragt_finn_nf6d` returned. `Q_s = |μ_s| / 2π`. The detuning matrix is stored as `h4 = Σ_{i≤j} c_ij J_i J_j`, with zeros below the diagonal. `μ_k(J) = μ_k − ∂h4/∂J_k`, and `∂h4/∂J_k = 2 c_kk J_k + Σ_{j≠k} c_{min,max} J_j`.")
        push!(lines, "")
        push!(lines, _md_table(["quantity", "value"], [
            ["mu_x [rad]", _fmt(nf.mu[1])],
            ["mu_y [rad]", _fmt(nf.mu[2])],
            ["mu_s [rad]", _fmt(nf.mu[3])],
            ["Q_s", _fmt(nf_result.Qs)],
            ["Q_s − Q_s note", _fmt6(nf_result.Qs - rf.Qs_note)],
            ["wrapped |Q_x − Q_x magnetic|", _fmt6(_frac_distance(nf.mu[1] / (2π), mag.po.x.mu / (2π)))],
            ["wrapped |Q_y − Q_y magnetic|", _fmt6(_frac_distance(nf.mu[2] / (2π), mag.po.y.mu / (2π)))],
            ["|A^-1 R A − rotation|", _fmt6(nf.Rn_residual)],
            ["|A^T S A − S|", _fmt6(nf.symp_residual)],
            ["max |coef| of f3", _fmt6(nf.f3_before)],
            ["cubic residual after F3", _fmt6(nf.f3_after)],
            ["|f4 − h4| after F4", _fmt6(nf.f4_nonres_after)],
            ["sigma_delta per sqrt(2 Js)", _fmt(nf.sigma_delta_s)],
            ["sigma_l per sqrt(2 Js)", _fmt(nf.sigma_l_s)],
        ]))
        push!(lines, "")
        c_rows = Vector{Vector{String}}()
        for i in 1:3
            push!(c_rows, [_fmt(nf.c[i, 1]), _fmt(nf.c[i, 2]), _fmt(nf.c[i, 3])])
        end
        push!(lines, "Detuning matrix `c`, rows and columns `(x, y, s)`:")
        push!(lines, "")
        push!(lines, _md_table(["c_i1", "c_i2", "c_i3"], c_rows))
        push!(lines, "")
        push!(lines, _md_table(["plane", "d_k"], [
            ["x", _fmt(nf.d[1])],
            ["y", _fmt(nf.d[2])],
            ["s", _fmt(nf.d[3])],
        ]))
        push!(lines, "")
        push!(lines, "`d` is the normal-form coefficient in `δ_eq = z*_6 + d · J`. `z*` from this call is the origin.")
        push!(lines, "")
    end
    push!(lines, "## Tracking")
    push!(lines, "")
    if track === nothing
        if !nf_result.ok
            push!(lines, "No turn-by-turn table. The normal form did not return, so the normal-form actions were not evaluated.")
        elseif isempty(track_error)
            push!(lines, "Tracking has not been written yet.")
        else
            push!(lines, "Tracking threw. The magnetic tables and the normal form above are kept.")
            push!(lines, "")
            push!(lines, "```")
            push!(lines, track_error)
            push!(lines, "```")
        end
        push!(lines, "")
    else
        push!(lines, "Each magnet is `exact_track_element` with step $(track.ds) m. The thin cavity is `rf_kick`. The observation point is the start of the element list, which is the entrance of the one-turn map. The launch point is the periodic dispersion orbit at `δ = δ_m` plus a betatron offset. The horizontal offset is twice `|D_x δ_m|` at this location (`a_x = $( _fmt(track.ax) )` m, `a_y = $( _fmt(track.ay) )` m, `δ_m = $( _fmt(rf.regime.delta_m) )`). `l = 0`, so the particle starts at the momentum crest of the stationary bucket.")
        push!(lines, "")
        push!(lines, "Normal-form actions use `u = A⁻¹ (z − z*)` and `w = Tinv(u)`, then `J_k = (w_{2k−1}² + w_{2k}²) / 2`. Courant–Snyder `J` uses the on-momentum Twiss at this location. The geometric Courant–Snyder action keeps the dispersive piece of `(x, p_x)`. The dispersion-free column subtracts `D δ` first.")
        push!(lines, "")
        push!(lines, @sprintf("%d turns, including the launch point as turn 0. Elapsed time %.1f s. Synchrotron periods at the normal-form Q_s: %.3f.",
                              track.nturns, track.elapsed, track.nturns * nf_result.Qs))
        push!(lines, "")
        push!(lines, _md_table(["action", "mean", "min", "max", "max − min", "(max−min)/|mean|"],
                               _action_rows(track)))
        push!(lines, "")
        push!(lines, "Figure: `03_actions.svg`. Samples: `03_actions.csv`.")
        push!(lines, "")
    end
    write(path, join(lines, "\n"))
    return path
end

function _write_draft(path, kicks, nf_result, rf)
    lines = String[]
    push!(lines, "# one_turn_factor_rf appends a kick the lattice may already contain")
    push!(lines, "")
    push!(lines, "## Summary")
    push!(lines, "")
    push!(lines, "`one_turn_factor_rf(lattice, cav, ref, path_length)` calls `ring_lie_factor(lattice)` and then concatenates `dragt_finn_factor(rf_kick_taylor(cav, ...))`. `ring_lie_factor` builds a factor for every element, and a thin `RFCavity` is the Lie map of `rf_potential`. On a lattice that already contains that cavity, the function applies the voltage twice.")
    push!(lines, "")
    push!(lines, "## What this ring shows")
    push!(lines, "")
    push!(lines, "Closed `magnetic_2`, sextupoles at `k2 = 0`, one thin cavity appended at the end. Voltage $( _fmt(rf.cav.voltage_MV) ) MV, harmonic $(rf.cav.harmonic), `phi_s = π`. The study factor is `ring_lie_factor` of that lattice. `one_turn_factor_rf` was not called.")
    push!(lines, "")
    push!(lines, "Replacing the in-lattice cavity factor by one `rf_kick_taylor` factor changes the one-turn matrix by $( _fmt6(kicks.R_gap) ), f3 by $( _fmt6(kicks.f3_gap) ), and f4 by $( _fmt6(kicks.f4_gap) ). Concatenating that kick once more, which is the body of `one_turn_factor_rf` on this lattice, moves `|μ_s| / 2π` from $( _fmt(kicks.Qs_one) ) to $( _fmt(kicks.Qs_two) ). The ratio is $( _fmt(kicks.Qs_ratio) ).")
    push!(lines, "")
    if !kicks.mu_one.ok || !kicks.mu_two.ok
        push!(lines, "Eigenvalue status of the one-kick matrix: $(kicks.mu_one.ok ? "returned" : kicks.mu_one.error). Second kick: $(kicks.mu_two.ok ? "returned" : kicks.mu_two.error).")
        push!(lines, "")
    end
    if nf_result.ok
        push!(lines, "The normal form of the single in-lattice cavity returned `Q_s = $( _fmt(nf_result.Qs) )`.")
        push!(lines, "")
    else
        push!(lines, "`dragt_finn_nf6d` on the single-cavity factor threw: $(nf_result.error)")
        push!(lines, "")
    end
    push!(lines, "## Recommendation")
    push!(lines, "")
    push!(lines, "`one_turn_factor_rf` should take the magnetic lattice, via `magnetic_lattice`, before `ring_lie_factor`, or it should refuse a lattice that already contains an `RFCavity`. The docstring should state that a thin cavity already in the element list is a kick, so a second `rf_kick_taylor` doubles the voltage. Callers that have already appended the cavity should use `ring_lie_factor` and should not call `one_turn_factor_rf`.")
    push!(lines, "")
    push!(lines, "## Reproducer")
    push!(lines, "")
    push!(lines, "```julia")
    push!(lines, "using Bell")
    push!(lines, "using Bell.NormalForm: one_turn_factor_rf, symplectic_eigen_A")
    push!(lines, "lat = load_structure(raw\"$(replace(closed_lattice_path(), "\\" => "/"))\")")
    push!(lines, "ref = lat.ref")
    push!(lines, "C = path_length(lat)")
    push!(lines, "cav = RFCavity(\"RF\", $(rf.cav.voltage_MV), $(rf.cav.harmonic), π)")
    push!(lines, "elems = AbstractElement[e for e in lat.elements]")
    push!(lines, "push!(elems, cav)")
    push!(lines, "lat_rf = Lattice(lat.name * \"+RF\", elems, ref, lat.harmonic)")
    push!(lines, "F = ring_lie_factor(lat_rf)                          # one kick")
    push!(lines, "F2 = one_turn_factor_rf(lat_rf, cav, ref, C)         # two kicks")
    push!(lines, "_, mu = symplectic_eigen_A(F.R)")
    push!(lines, "_, mu2 = symplectic_eigen_A(F2.R)")
    push!(lines, "abs(mu[3]) / 2π, abs(mu2[3]) / 2π")
    push!(lines, "```")
    push!(lines, "")
    push!(lines, "This note records the double count. It does not change BELL.")
    push!(lines, "")
    write(path, join(lines, "\n"))
    return path
end

function _print_summary(mag, rf, kicks, nf_result, track, paths, track_error)
    println()
    @printf("eta1 = %.8e\n", mag.pe.eta1)
    @printf("eta1_formula = %.8e\n", mag.eta1_formula)
    @printf("eta1 - eta1_formula = %.6e\n", mag.eta1_minus_formula)
    @printf("xi1 generator x, y = %.8e  %.8e\n", _xi(mag.gen, :x, :xi1), _xi(mag.gen, :y, :xi1))
    @printf("xi1 orders    x, y = %.8e  %.8e\n",
            _orders_xi(mag.orders.fine, :x, :xi1), _orders_xi(mag.orders.fine, :y, :xi1))
    @printf("xi2 generator x, y = %.8e  %.8e\n", _xi(mag.gen, :x, :xi2), _xi(mag.gen, :y, :xi2))
    @printf("xi2 orders    x, y = %.8e  %.8e\n",
            _orders_xi(mag.orders.fine, :x, :xi2), _orders_xi(mag.orders.fine, :y, :xi2))
    @printf("V = %.8e MV   h = %d   phi_s = %.6f rad   Qs_note = %.8e\n",
            rf.cav.voltage_MV, rf.cav.harmonic, rf.cav.phi_s, rf.Qs_note)
    @printf("kick |R| gap = %.3e    Qs one = %.8e    Qs second kick = %.8e    ratio = %.6f\n",
            kicks.R_gap, kicks.Qs_one, kicks.Qs_two, kicks.Qs_ratio)
    if nf_result.ok
        nf = nf_result.nf
        println("NF returned")
        @printf("mu = %.8e  %.8e  %.8e\n", nf.mu[1], nf.mu[2], nf.mu[3])
        @printf("Qs = %.8e    Rn residual = %.3e    symp residual = %.3e    |f4-h4| = %.3e\n",
                nf_result.Qs, nf.Rn_residual, nf.symp_residual, nf.f4_nonres_after)
        @printf("c diagonal = %.8e  %.8e  %.8e\n", nf.c[1, 1], nf.c[2, 2], nf.c[3, 3])
    else
        println("NF failed: ", nf_result.error)
    end
    if track !== nothing
        st = track.stats
        @printf("turns = %d\n", track.nturns)
        @printf("Jx NF  mean %.8e  span %.8e  span/|mean| %.6e\n",
                st.Jx_nf.mean, st.Jx_nf.span, st.Jx_nf.span_over_abs_mean)
        @printf("Jy NF  mean %.8e  span %.8e  span/|mean| %.6e\n",
                st.Jy_nf.mean, st.Jy_nf.span, st.Jy_nf.span_over_abs_mean)
        @printf("Js NF  mean %.8e  span %.8e  span/|mean| %.6e\n",
                st.Js_nf.mean, st.Js_nf.span, st.Js_nf.span_over_abs_mean)
        @printf("Jx CS  mean %.8e  span %.8e  span/|mean| %.6e\n",
                st.Jx_cs.mean, st.Jx_cs.span, st.Jx_cs.span_over_abs_mean)
        @printf("Jy CS  mean %.8e  span %.8e  span/|mean| %.6e\n",
                st.Jy_cs.mean, st.Jy_cs.span, st.Jy_cs.span_over_abs_mean)
    end
    println("outputs:")
    for p in paths
        println("  ", p)
    end
    flush(stdout)
    return nothing
end

main()
