# Step 2. Linear optics of the closed magnetic_2 ring, from the Lie factor
# and from the one-turn matrix. Choose two stationary RF regimes.
#
#   julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/02_linear_optics.jl

import Pkg

const STUDY_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(STUDY_DIR)

using Bell
using MagneticRing
using Printf

const DS = 0.02
const OUT_DIR = joinpath(STUDY_DIR, "output")
const TOL = 1.0e-8

function main()
    path = closed_lattice_path()
    println("Loading ", path)
    flush(stdout)
    lattice = load_structure(path)
    _require_closed_lattice(lattice)

    ref = lattice.ref
    C = path_length(lattice)
    println("One-turn matrix")
    flush(stdout)
    R = ring_matrix(lattice)
    po = periodic_optics(R, ref, C)
    po.x.stable && po.y.stable || error("periodic optics of the one-turn matrix is unstable")

    println("One-turn Lie factor and path expansion")
    flush(stdout)
    obs = observables(lattice)
    pe = path_expansion(obs.M, ref, C)
    _readings_agree(R, po, obs, pe, C)

    println("Optics along s, ds = ", DS, " m")
    flush(stdout)
    opt = optics_along(lattice, po; ds = DS)
    disp = sample_dispersion(lattice, po; ds = DS)
    _transport_agrees(opt, disp, po)
    alpha_geom = _compaction_integral(opt)
    d_alpha = _must_agree("alpha0 versus (1/C) int D/rho ds", pe.alpha0, alpha_geom, 1.0e-4)

    cells = _arc_cells(lattice)
    length(cells) == 8 || error("expected 8 arc cells, found $(length(cells))")
    _bends_inside_cells(lattice, cells)
    cell_rows = _cell_table(lattice, opt, cells)
    gap_rows = _gap_table(lattice, opt, cell_rows)
    _phases_cover_the_ring(opt, cell_rows, gap_rows)
    half = _superperiod(lattice, opt)

    api = find_optics_maxima(opt)
    peaks = _curve_peaks(opt, disp, api)
    println("Curve maxima: ", length(peaks), "    find_optics_maxima: ", length(api))
    flush(stdout)

    xi_x = obs.xi.x[2]
    xi_y = obs.xi.y[2]
    regimes = build_rf_regimes(ref, C, pe.eta0, pe.eta1, opt.Qx, opt.Qy,
                               maximum(opt.data.betx), maximum(opt.data.bety),
                               xi_x, xi_y, lattice.harmonic)
    qs_smooth = _smooth_model_tunes(ref, pe.eta0, regimes)

    mkpath(OUT_DIR)
    paths = (
        report = joinpath(OUT_DIR, "02_linear_optics.md"),
        numbers = joinpath(OUT_DIR, "02_numbers.csv"),
        cells = joinpath(OUT_DIR, "02_cells.csv"),
        maxima = joinpath(OUT_DIR, "02_maxima.csv"),
        optics = joinpath(OUT_DIR, "02_optics.csv"),
        figure = joinpath(OUT_DIR, "02_optics.svg"),
        rf_json = joinpath(OUT_DIR, "02_rf_regimes.json"),
        rf_md = joinpath(OUT_DIR, "02_rf_regimes.md"),
    )
    _write_optics_csv(paths.optics, opt, disp)
    _write_numbers_csv(paths.numbers, opt, po, pe, C, cell_rows, half, alpha_geom)
    _write_cells_csv(paths.cells, cell_rows, gap_rows)
    _write_maxima_csv(paths.maxima, peaks)
    _write_svg(paths.figure, opt, disp, peaks)
    _write_rf_json(paths.rf_json, ref, C, opt, po, pe, cell_rows, half, regimes, qs_smooth, alpha_geom, xi_x, xi_y)
    _write_rf_md(paths.rf_md, ref, pe, regimes, qs_smooth)
    _write_report(paths.report, lattice, path, opt, po, pe, obs, cell_rows, gap_rows, half,
                  peaks, api, regimes, qs_smooth, alpha_geom, d_alpha, disp)

    _print_summary(opt, pe, cell_rows, regimes, paths)
    return nothing
end

function _require_closed_lattice(lattice)
    ref = lattice.ref
    ref.particle.name == "deuteron" || error("reference particle is $(ref.particle.name)")
    ref.Ekin_MeV == 270.0 || error("kinetic energy is $(ref.Ekin_MeV) MeV")
    abs(total_bend_angle(lattice) - 2π) < 1.0e-9 ||
        error("total bend is not 2π: $(total_bend_angle(lattice))")
    for e in lattice.elements
        if e isa Sextupole && e.k2 != 0
            error("sextupole $(e.name) has k2 = $(e.k2); this step reads k2 = 0")
        end
    end
    n_mag = length(magnetic_lattice(lattice).elements)
    n_mag == length(lattice.elements) ||
        error("magnetic_lattice dropped $(length(lattice.elements) - n_mag) elements")
    return nothing
end

function _must_agree(label, a, b, tol)
    d = abs(float(a) - float(b))
    @printf("  %-44s |diff| = %.3e\n", label, d)
    flush(stdout)
    if d > tol
        error("Readings disagree for $label: $a versus $b (|diff| = $d, tolerance $tol). No value was selected.")
    end
    return d
end

function _frac_distance(a, b)
    d = mod(float(a) - float(b), 1.0)
    return min(d, 1.0 - d)
end

function _readings_agree(R, po, obs, pe, C)
    println("Comparing the matrix and the Lie factor")
    flush(stdout)
    dR = maximum(abs, R .- obs.F.R)
    @printf("  %-44s |diff| = %.3e\n", "max |ring_matrix - Lie R|", dR)
    flush(stdout)
    dR <= TOL || error("ring_matrix and the Lie factor R disagree by $dR. No map was selected.")
    _must_agree("eta0 versus eta_slip", pe.eta0, po.eta_slip, TOL)
    _must_agree("eta0 versus Lie periodic eta_slip", pe.eta0, obs.po.eta_slip, TOL)
    _must_agree("alpha0 versus alpha_c", pe.alpha0, po.alpha_c, TOL)
    _must_agree("alpha0 versus Lie periodic alpha_c", pe.alpha0, obs.po.alpha_c, TOL)
    _must_agree("path_expansion eta1 versus eta1_formula", pe.eta1, pe.eta1_formula, TOL)
    _must_agree("Lie path_expansion eta0", pe.eta0, obs.pe.eta0, TOL)
    _must_agree("Lie path_expansion eta1", pe.eta1, obs.pe.eta1, TOL)
    pe.alpha0 > 0 || error("alpha0 = $(pe.alpha0) is not positive; gamma_t is not real")
    gt = 1 / sqrt(pe.alpha0)
    _must_agree("gamma_t versus 1/sqrt(alpha0)", po.gamma_t, gt, TOL)
    _must_agree("path_expansion gamma_t", pe.gamma_t, gt, TOL)
    _must_agree("matrix mu_x versus Lie mu_x", po.x.mu, obs.po.x.mu, TOL)
    _must_agree("matrix mu_y versus Lie mu_y", po.y.mu, obs.po.y.mu, TOL)
    dx = _frac_distance(obs.xi.x[1], po.x.mu / 2π)
    dy = _frac_distance(obs.xi.y[1], po.y.mu / 2π)
    @printf("  %-44s |diff| = %.3e\n", "fractional Qx, Jacobian versus matrix", dx)
    @printf("  %-44s |diff| = %.3e\n", "fractional Qy, Jacobian versus matrix", dy)
    flush(stdout)
    dx <= TOL || error("fractional Qx disagrees: Jacobian $(obs.xi.x[1]), matrix $(po.x.mu / 2π)")
    dy <= TOL || error("fractional Qy disagrees: Jacobian $(obs.xi.y[1]), matrix $(po.y.mu / 2π)")
    sym = symplectic_error(R)
    @printf("  %-44s |diff| = %.3e\n", "symplectic error of ring_matrix", sym)
    flush(stdout)
    sym <= TOL || error("one-turn matrix is not symplectic at the working tolerance: $sym")
    C > 0 || error("path length is not positive")
    return nothing
end

function _transport_agrees(opt, disp, po)
    length(disp.Dx) == length(opt.data.s) ||
        error("dispersion samples $(length(disp.Dx)) versus optics samples $(length(opt.data.s))")
    dD = maximum(abs, disp.Dx .- opt.data.D)
    dDp = maximum(abs, disp.Dpx .- opt.data.Dp)
    _must_agree("sampled Dx versus optics_along", 0.0, dD, 1.0e-9)
    _must_agree("sampled Dpx versus optics_along", 0.0, dDp, 1.0e-9)
    qx = _frac_distance(mod(opt.Qx, 1.0), po.x.mu / 2π)
    qy = _frac_distance(mod(opt.Qy, 1.0), po.y.mu / 2π)
    _must_agree("accumulated fractional Qx", 0.0, qx, TOL)
    _must_agree("accumulated fractional Qy", 0.0, qy, TOL)
    for (k, v) in pairs(opt.closure)
        abs(v) <= TOL || error("optics_along did not close in $k: $v")
    end
    return nothing
end

function _compaction_integral(opt)
    acc = 0.0
    s = opt.data.s
    D = opt.data.D
    for k in 2:length(s)
        s_mid = (s[k] + s[k - 1]) / 2
        idx = findfirst(b -> b[1] <= s_mid < b[2], opt.bounds)
        if idx === nothing
            idx = findfirst(b -> b[1] <= s_mid <= b[2], opt.bounds)
        end
        idx === nothing && error("no element covers s = $s_mid m")
        e = opt.bounds[idx][3]
        acc += curvature(e) * (D[k] + D[k - 1]) / 2 * (s[k] - s[k - 1])
    end
    return acc / opt.path_length
end

function _short_012(e)
    L = elem_length(e)
    return (e isa Drift || e isa Sextupole) && 0.10 < L < 0.13
end

function _arc_cells(lattice)
    els = lattice.elements
    n = length(els)
    cells = NamedTuple[]
    i = 1
    while i <= n - 6
        block = els[i].name == "QF" && els[i + 1].name == "BM" && els[i + 2].name == "QD" &&
                _short_012(els[i + 3]) && els[i + 4].name == "QD" && els[i + 5].name == "BM" &&
                els[i + 6].name == "QF"
        if block && i + 7 <= n && _short_012(els[i + 7])
            push!(cells, (first = i, last = i + 7))
            i += 8
            continue
        elseif block && i + 8 <= n && els[i + 7].name == "OS2" && els[i + 8].name == "OS2"
            push!(cells, (first = i, last = i + 8))
            i += 9
            continue
        end
        i += 1
    end
    return cells
end

function _bends_inside_cells(lattice, cells)
    covered = falses(length(lattice.elements))
    for c in cells, i in c.first:c.last
        covered[i] = true
    end
    for (i, e) in enumerate(lattice.elements)
        e isa SBend && e.name == "BM" && !covered[i] &&
            error("BM at element $i is outside an arc cell")
    end
    return nothing
end

function _index_of_s(svec, s)
    k = findfirst(==(s), svec)
    if k === nothing
        k = findfirst(x -> abs(x - s) <= 1.0e-8, svec)
    end
    k === nothing && error("s = $s m is not on the optics table")
    return k
end

function _forward_dmu(mu, i0, i1)
    i1 >= i0 && return mu[i1] - mu[i0]
    return (mu[end] - mu[i0]) + (mu[i1] - mu[1])
end

function _cell_table(lattice, opt, cells)
    rows = NamedTuple[]
    for (n, c) in enumerate(cells)
        s0 = opt.bounds[c.first][1]
        s1 = opt.bounds[c.last][2]
        i0 = _index_of_s(opt.data.s, s0)
        i1 = _index_of_s(opt.data.s, s1)
        names = join((lattice.elements[i].name for i in c.first:c.last), " ")
        elems = lattice.elements[c.first:c.last]
        cell = Lattice("arc_cell", elems, lattice.ref, lattice.harmonic)
        poc = periodic_optics(ring_matrix(cell), lattice.ref, path_length(cell))
        push!(rows, (
            index = n,
            first = c.first,
            last = c.last,
            names = names,
            s0 = opt.data.s[i0],
            s1 = opt.data.s[i1],
            i0 = i0,
            i1 = i1,
            mu_x = _forward_dmu(opt.data.mux, i0, i1),
            mu_y = _forward_dmu(opt.data.muy, i0, i1),
            periodic_stable = poc.x.stable && poc.y.stable,
            mu_x_periodic = poc.x.stable ? poc.x.mu : NaN,
            mu_y_periodic = poc.y.stable ? poc.y.mu : NaN,
        ))
    end
    return rows
end

function _gap_table(lattice, opt, cells)
    rows = NamedTuple[]
    for n in eachindex(cells)
        j = mod1(n + 1, length(cells))
        i0 = cells[n].i1
        i1 = cells[j].i0
        s0 = opt.data.s[i0]
        s1 = opt.data.s[i1]
        length_m = i1 >= i0 ? s1 - s0 : (opt.path_length - s0) + s1
        kind = length_m > 1.0 ? "straight" : "gap"
        push!(rows, (
            index = n,
            kind = kind,
            s0 = s0,
            s1 = s1,
            length_m = length_m,
            mu_x = _forward_dmu(opt.data.mux, i0, i1),
            mu_y = _forward_dmu(opt.data.muy, i0, i1),
        ))
    end
    return rows
end

function _phases_cover_the_ring(opt, cells, gaps)
    sum_x = sum(c.mu_x for c in cells) + sum(g.mu_x for g in gaps)
    sum_y = sum(c.mu_y for c in cells) + sum(g.mu_y for g in gaps)
    target_x = opt.data.mux[end] - opt.data.mux[1]
    target_y = opt.data.muy[end] - opt.data.muy[1]
    _must_agree("sum of cell and gap phase x", sum_x, target_x, 1.0e-8)
    _must_agree("sum of cell and gap phase y", sum_y, target_y, 1.0e-8)
    return nothing
end

function _superperiod(lattice, opt)
    els = lattice.elements
    n = length(els)
    n % 2 == 0 || return (copies = 1, elements = 0, mu_x = NaN, mu_y = NaN)
    m = n ÷ 2
    for i in 1:m
        if !isequal(els[i], els[i + m])
            return (copies = 1, elements = 0, mu_x = NaN, mu_y = NaN)
        end
    end
    i_half = _index_of_s(opt.data.s, opt.bounds[m][2])
    mu_x = opt.data.mux[i_half] - opt.data.mux[1]
    mu_y = opt.data.muy[i_half] - opt.data.muy[1]
    _must_agree("superperiod phase x versus half tune", mu_x, opt.data.mux[end] / 2, 1.0e-8)
    _must_agree("superperiod phase y versus half tune", mu_y, opt.data.muy[end] / 2, 1.0e-8)
    return (copies = 2, elements = m, mu_x = mu_x, mu_y = mu_y)
end

function _periodic_peaks(s, y; frac = 0.8, sep = 0.5)
    n = length(y)
    thr = frac * maximum(y)
    idx = Int[]
    for k in 2:n - 1
        if y[k] >= y[k - 1] && y[k] >= y[k + 1] && y[k] >= thr
            push!(idx, k)
        end
    end
    if y[1] >= y[2] && y[1] >= y[n - 1] && y[1] >= thr
        pushfirst!(idx, 1)
    end
    kept = Int[]
    for k in idx
        if !isempty(kept) && s[k] - s[kept[end]] < sep
            continue
        end
        push!(kept, k)
    end
    if length(kept) >= 2
        wrap = (s[end] - s[kept[end]]) + s[kept[1]]
        wrap < sep && pop!(kept)
    end
    return kept
end

function _successive_phase(idx, mu)
    n = length(idx)
    dmu = Vector{Float64}(undef, n)
    for i in 1:n
        dmu[i] = i < n ? mu[idx[i + 1]] - mu[idx[i]] :
                 (mu[end] - mu[idx[i]]) + (mu[idx[1]] - mu[1])
    end
    return dmu
end

function _in_api(api, quantity, s)
    for r in api
        if r.quantity == quantity && abs(r.s - s) <= 1.0e-6
            return true
        end
    end
    return false
end

function _curve_peaks(opt, disp, api)
    d = opt.data
    specs = (
        (quantity = "D_x", y = d.D, mu = d.mux),
        (quantity = "beta_x", y = d.betx, mu = d.mux),
        (quantity = "beta_y", y = d.bety, mu = d.muy),
    )
    rows = NamedTuple[]
    for spec in specs
        idx = _periodic_peaks(d.s, spec.y)
        isempty(idx) && continue
        dmu = _successive_phase(idx, spec.mu)
        for (j, k) in enumerate(idx)
            push!(rows, (
                quantity = spec.quantity,
                s = d.s[k],
                value = spec.y[k],
                betx = d.betx[k],
                bety = d.bety[k],
                D = d.D[k],
                Dy = disp.Dy[k],
                mux = d.mux[k],
                muy = d.muy[k],
                dmu_next = dmu[j],
                in_find_optics_maxima = _in_api(api, spec.quantity, d.s[k]),
            ))
        end
    end
    return rows
end

function _smooth_phi_s(eta0)
    # The smooth-model tune is stable when eta0 * cos(phi_s) > 0.
    # That phase origin is opposite the notes. Used only to cross-check Q_s.
    return eta0 > 0 ? 0.0 : float(π)
end

function _smooth_model_tunes(ref, eta0, regimes)
    phi = _smooth_phi_s(eta0)
    out = NamedTuple[]
    for reg in (regimes.linear, regimes.cubic)
        qs = synchrotron_tune(ref, eta0, reg.voltage_MV, reg.harmonic, phi; T = 1.0)
        isfinite(qs) || error("synchrotron_tune returned a non-finite Q_s for $(reg.name)")
        _must_agree("Q_s note versus smooth model ($(reg.name))", reg.Qs, qs, 1.0e-9)
        push!(out, (name = reg.name, Qs = qs, phi_s_rad = phi))
    end
    return out
end

function _fmt(x)
    return @sprintf("%.8e", float(x))
end

function _csv_escape(s)
    t = string(s)
    if occursin(',', t) || occursin('"', t)
        return "\"" * replace(t, '"' => "\"\"") * "\""
    end
    return t
end

function _write_optics_csv(path, opt, disp)
    d = opt.data
    open(path, "w") do io
        println(io, "s_m,betx_m,alfx,bety_m,alfy,Dx_m,Dpx,Dy_m,Dpy,mux_rad,muy_rad")
        for k in eachindex(d.s)
            @printf(io, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                    d.s[k], d.betx[k], d.alfx[k], d.bety[k], d.alfy[k],
                    d.D[k], d.Dp[k], disp.Dy[k], disp.Dpy[k], d.mux[k], d.muy[k])
        end
    end
    return path
end

function _write_numbers_csv(path, opt, po, pe, C, cells, half, alpha_geom)
    mux = [c.mu_x for c in cells]
    muy = [c.mu_y for c in cells]
    rows = [
        ("C_m", C, "path_length"),
        ("Qx", opt.Qx, "optics_along"),
        ("Qy", opt.Qy, "optics_along"),
        ("Qx_frac", po.x.mu / 2π, "periodic_optics"),
        ("Qy_frac", po.y.mu / 2π, "periodic_optics"),
        ("alpha0", pe.alpha0, "path_expansion"),
        ("alpha1", pe.alpha1, "path_expansion"),
        ("eta0", pe.eta0, "path_expansion"),
        ("eta1", pe.eta1, "path_expansion"),
        ("gamma_t", pe.gamma_t, "path_expansion"),
        ("alpha0_dispersion_integral", alpha_geom, "(1/C) int D/rho ds"),
        ("arc_cell_mu_x_mean_rad", sum(mux) / length(mux), "optics_along across arc cells"),
        ("arc_cell_mu_y_mean_rad", sum(muy) / length(muy), "optics_along across arc cells"),
        ("arc_cell_mu_x_spread_rad", maximum(mux) - minimum(mux), "max minus min"),
        ("arc_cell_mu_y_spread_rad", maximum(muy) - minimum(muy), "max minus min"),
        ("superperiod_mu_x_rad", half.mu_x, "optics_along over one copy"),
        ("superperiod_mu_y_rad", half.mu_y, "optics_along over one copy"),
    ]
    open(path, "w") do io
        println(io, "quantity,value,source")
        for (q, v, src) in rows
            println(io, q, ",", _fmt(v), ",", _csv_escape(src))
        end
    end
    return path
end

function _write_cells_csv(path, cells, gaps)
    open(path, "w") do io
        println(io, "kind,index,s_start_m,s_end_m,length_m,mu_x_rad,mu_y_rad,mu_x_deg,mu_y_deg,mu_x_periodic_rad,mu_y_periodic_rad,names")
        for c in cells
            @printf(io, "arc_cell,%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%s\n",
                    c.index, c.s0, c.s1, c.s1 - c.s0, c.mu_x, c.mu_y,
                    rad2deg(c.mu_x), rad2deg(c.mu_y), c.mu_x_periodic, c.mu_y_periodic,
                    _csv_escape(c.names))
        end
        for g in gaps
            g.length_m > 1.0e-9 || continue
            @printf(io, "%s,%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,,,\n",
                    g.kind, g.index, g.s0, g.s1, g.length_m, g.mu_x, g.mu_y,
                    rad2deg(g.mu_x), rad2deg(g.mu_y))
        end
    end
    return path
end

function _write_maxima_csv(path, peaks)
    open(path, "w") do io
        println(io, "quantity,s_m,value,betx_m,bety_m,D_m,Dy_m,mux_rad,muy_rad,dmu_to_next_rad,dmu_to_next_deg,in_find_optics_maxima")
        for r in peaks
            println(io, r.quantity, ",", _fmt(r.s), ",", _fmt(r.value), ",",
                    _fmt(r.betx), ",", _fmt(r.bety), ",", _fmt(r.D), ",", _fmt(r.Dy), ",",
                    _fmt(r.mux), ",", _fmt(r.muy), ",", _fmt(r.dmu_next), ",",
                    _fmt(rad2deg(r.dmu_next)), ",", r.in_find_optics_maxima)
        end
    end
    return path
end

function _xml(s)
    return replace(replace(replace(string(s), "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")
end

function _limits(v)
    lo = minimum(v)
    hi = maximum(v)
    span = hi - lo
    if span < 1.0e-12
        span = max(1.0e-6, abs(lo) * 0.2)
        mid = (lo + hi) / 2
        return mid - span / 2, mid + span / 2
    end
    pad = 0.08 * span
    return lo - pad, hi + pad
end

function _tick_label(x, span)
    if span >= 20
        return @sprintf("%.0f", x)
    elseif span >= 2
        return @sprintf("%.1f", x)
    elseif span >= 0.2
        return @sprintf("%.2f", x)
    elseif span >= 1.0e-3
        return @sprintf("%.3f", x)
    else
        return @sprintf("%.1e", x)
    end
end

function _polyline(xs, ys, xp, yp)
    parts = String[]
    for i in eachindex(xs)
        push!(parts, @sprintf("%.2f,%.2f", xp(xs[i]), yp(ys[i])))
    end
    return join(parts, " ")
end

function _write_svg(path, opt, disp, peaks)
    d = opt.data
    W, H = 980.0, 860.0
    left, right, top = 78.0, 24.0, 42.0
    pw = W - left - right
    ph, gap = 150.0, 16.0
    panels = (
        (title = "beta", ylabel = "beta [m]",
         curves = (("beta_x", d.betx, "#1d4e89"), ("beta_y", d.bety, "#c0392b")),
         marks = ("beta_x", "beta_y")),
        (title = "alpha", ylabel = "alpha",
         curves = (("alpha_x", d.alfx, "#1d4e89"), ("alpha_y", d.alfy, "#c0392b")),
         marks = ()),
        (title = "D", ylabel = "D [m]",
         curves = (("D_x", d.D, "#1e7f4f"), ("D_y", disp.Dy, "#7f8c8d")),
         marks = ("D_x",)),
        (title = "Dp", ylabel = "D'",
         curves = (("D_x'", d.Dp, "#1e7f4f"), ("D_y'", disp.Dpy, "#7f8c8d")),
         marks = ()),
    )
    open(path, "w") do io
        println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
        @printf(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" viewBox=\"0 0 %.0f %.0f\">\n",
                W, H, W, H)
        println(io, "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>")
        println(io, "<text x=\"78\" y=\"26\" font-family=\"sans-serif\" font-size=\"16\" fill=\"#222\">closed magnetic_2 linear optics</text>")
        for (p, panel) in enumerate(panels)
            y0 = top + (p - 1) * (ph + gap)
            ys = Float64[]
            for c in panel.curves
                append!(ys, c[2])
            end
            ylo, yhi = _limits(ys)
            span = yhi - ylo
            xp = x -> left + (x - 0.0) / opt.path_length * pw
            yp = y -> y0 + (yhi - y) / (yhi - ylo) * ph
            @printf(io, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"#fcfcfc\" stroke=\"#dddddd\"/>\n",
                    left, y0, pw, ph)
            println(io, "<text x=\"8\" y=\"", y0 + ph / 2, "\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#333\" transform=\"rotate(-90 14,",
                    y0 + ph / 2, ")\">", _xml(panel.ylabel), "</text>")
            for tick in 0:4
                xv = tick * opt.path_length / 4
                px = xp(xv)
                @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#eeeeee\"/>\n",
                        px, y0, px, y0 + ph)
                if p == length(panels)
                    @printf(io, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#444\">%s</text>\n",
                            px, y0 + ph + 14, _tick_label(xv, opt.path_length))
                end
            end
            for tick in 0:3
                yv = ylo + tick * span / 3
                py = yp(yv)
                @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#f2f2f2\"/>\n",
                        left, py, left + pw, py)
                @printf(io, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#444\">%s</text>\n",
                        left - 6, py + 4, _xml(_tick_label(yv, span)))
            end
            for (name, y, color) in panel.curves
                println(io, "<polyline fill=\"none\" stroke=\"", color, "\" stroke-width=\"1.4\" points=\"",
                        _polyline(d.s, y, xp, yp), "\"/>")
            end
            for (j, (name, y, color)) in enumerate(panel.curves)
                lx = left + pw - 150
                ly = y0 + 16 + (j - 1) * 14
                @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-width=\"2\"/>\n",
                        lx, ly - 3, lx + 16, ly - 3, color)
                @printf(io, "<text x=\"%.2f\" y=\"%.2f\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#222\">%s</text>\n",
                        lx + 20, ly, _xml(name))
            end
            for quantity in panel.marks
                for r in peaks
                    r.quantity == quantity || continue
                    mx = xp(r.s)
                    my = yp(r.value)
                    if quantity == "beta_y"
                        @printf(io, "<rect x=\"%.2f\" y=\"%.2f\" width=\"6\" height=\"6\" fill=\"#c0392b\" stroke=\"#222\" stroke-width=\"0.4\"/>\n",
                                mx - 3, my - 3)
                    elseif quantity == "D_x"
                        println(io, @sprintf("<polygon points=\"%.2f,%.2f %.2f,%.2f %.2f,%.2f\" fill=\"#1e7f4f\" stroke=\"#222\" stroke-width=\"0.4\"/>",
                                             mx, my - 4, mx - 4, my + 3, mx + 4, my + 3))
                    else
                        @printf(io, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"3.2\" fill=\"#1d4e89\" stroke=\"#222\" stroke-width=\"0.4\"/>\n",
                                mx, my)
                    end
                end
            end
        end
        y_strip = top + length(panels) * (ph + gap) + 8
        @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#222\" stroke-width=\"1\"/>\n",
                left, y_strip + 18, left + pw, y_strip + 18)
        xp = x -> left + x / opt.path_length * pw
        for (s0, s1, e) in opt.bounds
            x0 = xp(s0)
            w = max(xp(s1) - x0, 0.4)
            if e isa SBend
                @printf(io, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"16\" fill=\"#b0b0b0\"/>\n", x0, y_strip + 10, w)
            elseif e isa Quadrupole && e.k1 > 0
                @printf(io, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"14\" fill=\"#d64545\"/>\n", x0, y_strip + 4, w)
            elseif e isa Quadrupole
                @printf(io, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"14\" fill=\"#3a6ea5\"/>\n", x0, y_strip + 18, w)
            elseif e isa Sextupole
                @printf(io, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"10\" fill=\"#3d8c5a\"/>\n", x0, y_strip + 13, w)
            end
        end
        println(io, "<text x=\"78\" y=\"", y_strip + 48, "\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#444\">s [m]. Markers: beta_x circles, beta_y squares, D_x triangles. Strip: bends grey, focusing quadrupoles red, defocusing quadrupoles blue, sextupoles green.</text>")
        println(io, "</svg>")
    end
    return path
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
    isfinite(x) || error("refusing a non-finite number in the regime file: $x")
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

function _write_rf_json(path, ref, C, opt, po, pe, cells, half, regimes, qs_smooth, alpha_geom, xi_x, xi_y)
    mux = [c.mu_x for c in cells]
    muy = [c.mu_y for c in cells]
    doc = (
        lattice = "lattices/magnetic_2_closed.bell",
        step = "linear_optics",
        reference = (
            particle = ref.particle.name,
            charge_number = ref.particle.charge,
            kinetic_MeV = ref.Ekin_MeV,
            mass_MeV = ref.mass_MeV,
            gamma = ref.gamma,
            beta = ref.beta,
            energy_MeV = ref.gamma * ref.mass_MeV,
        ),
        optics = (
            C_m = C,
            Qx = opt.Qx,
            Qy = opt.Qy,
            Qx_frac = po.x.mu / 2π,
            Qy_frac = po.y.mu / 2π,
            alpha0 = pe.alpha0,
            alpha1 = pe.alpha1,
            eta0 = pe.eta0,
            eta1 = pe.eta1,
            gamma_t = pe.gamma_t,
            below_transition = pe.eta0 < 0,
            arc_cell_mu_x_mean_rad = sum(mux) / length(mux),
            arc_cell_mu_y_mean_rad = sum(muy) / length(muy),
            arc_cell_mu_x_spread_rad = maximum(mux) - minimum(mux),
            arc_cell_mu_y_spread_rad = maximum(muy) - minimum(muy),
            arc_cell_mu_x_periodic_rad = isfinite(cells[1].mu_x_periodic) ? cells[1].mu_x_periodic : nothing,
            arc_cell_mu_y_periodic_rad = isfinite(cells[1].mu_y_periodic) ? cells[1].mu_y_periodic : nothing,
            superperiod_copies = half.copies,
            superperiod_mu_x_rad = isfinite(half.mu_x) ? half.mu_x : nothing,
            superperiod_mu_y_rad = isfinite(half.mu_y) ? half.mu_y : nothing,
            alpha0_dispersion_integral = alpha_geom,
        ),
        xi_for_lambda = (
            xi_x = xi_x,
            xi_y = xi_y,
            role = "Natural first chromaticity of the k2=0 Lie polynomial, used only to size Lambda. Not a chromaticity result of this step.",
        ),
        selection_cuts = regimes.cuts,
        placeholder_action = regimes.placeholder,
        stationary_bucket = (
            sin_phi_s = 0.0,
            phi_s_rad = regimes.phi_s_rad,
            cos_phi_s = cos(regimes.phi_s_rad),
            T_s = regimes.T_s,
            harmonic = regimes.harmonic,
            delta_m = regimes.delta_m,
        ),
        regimes = (
            linear = regimes.linear,
            cubic_phase = regimes.cubic,
        ),
        Qs_smooth_model = qs_smooth,
    )
    open(path, "w") do io
        _json(io, doc, 0)
        println(io)
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

function _ineq_rows(reg)
    return [
        ("abs(eta1*delta_m/eta0)", _fmt(reg.abs_eta1_delta_m_over_eta0), "much less than 1"),
        ("abs(Lambda/(eta0*delta_m))", _fmt(reg.abs_Lambda_over_eta0_delta_m), "much less than 1"),
        ("phi_hat", _fmt(reg.phi_hat), reg.name == "linear" ? "much less than 1" : "not the linear condition"),
        ("phi_hat^2", _fmt(reg.phi_hat_squared), reg.name == "linear" ? "follows from phi_hat" : "much less than 1"),
        ("delta_m", _fmt(reg.delta_m), "much less than 1"),
        ("omega_s T_s", _fmt(reg.omega_s_Ts), "much less than 1"),
        ("2 pi {nu_x}", _fmt(reg.two_pi_frac_nu_x), "printed"),
        ("2 pi {nu_y}", _fmt(reg.two_pi_frac_nu_y), "printed"),
        ("omega_s T_s / (2 pi {nu_x})", _fmt(reg.omega_s_Ts_over_two_pi_frac_nu_x), "much less than 1"),
        ("omega_s T_s / (2 pi {nu_y})", _fmt(reg.omega_s_Ts_over_two_pi_frac_nu_y), "much less than 1"),
        ("2 pi dist(nu_x, integer)", _fmt(reg.two_pi_nearest_nu_x), "phase to the nearest integer"),
        ("2 pi dist(nu_y, integer)", _fmt(reg.two_pi_nearest_nu_y), "phase to the nearest integer"),
        ("omega_s T_s / that x distance", _fmt(reg.omega_s_Ts_over_two_pi_nearest_nu_x), "much less than 1"),
        ("omega_s T_s / that y distance", _fmt(reg.omega_s_Ts_over_two_pi_nearest_nu_y), "much less than 1"),
        ("Q_s", _fmt(reg.Qs), "omega_s T_s / 2 pi"),
        ("V_MV", _fmt(reg.voltage_MV), "stationary bucket"),
        ("formula delta_eq", _fmt(reg.formula_delta_eq), "notes, with placeholder Lambda"),
    ]
end

function _write_rf_md(path, ref, pe, regimes, qs_smooth)
    ph = regimes.placeholder
    lines = String[]
    push!(lines, "# Stationary RF regimes for the closed magnetic ring")
    push!(lines, "")
    push!(lines, "Chosen after `eta0`, `eta1`, and the tunes of this step. The lattice harmonic is kept. `V = 0.2` MV and `h = 5` are not reused.")
    push!(lines, "")
    push!(lines, "The bucket is stationary: `sin phi_s = 0` and `cos phi_s = -sign(eta0)`, so a positive voltage is stable. `phi_s = $(regimes.phi_s_rad)` rad.")
    push!(lines, "")
    push!(lines, @sprintf("Reference energy %.6f MeV, beta %.8f, gamma %.8f, eta0 %.8e, eta1 %.8e.",
                          ref.gamma * ref.mass_MeV, ref.beta, ref.gamma, pe.eta0, pe.eta1))
    push!(lines, @sprintf("Revolution period T_s = %.8e s.", regimes.T_s))
    push!(lines, "")
    push!(lines, "Selection cuts, used as a concrete stand-in for much-less-than-one:")
    push!(lines, "")
    push!(lines, _md_table(["cut", "value"], [
        ["abs(eta1*delta_m/eta0)", _fmt(regimes.cuts.eta_ratio_max)],
        ["abs(Lambda/(eta0*delta_m))", _fmt(regimes.cuts.lambda_ratio_max)],
        ["delta_m", _fmt(regimes.cuts.delta_m_max)],
        ["omega_s T_s", _fmt(regimes.cuts.omega_s_Ts_max)],
        ["omega_s T_s / betatron phase", _fmt(regimes.cuts.omega_over_betatron_max)],
        ["linear phi_hat target", _fmt(regimes.cuts.linear_phi_hat)],
        ["cubic phi_hat target", _fmt(regimes.cuts.cubic_phi_hat)],
    ]))
    push!(lines, "")
    push!(lines, @sprintf("Shared harmonic h = %d, shared delta_m = %.8e.", regimes.harmonic, regimes.delta_m))
    push!(lines, "")
    push!(lines, "Placeholder action, not a measured amplitude. Peak displacement `a` at the beta maximum, `J = a^2 / (2 beta_max)`. Lambda uses the natural first chromaticity of the same Lie polynomial (`xi_x = $( _fmt(ph.xi_x) )`, `xi_y = $( _fmt(ph.xi_y) )`). That chromaticity is only a size for Lambda. The chromaticity table is a later step.")
    push!(lines, "")
    push!(lines, _md_table(["quantity", "value"], [
        ["amplitude_m", _fmt(ph.amplitude_m)],
        ["beta_x_max_m", _fmt(ph.beta_x_m)],
        ["beta_y_max_m", _fmt(ph.beta_y_m)],
        ["J_x_m", _fmt(ph.Jx)],
        ["J_y_m", _fmt(ph.Jy)],
        ["Lambda", _fmt(ph.Lambda)],
    ]))
    push!(lines, "")
    for (reg, qs) in zip((regimes.linear, regimes.cubic), qs_smooth)
        title = reg.name == "linear" ? "A. Linear regime" : "B. Cubic-phase regime"
        push!(lines, "## $title")
        push!(lines, "")
        if reg.name == "linear"
            push!(lines, "Target phase half-range 0.05, so `sin phi ≈ phi` applies.")
        else
            push!(lines, "Target phase half-range 0.30. `phi_hat` is not small. `phi_hat^2 = 0.09` is the condition that replaces `phi_hat ≪ 1`. The other conditions stay.")
        end
        push!(lines, "")
        push!(lines, _md_table(["quantity", "value", "role"], [collect(r) for r in _ineq_rows(reg)]))
        push!(lines, "")
        push!(lines, @sprintf("Smooth-model Q_s at its own stable phase (%.6f rad) is %.8e. The absolute difference from the note tune was required to be at most 1e-9.",
                              qs.phi_s_rad, qs.Qs))
        push!(lines, "")
    end
    push!(lines, "`formula delta_eq` evaluates `-(eta1 * delta_m^2 / 2 * f + Lambda) / eta0`, with `f = 1` in the linear regime and `f = 1 - phi_hat^2/32` in the cubic regime. It is not a tracking result.")
    push!(lines, "")
    write(path, join(lines, "\n"))
    return path
end

function _write_report(path, lattice, latpath, opt, po, pe, obs, cells, gaps, half, peaks, api,
                       regimes, qs_smooth, alpha_geom, d_alpha, disp)
    ref = lattice.ref
    mux = [c.mu_x for c in cells]
    muy = [c.mu_y for c in cells]
    lines = String[]
    push!(lines, "# Linear optics of closed magnetic_2")
    push!(lines, "")
    push!(lines, "Lattice file: `$latpath`. Every sextupole has `k2 = 0`. Reference: $(ref.particle.name), kinetic energy $(ref.Ekin_MeV) MeV, charge number $(ref.particle.charge), harmonic $(lattice.harmonic).")
    push!(lines, "")
    push!(lines, @sprintf("gamma = %.8f, beta = %.8f, total energy = %.6f MeV, brho = %.6f T m.",
                          ref.gamma, ref.beta, ref.gamma * ref.mass_MeV, ref.brho_Tm))
    push!(lines, "")
    push!(lines, "The one-turn matrix is `ring_matrix` (`linear_matrix` of each element, `R = exp(L S Hess H)`). The polynomial is `path_expansion` of `observables`, which builds `ring_lie_factor`. `periodic_optics` is taken on the matrix and on the Lie factor `R`. `optics_along` transports the matrix Twiss. Where two readings name the same quantity, both are printed and the run stops if they differ by more than $TOL. Nothing was averaged and nothing was tuned.")
    push!(lines, "")
    push!(lines, @sprintf("Largest |matrix − Lie R| and the slip, compaction, and tune residuals are below %.1e. Symplectic error of the matrix is %.3e. `eta1` matches `eta1_formula` from the same expansion. `(1/C) ∫ D/ρ ds` differs from `alpha0` by %.3e.",
                          TOL, symplectic_error(ring_matrix(lattice)), d_alpha))
    push!(lines, "")
    push!(lines, "## One-turn numbers")
    push!(lines, "")
    push!(lines, "`Qx` and `Qy` include the integer from `optics_along`. Fractional tunes are `mu/2π` from `periodic_optics` of the matrix. `alpha0`, `eta0`, `eta1`, and `gamma_t` are from `path_expansion`.")
    push!(lines, "")
    push!(lines, _md_table(["quantity", "value"], [
        ["C [m]", _fmt(opt.path_length)],
        ["Qx", _fmt(opt.Qx)],
        ["Qy", _fmt(opt.Qy)],
        ["Qx fractional", _fmt(po.x.mu / 2π)],
        ["Qy fractional", _fmt(po.y.mu / 2π)],
        ["alpha0", _fmt(pe.alpha0)],
        ["alpha1", _fmt(pe.alpha1)],
        ["eta0", _fmt(pe.eta0)],
        ["eta1", _fmt(pe.eta1)],
        ["gamma_t", _fmt(pe.gamma_t)],
        ["1/gamma^2", _fmt(1 / ref.gamma^2)],
        ["(1/C) int D/rho ds", _fmt(alpha_geom)],
        ["beta_x(0) [m]", _fmt(po.x.beta)],
        ["alpha_x(0)", _fmt(po.x.alpha)],
        ["beta_y(0) [m]", _fmt(po.y.beta)],
        ["alpha_y(0)", _fmt(po.y.alpha)],
        ["D_x(0) [m]", _fmt(po.D4[1])],
        ["D_x'(0)", _fmt(po.D4[2])],
        ["D_y(0) [m]", _fmt(po.D4[3])],
        ["D_y'(0)", _fmt(po.D4[4])],
        ["max abs(D_y) [m]", _fmt(maximum(abs, disp.Dy))],
        ["max abs(D_y')", _fmt(maximum(abs, disp.Dpy))],
        ["max beta_x [m]", _fmt(maximum(opt.data.betx))],
        ["max beta_y [m]", _fmt(maximum(opt.data.bety))],
        ["max D_x [m]", _fmt(maximum(opt.data.D))],
        ["min D_x [m]", _fmt(minimum(opt.data.D))],
    ]))
    push!(lines, "")
    push!(lines, "Closure of `optics_along` (relative beta, absolute alpha, absolute dispersion, fractional tune):")
    push!(lines, "")
    cl = opt.closure
    push!(lines, _md_table(["quantity", "residual"], [
        ["betx", _fmt(cl.betx)], ["alfx", _fmt(cl.alfx)],
        ["bety", _fmt(cl.bety)], ["alfy", _fmt(cl.alfy)],
        ["D", _fmt(cl.D)], ["Dp", _fmt(cl.Dp)],
        ["Qx_frac", _fmt(cl.Qx_frac)], ["Qy_frac", _fmt(cl.Qy_frac)],
    ]))
    push!(lines, "")
    push!(lines, "## Cell phase advance")
    push!(lines, "")
    push!(lines, "An arc cell is the repeated block `QF BM QD (L = 0.12) QD BM QF` closed by one `L = 0.12` element, or by two `OS2` drifts where that element is split. Sextupoles are at `k2 = 0`. The in-ring phase is the advance of `optics_along` from the cell entrance to the cell exit. It is not set to 90 degrees. The periodic column is `mu` of `periodic_optics` on that cell's own matrix. Those two phases answer different questions: the periodic column assumes the cell is repeated by itself, and the in-ring column is the advance inside this ring, straights included in the periodic Twiss.")
    push!(lines, "")
    push!(lines, @sprintf("In-ring mean: mu_x = %.8f rad (%.4f deg), mu_y = %.8f rad (%.4f deg). Spread max−min: %.3e rad in x, %.3e rad in y.",
                          sum(mux) / length(mux), rad2deg(sum(mux) / length(mux)),
                          sum(muy) / length(muy), rad2deg(sum(muy) / length(muy)),
                          maximum(mux) - minimum(mux), maximum(muy) - minimum(muy)))
    push!(lines, "")
    cell_md = Vector{Vector{String}}()
    for c in cells
        push!(cell_md, [string(c.index), @sprintf("%.4f", c.s0), @sprintf("%.4f", c.s1),
                        @sprintf("%.4f", rad2deg(c.mu_x)), @sprintf("%.4f", rad2deg(c.mu_y)),
                        c.periodic_stable ? @sprintf("%.4f", rad2deg(c.mu_x_periodic)) : "unstable",
                        c.periodic_stable ? @sprintf("%.4f", rad2deg(c.mu_y_periodic)) : "unstable",
                        c.names])
    end
    push!(lines, _md_table(["cell", "s_start [m]", "s_end [m]", "mu_x [deg]", "mu_y [deg]",
                            "periodic mu_x [deg]", "periodic mu_y [deg]", "elements"], cell_md))
    push!(lines, "")
    if half.copies == 2
        push!(lines, @sprintf("The element list is two copies of %d elements. Phase over the first copy: mu_x = %.8f rad (%.4f deg), mu_y = %.8f rad (%.4f deg), equal to half the one-turn phase within 1e-8 rad.",
                              half.elements, half.mu_x, rad2deg(half.mu_x), half.mu_y, rad2deg(half.mu_y)))
    else
        push!(lines, "The element list does not split into two identical copies.")
    end
    push!(lines, "")
    gap_md = Vector{Vector{String}}()
    for g in gaps
        g.kind == "straight" || continue
        push!(gap_md, [string(g.index), @sprintf("%.4f", g.length_m),
                       @sprintf("%.4f", rad2deg(g.mu_x)), @sprintf("%.4f", rad2deg(g.mu_y))])
    end
    push!(lines, "Straights are the gaps longer than 1 m between arc cells. Their phase is included so that cells plus gaps sum to the one-turn phase.")
    push!(lines, "")
    push!(lines, _md_table(["after cell", "length [m]", "mu_x [deg]", "mu_y [deg]"], gap_md))
    push!(lines, "")
    push!(lines, "## Maxima")
    push!(lines, "")
    push!(lines, "Local maxima of `D_x`, `beta_x`, and `beta_y` on the `optics_along` table, at least 0.8 of that curve's global maximum, and at least 0.5 m apart. A maximum on the periodic boundary `s = 0` is kept when it is above both neighbours. `dmu_to_next` is the phase to the next maximum of the same quantity, using `mux` for `D_x` and `beta_x` and `muy` for `beta_y`. The step that crosses `s = C` is included, so the advances of one quantity sum to the one-turn phase.")
    push!(lines, "")
    push!(lines, "`find_optics_maxima` applies the same separation and fraction, and also requires `D_x > 0.5` m. The last column says whether that call returned the same peak. A step can exceed 360 degrees because the advance is betatron phase and the tune is above 2.")
    push!(lines, "")
    n_api = length(api)
    n_flag = count(r -> r.in_find_optics_maxima, peaks)
    push!(lines, "Curve maxima: $(length(peaks)). `find_optics_maxima` rows: $n_api. Curve maxima also returned by that call: $n_flag.")
    push!(lines, "")
    max_md = Vector{Vector{String}}()
    for r in peaks
        push!(max_md, [r.quantity, @sprintf("%.4f", r.s), @sprintf("%.6f", r.value),
                       @sprintf("%.4f", r.betx), @sprintf("%.4f", r.bety), @sprintf("%.4f", r.D),
                       @sprintf("%.4f", rad2deg(r.dmu_next)), r.in_find_optics_maxima ? "yes" : "no"])
    end
    push!(lines, _md_table(["quantity", "s [m]", "value", "beta_x [m]", "beta_y [m]", "D_x [m]",
                            "phase to next [deg]", "find_optics_maxima"], max_md))
    push!(lines, "")
    push!(lines, "Figure: `02_optics.svg`. Samples: `02_optics.csv` (ds = $DS m).")
    push!(lines, "")
    push!(lines, "## RF regimes")
    push!(lines, "")
    push!(lines, "Two stationary buckets share `h = $(regimes.harmonic)` and `delta_m = $(_fmt(regimes.delta_m))`. The linear voltage realises `phi_hat = 0.05`. The cubic-phase voltage realises `phi_hat = 0.30`. Full tables, the placeholder action, and the inequality values are in `02_rf_regimes.md` and `02_rf_regimes.json`.")
    push!(lines, "")
    push!(lines, _md_table(["regime", "h", "V [MV]", "delta_m", "phi_hat", "abs(eta1*dm/eta0)",
                            "abs(Lambda/(eta0*dm))", "omega_s T_s", "Q_s"],
        [
            ["linear", string(regimes.linear.harmonic), _fmt(regimes.linear.voltage_MV),
             _fmt(regimes.linear.delta_m), _fmt(regimes.linear.phi_hat),
             _fmt(regimes.linear.abs_eta1_delta_m_over_eta0),
             _fmt(regimes.linear.abs_Lambda_over_eta0_delta_m),
             _fmt(regimes.linear.omega_s_Ts), _fmt(regimes.linear.Qs)],
            ["cubic_phase", string(regimes.cubic.harmonic), _fmt(regimes.cubic.voltage_MV),
             _fmt(regimes.cubic.delta_m), _fmt(regimes.cubic.phi_hat),
             _fmt(regimes.cubic.abs_eta1_delta_m_over_eta0),
             _fmt(regimes.cubic.abs_Lambda_over_eta0_delta_m),
             _fmt(regimes.cubic.omega_s_Ts), _fmt(regimes.cubic.Qs)],
        ]))
    push!(lines, "")
    push!(lines, @sprintf("Natural chromaticity used only inside Lambda: xi_x = %.6f, xi_y = %.6f (first derivative from `observables`, k2 = 0).",
                          obs.xi.x[2], obs.xi.y[2]))
    push!(lines, "")
    write(path, join(lines, "\n"))
    return path
end

function _print_summary(opt, pe, cells, regimes, paths)
    mux = sum(c.mu_x for c in cells) / length(cells)
    muy = sum(c.mu_y for c in cells) / length(cells)
    println()
    @printf("Qx = %.10f    Qy = %.10f\n", opt.Qx, opt.Qy)
    @printf("alpha0 = %.8e    eta0 = %.8e    eta1 = %.8e    gamma_t = %.8f    C = %.8f m\n",
            pe.alpha0, pe.eta0, pe.eta1, pe.gamma_t, opt.path_length)
    @printf("arc-cell phase  mu_x = %.6f rad (%.3f deg)    mu_y = %.6f rad (%.3f deg)\n",
            mux, rad2deg(mux), muy, rad2deg(muy))
    for reg in (regimes.linear, regimes.cubic)
        @printf("%s  h = %d  V = %.8e MV  delta_m = %.8e  phi_hat = %.6f  |eta1 dm/eta0| = %.6e  |Lambda/(eta0 dm)| = %.6e  omega_s Ts = %.6e\n",
                reg.name, reg.harmonic, reg.voltage_MV, reg.delta_m, reg.phi_hat,
                reg.abs_eta1_delta_m_over_eta0, reg.abs_Lambda_over_eta0_delta_m, reg.omega_s_Ts)
    end
    println("outputs:")
    for p in paths
        println("  ", p)
    end
    flush(stdout)
    return nothing
end

main()
