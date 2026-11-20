# Step 1. Write the magnetic_2 element table and close the sector bends.
#
#   julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/01_lattice.jl

import Pkg

const STUDY_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(STUDY_DIR)

using Bell
using MagneticRing
using Printf

const OUT_DIR = joinpath(STUDY_DIR, "output")

function main()
    println("Loading ", BELL_LATTICE)
    flush(stdout)
    source = load_magnetic_2()
    layout = source_layout(source)
    token = source_bm_angle_token()
    println("BM file token = ", token)
    println("BM loaded angle = ", layout.bm_angle)
    println("Reading open ring (one-turn matrix and optics along s)")
    flush(stdout)
    open_r = ring_readings(source)

    closed = close_bm_angles(source)
    source_layout(closed)
    only_bm_angle_changed(source, closed) ||
        error("closing the bends changed something other than the BM angle")
    all(e -> !(e isa SBend && e.name == "BM") || e.angle == TARGET_BM_ANGLE, closed.elements) ||
        error("a BM angle is not pi/8")
    path_length_unchanged = open_r.path_length == path_length(closed)
    path_length_unchanged || error("path length changed when only the bend angle changed")

    println("Reading closed ring")
    flush(stdout)
    closed_r = ring_readings(closed)

    written = write_closed_lattice()
    reloaded = load_closed_lattice()
    same_elements(closed, reloaded) ||
        error("reloaded closed lattice differs from the in-memory ring: $written")
    println("Wrote ", written)
    flush(stdout)

    roundtrip_note = _write_structure_roundtrip(closed)

    mkpath(OUT_DIR)
    elements_csv = joinpath(OUT_DIR, "unique_elements.csv")
    sequence_csv = joinpath(OUT_DIR, "sequence.csv")
    closure_csv = joinpath(OUT_DIR, "closure.csv")
    report_md = joinpath(OUT_DIR, "01_lattice.md")
    figure = joinpath(OUT_DIR, "bend_deficit.svg")

    _write_unique_csv(elements_csv, source, closed)
    _write_sequence_csv(sequence_csv, source, closed)
    rows = _closure_rows(layout, open_r, closed_r)
    _write_closure_csv(closure_csv, rows)
    _write_report(report_md, layout, token, source, closed, open_r, closed_r, rows,
                  written, roundtrip_note)
    _write_deficit_svg(figure, open_r.bend_minus_2pi, closed_r.bend_minus_2pi)

    _check_tune_closure(open_r, "open")
    _check_tune_closure(closed_r, "closed")

    println()
    println("BM angle: ", layout.bm_angle, " -> ", TARGET_BM_ANGLE)
    println("total bend: ", open_r.total_bend_angle, " -> ", closed_r.total_bend_angle)
    println("bend - 2pi: ", open_r.bend_minus_2pi, " -> ", closed_r.bend_minus_2pi)
    println("path length: ", open_r.path_length, " -> ", closed_r.path_length,
            "  delta = ", closed_r.path_length - open_r.path_length)
    println("Qx: ", open_r.Qx, " -> ", closed_r.Qx)
    println("Qy: ", open_r.Qy, " -> ", closed_r.Qy)
    println("Qx_frac: ", open_r.Qx_frac, " -> ", closed_r.Qx_frac)
    println("Qy_frac: ", open_r.Qy_frac, " -> ", closed_r.Qy_frac)
    shift = TARGET_BM_ANGLE - layout.bm_angle
    println("angle shift = ", shift, "   16 * shift = ", 16 * shift)
    println("bend delta = ", closed_r.total_bend_angle - open_r.total_bend_angle)
    println("outputs:")
    for path in (report_md, elements_csv, sequence_csv, closure_csv, figure, written)
        println("  ", path)
    end
    flush(stdout)
    return nothing
end

function _check_tune_closure(r, label)
    limit = 1e-8
    if abs(r.closure_Qx) > limit || abs(r.closure_Qy) > limit
        error("tune readings disagree on the $label ring: " *
              "closure_Qx=$(r.closure_Qx), closure_Qy=$(r.closure_Qy). " *
              "Accumulated Qx=$(r.Qx), Qy=$(r.Qy); fractional Qx=$(r.Qx_frac), Qy=$(r.Qy_frac).")
    end
    return nothing
end

function _write_structure_roundtrip(closed::Lattice)
    tmp = joinpath(OUT_DIR, "_write_structure_roundtrip.bell")
    mkpath(OUT_DIR)
    write_structure(tmp, closed)
    back = load_structure(tmp)
    rm(tmp)
    bm_angle_delta = 0.0
    n_bm = 0
    n_other = length(closed.elements) == length(back.elements) ? 0 : -1
    for (a, b) in zip(closed.elements, back.elements)
        if a.name == "BM" && a isa SBend && b isa SBend
            n_bm += 1
            bm_angle_delta = b.angle - a.angle
            restored = SBend(a.name, a.L, b.angle, a.k1, a.e1, a.e2, a.gap, a.fint)
            isequal(restored, b) || (n_other += 1)
        elseif !isequal(a, b)
            n_other += 1
        end
    end
    return (
        match = same_elements(closed, back),
        bend_delta = total_bend_angle(back) - total_bend_angle(closed),
        path_delta = path_length(back) - path_length(closed),
        bm_angle_delta = bm_angle_delta,
        n_bm = n_bm,
        n_other = n_other,
        ref_match = isequal(closed.ref, back.ref),
    )
end

function _blank(x)
    return x === nothing ? "" : string(x)
end

function _write_unique_csv(path, source, closed)
    open_rows = unique_rows(source)
    closed_rows = unique_rows(closed)
    length(open_rows) == length(closed_rows) || error("unique-element count changed")
    open(path, "w") do io
        println(io, "name,count,type,length_m,angle_open_rad,angle_closed_rad,k1,k2")
        for (a, b) in zip(open_rows, closed_rows)
            a.name == b.name && a.count == b.count && a.type == b.type && a.length == b.length &&
                a.k1 == b.k1 && a.k2 == b.k2 ||
                error("unique element $(a.name) changed beyond the bend angle")
            println(io, join((
                a.name, a.count, a.type, a.length, _blank(a.angle), _blank(b.angle),
                _blank(a.k1), _blank(a.k2),
            ), ","))
        end
    end
    return path
end

function _write_sequence_csv(path, source, closed)
    length(source.elements) == length(closed.elements) || error("sequence length changed")
    open(path, "w") do io
        println(io, "index,name,type,length_m,angle_open_rad,angle_closed_rad,k1,k2")
        for (i, (a, b)) in enumerate(zip(source.elements, closed.elements))
            ra = element_parameters(a)
            rb = element_parameters(b)
            ra.name == rb.name && ra.type == rb.type && ra.length == rb.length &&
                ra.k1 == rb.k1 && ra.k2 == rb.k2 ||
                error("sequence element $i changed beyond the bend angle")
            println(io, join((
                i, ra.name, ra.type, ra.length, _blank(ra.angle), _blank(rb.angle),
                _blank(ra.k1), _blank(ra.k2),
            ), ","))
        end
    end
    return path
end

function _closure_rows(layout, open_r, closed_r)
    shift = TARGET_BM_ANGLE - layout.bm_angle
    return [
        ("bm_angle_rad", layout.bm_angle, TARGET_BM_ANGLE, shift),
        ("total_bend_angle_rad", open_r.total_bend_angle, closed_r.total_bend_angle,
         closed_r.total_bend_angle - open_r.total_bend_angle),
        ("total_bend_minus_2pi_rad", open_r.bend_minus_2pi, closed_r.bend_minus_2pi,
         closed_r.bend_minus_2pi - open_r.bend_minus_2pi),
        ("path_length_m", open_r.path_length, closed_r.path_length,
         closed_r.path_length - open_r.path_length),
        ("Qx", open_r.Qx, closed_r.Qx, closed_r.Qx - open_r.Qx),
        ("Qy", open_r.Qy, closed_r.Qy, closed_r.Qy - open_r.Qy),
        ("Qx_frac", open_r.Qx_frac, closed_r.Qx_frac, closed_r.Qx_frac - open_r.Qx_frac),
        ("Qy_frac", open_r.Qy_frac, closed_r.Qy_frac, closed_r.Qy_frac - open_r.Qy_frac),
    ]
end

function _write_closure_csv(path, rows)
    open(path, "w") do io
        println(io, "quantity,before,after,delta")
        for (name, before, after, delta) in rows
            println(io, name, ",", before, ",", after, ",", delta)
        end
    end
    return path
end

function _md_table(headers, rows)
    esc(s) = replace(string(s), "|" => "\\|")
    head = "| " * join(esc.(headers), " | ") * " |"
    sep = "| " * join(["---" for _ in headers], " | ") * " |"
    body = ["| " * join(esc.(r), " | ") * " |" for r in rows]
    return join(vcat(head, sep, body), "\n")
end

function _write_report(path, layout, token, source, closed, open_r, closed_r, rows,
                       written, roundtrip_note)
    unique_body = Any[]
    for (a, b) in zip(unique_rows(source), unique_rows(closed))
        push!(unique_body, [
            a.name, a.count, a.type, a.length, _blank(a.angle), _blank(b.angle),
            _blank(a.k1), _blank(a.k2),
        ])
    end
    sequence_body = Any[]
    for (i, (a, b)) in enumerate(zip(source.elements, closed.elements))
        ra = element_parameters(a)
        rb = element_parameters(b)
        push!(sequence_body, [
            i, ra.name, ra.type, ra.length, _blank(ra.angle), _blank(rb.angle),
            _blank(ra.k1), _blank(ra.k2),
        ])
    end
    closure_body = [[name, before, after, delta] for (name, before, after, delta) in rows]
    io = IOBuffer()
    println(io, "# magnetic_2 lattice closed to 2π")
    println(io)
    println(io, "Source file (read-only): `", BELL_LATTICE, "`.")
    println(io)
    println(io, "The file token on the BM card is `", token, "`. Loaded as a Float64 it is ",
            layout.bm_angle, ". Sixteen equal sector bends therefore sum to ",
            open_r.total_bend_angle, ", which is ", open_r.bend_minus_2pi,
            " relative to 2π. `OBM` is a drift of the same length and is left unchanged. Every BM angle in the closed ring is the Float64 value of π/8, ",
            TARGET_BM_ANGLE, ". `total_bend_angle` of that ring is ",
            closed_r.total_bend_angle, ", which differs from 2π by ",
            closed_r.bend_minus_2pi, ". The product 16 * (π/8) differs from 2π by ",
            16 * TARGET_BM_ANGLE - 2π, ".")
    println(io)
    println(io, "Closed lattice: `", written, "`. Reloading that file reproduces the in-memory closed ring, including the reference particle. A separate `write_structure` round trip does not: element-wise match ",
            roundtrip_note.match, ", BM angle change ", roundtrip_note.bm_angle_delta,
            " rad on each of ", roundtrip_note.n_bm, " bends, total_bend_angle change ",
            roundtrip_note.bend_delta, " rad. That bend change is not 16 times the angle change: total_bend_angle reduces over the whole element list. Path-length change ", roundtrip_note.path_delta,
            " m, other element mismatches ", roundtrip_note.n_other,
            ", reference match ", roundtrip_note.ref_match,
            ". The saved lattice is the edited copy, not the `write_structure` text.")
    println(io)
    println(io, "Counts: ", layout.n_bm, " BM, ", layout.n_obm,
            " OBM, shared length ", layout.length_m, " m, ", layout.n_elements,
            " elements in sequence. Reference: ", source.ref.particle.name,
            ", kinetic energy ", source.ref.Ekin_MeV, " MeV, harmonic ", source.harmonic, ".")
    println(io)
    println(io, "`Qx` and `Qy` are accumulated tunes from `optics_along` (integer part included). `Qx_frac` and `Qy_frac` are μ/2π from `periodic_optics` of `ring_matrix`. Open-ring fractional closure residuals: Qx ",
            open_r.closure_Qx, ", Qy ", open_r.closure_Qy,
            ". Closed-ring residuals: Qx ", closed_r.closure_Qx, ", Qy ", closed_r.closure_Qy,
            ". Symplectic error of the one-turn matrix: open ", open_r.symplectic_error,
            ", closed ", closed_r.symplectic_error, ".")
    println(io)
    println(io, "Path length is the sum of element lengths. The arc length of BM is not changed, so the circumference does not change.")
    println(io)
    println(io, "## Closure")
    println(io)
    println(io, _md_table(["quantity", "before", "after", "delta"], closure_body))
    println(io)
    println(io, "## Unique elements")
    println(io)
    println(io, _md_table(["name", "count", "type", "length_m", "angle_open_rad", "angle_closed_rad", "k1", "k2"], unique_body))
    println(io)
    println(io, "## Sequence")
    println(io)
    println(io, _md_table(["index", "name", "type", "length_m", "angle_open_rad", "angle_closed_rad", "k1", "k2"], sequence_body))
    println(io)
    write(path, String(take!(io)))
    return path
end

function _write_deficit_svg(path, deficit_before, deficit_after)
    vals = [deficit_before, deficit_after]
    peak = max(maximum(abs.(vals)), eps(Float64))
    width = 640
    height = 360
    left = 90
    right = 40
    top = 70
    bottom = 70
    plot_w = width - left - right
    plot_h = height - top - bottom
    baseline = top + plot_h / 2
    scale = (plot_h / 2 - 16) / peak
    labels = ["before", "after"]
    colors = ["#9a3412", "#1d4ed8"]
    bars = String[]
    for i in 1:2
        cx = left + (i - 0.5) * plot_w / 2
        bh = vals[i] * scale
        y = bh >= 0 ? baseline - bh : baseline
        h = max(abs(bh), vals[i] == 0 ? 0 : 1.5)
        push!(bars, """<rect x="$(cx - 36)" y="$y" width="72" height="$h" fill="$(colors[i])"/>""")
        push!(bars, """<text x="$cx" y="$(baseline + plot_h / 2 + 22)" text-anchor="middle" font-size="14">$(labels[i])</text>""")
        push!(bars, """<text x="$cx" y="$(y - 8)" text-anchor="middle" font-size="13">$(@sprintf("%.3e", vals[i]))</text>""")
    end
    zero_label = @sprintf("0  (peak %.3e rad)", peak)
    svg = """
    <svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">
      <rect width="$width" height="$height" fill="#ffffff"/>
      <text x="32" y="36" font-family="sans-serif" font-size="18">Sector-bend deficit, total angle minus 2π</text>
      <line x1="$left" y1="$baseline" x2="$(width - right)" y2="$baseline" stroke="#334155" stroke-width="1.5"/>
      <text x="$(left - 8)" y="$(baseline + 4)" text-anchor="end" font-family="sans-serif" font-size="12">$zero_label</text>
      $(join(bars, "\n  "))
    </svg>
    """
    write(path, svg)
    return path
end

main()
