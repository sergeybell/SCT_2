# Magnetic ring study on the current BELL API.
# The BELL tree is read-only. This module loads examples/lattices/magnetic_2.bell
# and writes a closed copy under this study.

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

export BELL_LATTICE, TARGET_BM_ANGLE, closed_lattice_path,
       load_magnetic_2, close_bm_angles, write_closed_lattice, load_closed_lattice,
       element_parameters, unique_rows, source_layout, only_bm_angle_changed,
       same_elements, ring_readings, source_bm_angle_token

end # module
