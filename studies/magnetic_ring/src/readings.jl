# Readings for the magnetic one-turn factor and for a thin RF cavity.
#
# Slip coefficients come from path_expansion. The second-order slip in the
# notes is η1 = α1 − α0/γ² + (3γ² − 1)/(2γ⁴). The normal-form actions of a
# laboratory point are J_k = (w_{2k−1}² + w_{2k}²) / 2 after
# u = A⁻¹ (z − z*) and w = T⁻¹(u). Courant–Snyder J at one location uses the
# on-momentum Twiss of that location and does not remove dispersion.

module LieReadings

using Bell

export eta1_formula, read_linear_regime, append_thin_cavity,
       normal_form_actions, normal_form_angles, courant_snyder_action,
       dispersion_free_transverse, track_turn_exact

"""
    eta1_formula(α0, α1, γ) -> η1

Second-order phase slip from the notes,

    η1 = α1 − α0/γ² + (3γ² − 1)/(2γ⁴).
"""
function eta1_formula(alpha0::Real, alpha1::Real, gamma::Real)
    gamma > 1 || error("gamma must be greater than 1, got $gamma")
    g2 = float(gamma)^2
    return float(alpha1) - float(alpha0) / g2 + (3 * g2 - 1) / (2 * g2^2)
end

"""
    read_linear_regime(path) -> NamedTuple

`voltage_MV`, `harmonic`, `delta_m`, `Qs`, and `phi_s_rad` of the `linear`
object in the step-2 regime file. `phi_s_rad` there is the notes' stationary
phase. The cavity built for BELL uses `π` below transition, not this value.
"""
function read_linear_regime(path::AbstractString)
    isfile(path) || error("regime file is missing: $path")
    text = read(path, String)
    block = _object_after(text, "linear")
    return (
        voltage_MV = _json_number(block, "voltage_MV"),
        harmonic = Int(_json_number(block, "harmonic")),
        delta_m = _json_number(block, "delta_m"),
        Qs = _json_number(block, "Qs"),
        phi_s_rad = _json_number(block, "phi_s_rad"),
    )
end

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

"""
    append_thin_cavity(lattice, voltage_MV, harmonic, phi_s) -> (lattice, cavity)

Copy `lattice` and append one thin `RFCavity`. The ring harmonic is unchanged.
`phi_s` is the phase stored on the cavity.
"""
function append_thin_cavity(lattice::Lattice, voltage_MV::Real, harmonic::Integer, phi_s::Real)
    voltage_MV > 0 || error("voltage must be positive, got $voltage_MV MV")
    harmonic > 0 || error("harmonic must be positive, got $harmonic")
    cav = RFCavity("RF", voltage_MV, Int(harmonic), float(phi_s))
    cav.L == 0 || error("appended cavity is not thin: L = $(cav.L)")
    elems = AbstractElement[e for e in lattice.elements]
    push!(elems, cav)
    lat = Lattice(lattice.name * "+RF", elems, lattice.ref, lattice.harmonic)
    return lat, cav
end

"""
    normal_form_actions(nf, z) -> (Jx, Jy, Js)

Actions of a canonical point. `u = A⁻¹ (z − z*)`, then `w = Tinv(u)`,
and `J_k = (w_{2k−1}² + w_{2k}²) / 2`. The angle on that plane is
`atan2(−w_{2k}, w_{2k−1})`.
"""
function normal_form_actions(nf, z)
    u = nf.A \ (collect(Float64, z) .- nf.zstar)
    w = ntuple(i -> nf.Tinv[i](u), 6)
    return ntuple(k -> (w[2k - 1]^2 + w[2k]^2) / 2, 3)
end

"""
    normal_form_angles(nf, z) -> (θx, θy, θs)

Angles paired with [`normal_form_actions`](@ref):
`w_{2k−1} = √(2J) cos θ`, `w_{2k} = −√(2J) sin θ`.
"""
function normal_form_angles(nf, z)
    u = nf.A \ (collect(Float64, z) .- nf.zstar)
    w = ntuple(i -> nf.Tinv[i](u), 6)
    return ntuple(k -> atan(-w[2k], w[2k - 1]), 3)
end

"""
    courant_snyder_action(x, px, β, α) -> J

Courant–Snyder action at one location,

    J = (γ x² + 2 α x p_x + β p_x²) / 2,   γ = (1 + α²) / β.

`x` and `p_x` are the canonical coordinates at that location. Dispersion is
not removed.
"""
function courant_snyder_action(x::Real, px::Real, beta::Real, alpha::Real)
    beta > 0 || error("beta must be positive, got $beta")
    gamma = (1 + float(alpha)^2) / float(beta)
    return (gamma * x^2 + 2 * float(alpha) * x * px + float(beta) * px^2) / 2
end

"""Betatron coordinates at one location: `(x − D δ, p_x − D' δ)`."""
function dispersion_free_transverse(x::Real, px::Real, D::Real, Dp::Real, delta::Real)
    return (x - D * delta, px - Dp * delta)
end

"""
    track_turn_exact(z, lattice, ref, path_length; ds = 0.1) -> Vector

One turn. Every magnet is `exact_track_element`. A thin cavity is `rf_kick`.
"""
function track_turn_exact(z, lattice::Lattice, ref::RefParticle, path_length::Real; ds::Real = 0.1)
    w = collect(Float64, z)
    for e in lattice.elements
        if e isa RFCavity
            e.L == 0 || error("track_turn_exact expected a thin cavity, L = $(e.L)")
            w = rf_kick(w, e, ref, path_length)
        else
            w = exact_track_element(w, e, ref, path_length; ds = ds)
        end
    end
    return w
end

end # module
