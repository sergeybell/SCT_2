# =============================================================================
# EquilibriumStudy.jl - δ_eq: formula, Lie product, truncated map, exact track
# =============================================================================
# Reported phi is the NOTE convention (phi>0 = LATE) = − Bell phi.
# Tracking uses Bell canonical variables; conversion only for reporting.
# The thin RF cavity is an element of the lattice, not a kick after the turn.
# =============================================================================
module EquilibriumStudy

using Bell, LinearAlgebra, Statistics, Printf
# TaylorSeries is Bell's backend; evaluate() is used on NF polynomials.
using TaylorSeries

export load_magnetic, setup_thin_rf_ring, ThinRFStudy,
       note_phi, lie_phi, measure_synchrotron,
       formula_A, formula_B_Lambda, formula_C, formula_C_factor, formula_C_Omega,
       formula_delta_eq, action_from_amp, action_J, cs_action, cs_J_turn, dE_eq,
       init_canonical, step_exact, step_map, step_lie, run_turns,
       deltas_of, phis_note_of, running_mean, one_turn_lie, sample_along,
       nf_angles, setup_ring, track_longitudinal

"Bell structure written from the OptiM file on each load."
const STRUCTURE_FILE = normpath(joinpath(@__DIR__, "..", "lattices", "magnetic_2.bell"))

"Bell phi (early+) -> note phi (late+)."
note_phi(phi_lie) = -phi_lie
lie_phi(phi_note) = -phi_note

"""
    load_magnetic(optim_path; structure_path, Ekin, h) -> (lat, ref, of)

Convert an OptiM file to a Bell structure, load it, then set the study
kinetic energy and the harmonic used for φ. The file written by the
converter keeps harmonic 1; the returned lattice uses `h`.
"""
function load_magnetic(optim_path; structure_path = STRUCTURE_FILE, Ekin = 270.0, h = 5)
    mkpath(dirname(structure_path))
    of = convert_optim(optim_path, structure_path)
    loaded = load_structure(structure_path)
    ref = reference_from_kinetic(DEUTERON, Ekin)
    lat = Lattice(loaded.name, loaded.elements, ref, h)
    return lat, ref, of
end

"""
One-turn Lie factor of a lattice that already contains a thin RFCavity.
Magnetic elements go through Bell `build_factors`; the cavity becomes a
Dragt–Finn factor of `rf_kick_taylor`. Do not call `one_turn_factor_rf`
on this lattice: that would add a second kick.
"""
function one_turn_lie(lat::Lattice, ref::RefParticle, C::Real)
    parts = LieFactor[]
    for F in build_factors(lat)
        if F isa RFCavity
            push!(parts, dragt_finn_factor(rf_kick_taylor(F, ref, C); name = F.name))
        else
            push!(parts, F)
        end
    end
    return product_factors(parts; name = lat.name)
end

"Ring with a thin RFCavity at the end of the element list, plus the maps used by M1–M4."
struct ThinRFStudy
    lat::Lattice
    lat_mag::Lattice
    cav::RFCavity
    ref::RefParticle
    C::Float64
    h::Int
    V::Float64
    pe::NamedTuple
    xi::NamedTuple
    po::NamedTuple
    F::LieFactor
    nf::NF6D
    M::Vector
    factors::Vector
    Qs::Float64
    of::Any
end

"""
Magnetic lattice plus a thin `RFCavity` in `lat.elements`. Voltage is on.
η, ξ come from the magnetic elements of the same list: the thin kick does
not lengthen the orbit. The cavity is present and fires every turn.
"""
function setup_thin_rf_ring(optim_path; Ekin = 270.0, V = 0.2, h = 5, phi_s = Float64(pi))
    lat_mag, ref, of = load_magnetic(optim_path; Ekin = Ekin, h = h)
    C = circumference(lat_mag)
    cav = RFCavity("RF", V, h, phi_s)
    elems = copy(lat_mag.elements)
    push!(elems, cav)
    lat = Lattice(lat_mag.name * "+RF", elems, ref, h)
    Fmag = ring_lie_factor(lat_mag)
    Mmag = factor_to_map(Fmag)
    pe = path_expansion(Mmag, ref, C)
    xi = chromaticity_orders(Mmag)
    po = periodic_optics(ring_matrix(lat_mag), ref, C)
    factors = collect(build_factors(lat))
    F = one_turn_lie(lat, ref, C)
    nf = dragt_finn_nf6d(F)
    M = factor_to_map(F)
    Qs = abs(nf.mu[3] / 2pi)
    return ThinRFStudy(lat, lat_mag, cav, ref, C, h, V, pe, xi, po, F, nf, M, factors, Qs, of)
end

# --- one-turn steps (Bell only) ---------------------------------------------

"M4: exact integrator of each magnetic Hamiltonian, then `rf_kick` at the cavity."
function step_exact(S::ThinRFStudy, z; ds = 0.1)
    w = collect(z)
    for e in S.lat.elements
        w = e isa RFCavity ? rf_kick(w, e, S.ref, S.C) :
            exact_track_element(w, e, S.ref, S.C; ds = ds)
    end
    return w
end

"M3: degree-3 polynomial of the same Lie product, `evaluate_map`."
step_map(S::ThinRFStudy, z) = evaluate_map(S.M, z)

"One turn of the Lie product itself (`track` / `apply_factor` + `rf_kick`)."
step_lie(S::ThinRFStudy, z) = track(S.factors, z; turns = 1, ref = S.ref, C = S.C)

function run_turns(step, z0, N)
    Zs = Matrix{Float64}(undef, 6, N)
    z = collect(z0)
    for n in 1:N
        z = step(z)
        Zs[:, n] = z
    end
    return Zs
end

deltas_of(Zs) = Zs[6, :]
phis_note_of(Zs, S::ThinRFStudy) =
    [note_phi(canonical_to_user(Zs[:, n], S.ref, S.C, S.h)[5]) for n in 1:size(Zs, 2)]

running_mean(v) = cumsum(v) ./ (1:length(v))

"Courant–Snyder action of a particle that starts at (x=A, x'=0)."
cs_action(A, p) = (1 + p.alpha^2) * A^2 / (2 * p.beta)

"Turn-by-turn Courant–Snyder action from canonical (x, px) at the observation point."
function cs_J_turn(z, p)
    x, px = z[1], z[2]
    γ = (1 + p.alpha^2) / p.beta
    return (γ * x^2 + 2 * p.alpha * x * px + p.beta * px^2) / 2
end

"""
Walk the lattice element by element. Thin RFCavity is applied with `rf_kick`
and does not advance s. Returns s, x, x', path excess ΔL = −l.
"""
function sample_along(S::ThinRFStudy, z0; nturns = 2)
    z = collect(z0)
    s = Float64[0.0]
    x = Float64[z[1]]
    xp = Float64[canonical_to_user(z, S.ref, S.C, S.h)[2]]
    L = Float64[-z[5]]
    st = 0.0
    for _ in 1:nturns
        for e in S.lat.elements
            if e isa RFCavity
                z = rf_kick(z, e, S.ref, S.C)
            else
                z = exact_track_element(z, e, S.ref, S.C)
                st += elem_length(e)
            end
            u = canonical_to_user(z, S.ref, S.C, S.h)
            push!(s, st); push!(x, z[1]); push!(xp, u[2]); push!(L, -z[5])
        end
    end
    return (s = s, x = x, xp = xp, L = L, z = z)
end

function init_canonical(S::ThinRFStudy, dm, Ax, Ay)
    u = [S.po.D4[1] * dm + Ax, S.po.D4[2] * dm, Ay, 0.0, 0.0, dm]
    return user_to_canonical(u, S.ref, S.C, S.h)
end

"""
Normal-form angles and actions of a laboratory point.
w_{2k-1}=√(2J) cos θ, w_{2k}=−√(2J) sin θ, so θ = atan2(−w_{2k}, w_{2k-1}).
"""
function nf_angles(nf::NF6D, z)
    u = nf.A \ (collect(z) .- nf.zstar)
    w = [Float64(evaluate(nf.Tinv[i], u)) for i in 1:6]
    θ = ntuple(k -> atan(-w[2k], w[2k - 1]), 3)
    J = ntuple(k -> (w[2k - 1]^2 + w[2k]^2) / 2, 3)
    return θ, J
end

# --- formulas from the note (ch. 1–3) ----------------------------------------
formula_A(eta0, eta1, dm) = -eta1 * dm^2 / (2 * eta0)
formula_B_Lambda(C, Jx, xi_x, Jy, xi_y) = -2pi / C * (Jx * xi_x + Jy * xi_y)
formula_C_factor(phat) = 1 - phat^2 / 32
formula_C_Omega(phat) = 1 - phat^2 / 16
formula_C(eta0, eta1, dm, phat) = formula_A(eta0, eta1, dm) * formula_C_factor(phat)

function formula_delta_eq(pe, xi, C, dm, Jx, Jy; phat = 0.0)
    Λ = formula_B_Lambda(C, Jx, xi.x[2], Jy, xi.y[2])
    return -(pe.eta1 * dm^2 / 2 * formula_C_factor(phat) + Λ) / pe.eta0
end

"Courant-Snyder action J = A^2 / (2 beta) for a particle at crest (alpha~0)."
action_from_amp(A, beta) = A^2 / (2 * beta)
action_J(A, beta) = action_from_amp(A, beta)

function dE_eq(ref, deq, dm; phat = 0.0)
    fac = phat == 0 ? 1.0 : formula_C_factor(phat)
    beta2 = ref.beta^2; Es = ref.gamma * ref.mass_MeV
    beta2 * Es * (deq + dm^2 / (4 * ref.gamma^2) * fac)
end

"""
Measure synchrotron stats via harmonic least squares:
  delta_n = deq + a cos(2pi Q n) + b sin + c cos(4pi Q n) + d sin
Q refined on a grid. Removes O(dm/N_periods) window leakage.
"""
function measure_synchrotron(deltas, phis; n_drop = nothing, Qs = 0.01)
    N = length(deltas)
    nd = n_drop === nothing ? max(100, round(Int, 5 / max(Qs, 1e-4))) : n_drop
    nd = min(nd, N ÷ 4)
    d = deltas[nd+1:end]; p = phis[nd+1:end]
    n = length(d)
    Q0 = try; measure_tune(d); catch; abs(Qs); end
    (isnan(Q0) || Q0 == 0) && (Q0 = abs(Qs))
    function fit_at(Q)
        t = float.(0:n-1)
        A = hcat(ones(n), cos.(2pi .* Q .* t), sin.(2pi .* Q .* t),
                 cos.(4pi .* Q .* t), sin.(4pi .* Q .* t))
        coef = A \ d
        resid = d .- A * coef
        return coef, sqrt(mean(abs2, resid))
    end
    bestQ, bestrms = Q0, Inf
    for Q in range(Q0 - 0.002, Q0 + 0.002; length = 81)
        _, rms = fit_at(Q)
        rms < bestrms && ((bestrms, bestQ) = (rms, Q))
    end
    coef, _ = fit_at(bestQ)
    deq, a, b, c, dd = coef[1], coef[2], coef[3], coef[4], coef[5]
    dm = sqrt(a^2 + b^2)
    dm2h = sqrt(c^2 + dd^2)
    dm_hr = (maximum(d) - minimum(d)) / 2
    d2 = mean(abs2, d .- deq)
    phat = (maximum(p) - minimum(p)) / 2
    n2 = n ÷ 2
    function deq_of(seg)
        t = float.(0:length(seg)-1)
        A = hcat(ones(length(seg)), cos.(2pi .* bestQ .* t), sin.(2pi .* bestQ .* t),
                 cos.(4pi .* bestQ .* t), sin.(4pi .* bestQ .* t))
        return (A \ seg)[1]
    end
    deq1, deq2 = deq_of(d[1:n2]), deq_of(d[n2+1:end])
    return (delta_eq = deq, delta_m = dm_hr, delta_m_fit = dm, delta2 = d2,
            delta2_raw = mean(abs2, d), phi_hat = phat, Qs_meas = bestQ,
            n_used = n, n_drop = nd, n_periods = round(Int, n * bestQ),
            Tsyn = max(1, round(Int, 1 / bestQ)),
            err_split = abs(deq1 - deq2) / 2, delta_eq_h1 = deq1, delta_eq_h2 = deq2,
            rms = bestrms, a1 = dm, a2 = dm2h)
end

# --- leftover helpers (finite cavity; not used by the thin-RF study) ---------

"""
Lattice magnetic_2 + 1 m cavity, stationary bucket phi_s = pi (Bell),
below transition (eta0 < 0). Returns (lat, ref, C, h, V, T, pe, xi, po, info).
"""
function setup_ring(optim_path; Ekin = 270.0, V = 0.2, Lc = 1.0, h = 5)
    lat0, ref, of = load_magnetic(optim_path; Ekin = Ekin, h = h)
    cav = Cavity("RF", Lc, V, h, Float64(pi))
    lat, icav = insert_cavity(lat0, cav)
    lat = Lattice(lat.name, lat.elements, ref, h)
    C = circumference(lat)
    Ttf = transit_factor(cav, ref, C)
    lat_mag = Lattice(lat.name, [e isa Cavity ? Drift(e.name, e.L) : e for e in lat.elements], ref, h)
    M = ring_taylor_map(lat_mag)
    pe = path_expansion(M, ref, C)
    xi = chromaticity_orders(M)
    po = periodic_optics(ring_matrix(lat_mag), ref, C)
    Qs = synchrotron_tune(ref, pe.eta0, V, h, Float64(pi); T = Ttf)
    info = (Qs = Qs, Lc = Lc, icav = icav, mass_file = of.mass_MeV,
            nperiods_file = of.nperiods_file, nrepeat_used = of.nrepeat_used,
            warnings = of.warnings)
    return lat, ref, C, h, V, Ttf, pe, xi, po, info
end

"""
Track `turns` revolutions. A thin RFCavity is applied with `rf_kick`;
other elements go through the exact integrator.
"""
function track_longitudinal(lat, ref, C, h, u0; turns = 1000, ds = 0.1)
    z = user_to_canonical(u0, ref, C, h)
    deltas = Vector{Float64}(undef, turns)
    phis = Vector{Float64}(undef, turns)
    for n in 1:turns
        for e in lat.elements
            z = e isa RFCavity ? rf_kick(z, e, ref, C) :
                exact_track_element(z, e, ref, C; ds = ds)
        end
        u = canonical_to_user(z, ref, C, h)
        deltas[n] = u[6]
        phis[n] = note_phi(u[5])
    end
    return deltas, phis
end

end # module
