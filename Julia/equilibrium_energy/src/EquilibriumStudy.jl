# =============================================================================
# EquilibriumStudy.jl - delta_eq measurements vs note formulas (ch. 1-3)
# =============================================================================
# Convention in reported phi: NOTE convention (phi>0 = LATE) = - Bell phi.
# Tracking uses Bell canonical variables; conversion only for reporting.
# =============================================================================
module EquilibriumStudy

using Bell, LinearAlgebra, Statistics, Printf

export setup_ring, note_phi, lie_phi, track_longitudinal, measure_synchrotron,
       formula_A, formula_B_Lambda, formula_C, formula_C_factor, formula_C_Omega,
       action_from_amp, action_J, dE_eq, load_magnetic

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
Lattice magnetic_2 + 1 m cavity, stationary bucket phi_s = pi (Bell),
below transition (eta0 < 0). Returns (lat, ref, C, h, V, T, pe, xi, po, info).
"""
function setup_ring(optim_path; Ekin = 270.0, V = 0.2, Lc = 1.0, h = 5)
    lat0, ref, of = load_magnetic(optim_path; Ekin = Ekin, h = h)
    cav = Cavity("RF", Lc, V, h, Float64(pi))   # stationary; below transition
    lat, icav = insert_cavity(lat0, cav)
    lat = Lattice(lat.name, lat.elements, ref, h)
    C = circumference(lat)
    Ttf = transit_factor(cav, ref, C)
    # RF off for optics / eta
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
Track `turns` revolutions with the exact integrator starting from user state
`u0 = (x, x', y, y', phi_lie, delta)`. Returns (deltas, note_phis) per turn.
"""
function track_longitudinal(lat, ref, C, h, u0; turns = 1000, ds = 0.1)
    z = user_to_canonical(u0, ref, C, h)
    deltas = Vector{Float64}(undef, turns)
    phis = Vector{Float64}(undef, turns)
    for n in 1:turns
        z = exact_track_turn(z, lat, ref, C; ds = ds)
        u = canonical_to_user(z, ref, C, h)
        deltas[n] = u[6]
        phis[n] = note_phi(u[5])
    end
    return deltas, phis
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
            n_used = n, n_periods = round(Int, n * bestQ), Tsyn = max(1, round(Int, 1 / bestQ)),
            err_split = abs(deq1 - deq2) / 2, delta_eq_h1 = deq1, delta_eq_h2 = deq2,
            rms = bestrms, a1 = dm, a2 = dm2h)
end

# --- Formulas from the note (ch.1-3) ---
formula_A(eta0, eta1, dm) = -eta1 * dm^2 / (2 * eta0)
formula_B_Lambda(C, Jx, xi_x, Jy, xi_y) = -2pi / C * (Jx * xi_x + Jy * xi_y)
formula_C_factor(phat) = 1 - phat^2 / 32
formula_C_Omega(phat) = 1 - phat^2 / 16
formula_C(eta0, eta1, dm, phat) = formula_A(eta0, eta1, dm) * formula_C_factor(phat)

"Courant-Snyder action J = A^2 / (2 beta) for a particle at crest (alpha~0)."
action_from_amp(A, beta) = A^2 / (2 * beta)
action_J(A, beta) = action_from_amp(A, beta)

function dE_eq(ref, deq, dm; phat = 0.0)
    fac = phat == 0 ? 1.0 : formula_C_factor(phat)
    beta2 = ref.beta^2; Es = ref.gamma * ref.mass_MeV
    beta2 * Es * (deq + dm^2 / (4 * ref.gamma^2) * fac)
end

end # module
