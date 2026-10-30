# profile.jl - (4) where does Bell spend time in a DA ODE integration?
#   quad (exact sqrt Hamiltonian), nv = 6, DOP853 tol 1e-12, orders 5, 7, 10
#   - @timed: total time, GC fraction, bytes allocated per RHS evaluation
#   - RHS timed on the dense final state (identity-like states are sparse and much cheaper!)
#   - cost model: RHS = 4 mul + 1 sqrt + 1 inv + scalar*DA/adds, measured separately
#   - Profile (flat, order 7): top frames -> output/profile_flat_o7.txt
# Run: julia --project=julia julia/profile.jl  -> output/profile.csv
using Bell, Bell.DiffAlgebra, Bell.Integrators, BenchmarkTools, Printf, Profile
const OUT = joinpath(dirname(@__DIR__), "output")
const ETA = 0.2131578052277212; const KQ = 10 / 2.149635688032703
io = open(joinpath(OUT, "profile.csv"), "w")
println(io, "order,total_s,gc_frac,bytes_per_rhs,nfev,t_rhs_s,t_mul_s,t_sqrt_s,t_inv_s,t_axpy_s,rhs_model_s,frac_rhs_in_total,frac_lincomb_in_total")
for no in (5, 7, 10)
    sp = DASpace(6, no); X = da_vars(sp); dK = X[6]
    P2 = (1 + dK) * ((1 + dK) * ETA + 2) / (ETA + 2)
    f(z, s) = (ip = inv(sqrt(P2 - z[2] * z[2] - z[4] * z[4])); [z[2] * ip, -KQ * z[1], z[4] * ip, KQ * z[3]])
    z0 = X[1:4]
    run() = integrate(f, z0, 0.0, 0.5; method = DormandPrince853(), tol = 1e-12)
    run(); st = @timed run(); r = st.value
    zd = r.z                                  # dense state (as inside the integration)
    a = sqrt(P2 - zd[2] * zd[2] - zd[4] * zd[4]); b = zd[1]
    trhs = minimum(@benchmark($f($zd, 0.0), seconds = 2)).time * 1e-9
    tm = minimum(@benchmark($a * $b, seconds = 1)).time * 1e-9
    ts = minimum(@benchmark(sqrt($a), seconds = 1)).time * 1e-9
    ti = minimum(@benchmark(inv($a), seconds = 1)).time * 1e-9
    tx = minimum(@benchmark($a + 0.3 * $b, seconds = 1)).time * 1e-9
    model = 4tm + ts + ti + 4tx
    bytes = @allocated f(zd, 0.0)
    # outside RHS (dense-state RHS cost subtracted): the RK linear combinations z + h*sum(a_ij k_j): 13 stages, up to 12 terms, 4 comps
    flin = 1 - r.nfev * trhs / st.time
    @printf(io, "%d,%.4e,%.3f,%d,%d,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.3f,%.3f\n", no, st.time, st.gctime / st.time, bytes,
            r.nfev, trhs, tm, ts, ti, tx, model, r.nfev * trhs / st.time, flin)
    println("order $no: total $(st.time) s, nfev $(r.nfev), rhs $(trhs) s, model $(model) s, gc $(st.gctime/st.time)")
    if no == 7
        Profile.clear(); @profile (for _ in 1:3; run(); end)
        open(joinpath(OUT, "profile_flat_o7.txt"), "w") do pio
            Profile.print(IOContext(pio, :displaysize => (200, 220)); format = :flat, sortedby = :count, mincount = 20)
        end
    end
end
close(io); print(read(joinpath(OUT, "profile.csv"), String))