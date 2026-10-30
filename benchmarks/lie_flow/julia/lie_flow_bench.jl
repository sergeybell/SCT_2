# lie_flow_bench.jl - Bell.LieFlow (exp(-L:H:) with the full Hamiltonian) vs COSY INFINITY,
# vs the force-form DA flow prototype and vs DOP853/GBS (benchmarks/da_integrators).
# Environment: the sibling ../da_integrators/julia (Bell dev + BenchmarkTools).
#   julia --project=../da_integrators/julia julia/lie_flow_bench.jl
# Output (output/): flow_vs_cosy.csv  flow_symp.csv  flow_symp_L.csv  flow_speed.csv  ring_df_vs_flow.csv
using Bell, Bell.DiffAlgebra, BenchmarkTools, Printf, LinearAlgebra
using Bell.LieAlgebra.TaylorSeries: getcoeff
const ROOT = dirname(@__DIR__); const OUT = joinpath(ROOT, "output"); const CO = joinpath(ROOT, "cosy", "output")
const DAI = joinpath(dirname(ROOT), "da_integrators", "julia")
include(joinpath(DAI, "cosy_io.jl")); include(joinpath(DAI, "da_fast.jl"))
const CHIM, ETA = parse.(Float64, split(read(joinpath(CO, "refpart.txt"), String)))
# reference particle with exactly COSY's K/mc^2 (proton, 200 MeV)
const REF = reference_from_kinetic(Particle("proton (COSY mass)", 200.0 / ETA, 1.0, 0.0), 200.0)
const ELS = [("drift", Drift("D", 0.5)), ("quad", Quadrupole("Q", 0.5, (0.5 / 0.05) / CHIM)),
             ("dipole", SBend("B", 2pi / 6, pi / 6)), ("sext", Sextupole("S", 0.3, 2 * 0.2 / 0.05^2 / CHIM))]
const ORDERS = [3, 5, 7, 10]
# COSY variable dK -> Bell delta = dp/p0
onepd(X) = sqrt((1 + X[6]) * ((1 + X[6]) * ETA + 2) / (ETA + 2))
cosy_map(file, sp) = [(d = Dict(b); DA(sp, [get(d, monomial_exponents(sp, m), 0.0) for m in 1:nmonomials(sp)])) for b in read_cosy_da(file; nv = 6)]
tmin(f; s = 2) = minimum(@benchmark($f(), seconds = s, evals = 1)).time * 1e-9

ctime(no, k) = (l = split(readlines(joinpath(CO, "elem_time_o$no.txt"))[k]); parse(Float64, l[3]) / parse(Float64, l[2]))
acc = open(joinpath(OUT, "flow_vs_cosy.csv"), "w"); println(acc, "element,order,row,max_abs_diff,max_abs_coef,rel")
sym = open(joinpath(OUT, "flow_symp.csv"), "w"); println(sym, "element,order,bell_flow,cosy,bell_flow_rel,cosy_rel")
spd = open(joinpath(OUT, "flow_speed.csv"), "w"); println(spd, "element,order,bell_flow_s,cosy_s,forceflow_s,ratio_flow_cosy,ratio_force_cosy")
for no in ORDERS
    sp = flow_space(no); X = da_vars(sp)
    sub = [X[1], X[2], X[3], X[4], X[5], onepd(X) - 1]
    for (k, (nm, e)) in enumerate(ELS)
        M = lie_flow_map(e, REF, no)
        C = cosy_map(joinpath(CO, "elem_o$(no)_$k.txt"), sp)
        Mc = compose(M, sub)                       # Bell map with delta expressed through COSY dK
        for i in 1:4
            d = coeffnorm(Mc[i] - C[i]); c = coeffnorm(C[i])
            @printf(acc, "%s,%d,%d,%.3e,%.3e,%.3e\n", nm, no, i, d, c, d / c)
        end
        sb = da_symplectic_error(M); sc = da_symplectic_error(C)
        nb = maximum(coeffnorm, M)^2; nc = maximum(coeffnorm, C)^2
        @printf(sym, "%s,%d,%.3e,%.3e,%.3e,%.3e\n", nm, no, sb, sc, sb / nb, sc / nc)
        tb = tmin(() -> lie_flow_map(e, REF, no))
        tf = NaN
        if nm in ("quad", "dipole", "drift", "sext")      # force form z' = f(z), 4 states, order no
            P = onepd(X); P2 = P^2
            if nm == "quad"; kq = e.k1
                F = (ip = inv(sqrt(P2 - X[2]^2 - X[4]^2)); [X[2] * ip, -kq * X[1], X[4] * ip, kq * X[3]])
            elseif nm == "dipole"; h = 0.5
                F = (ee = 1 + h * X[1]; ps = sqrt(P2 - X[2]^2 - X[4]^2); ip = ee * inv(ps); [X[2] * ip, h * (ps - ee), X[4] * ip, 0 * X[4]])
            elseif nm == "drift"
                F = (ip = inv(sqrt(P2 - X[2]^2 - X[4]^2)); [X[2] * ip, 0 * X[1], X[4] * ip, 0 * X[1]])
            else; k2 = e.k2
                F = (ip = inv(sqrt(P2 - X[2]^2 - X[4]^2)); [X[2] * ip, -k2 / 2 * (X[1]^2 - X[3]^2), X[4] * ip, k2 * X[1] * X[3]])
            end
            zf = daflow(F, X, elem_length(e); nstate = 4)
            tf = tmin(() -> daflow(F, X, elem_length(e); nstate = 4))
            # force-form map is in (x,a,y,b) with dK as parameter: compare with COSY as a cross-check
            @printf(acc, "%s,%d,force%d,%.3e,%.3e,%.3e\n", nm, no, 1, maximum(coeffnorm(zf[i] - C[i]) for i in 1:4), maximum(coeffnorm(C[i]) for i in 1:4),
                    maximum(coeffnorm(zf[i] - C[i]) for i in 1:4) / maximum(coeffnorm(C[i]) for i in 1:4))
        end
        tc = ctime(no, k)
        @printf(spd, "%s,%d,%.4e,%.4e,%.4e,%.3f,%.3f\n", nm, no, tb, tc, tf, tb / tc, tf / tc)
        flush(acc); flush(sym); flush(spd); println("$nm N=$no: flow $(tb) s, COSY $(tc) s")
    end
end
close(acc); close(sym); close(spd)

# symplecticity vs length (scaling and squaring), quadrupole k = 4.65, and timing
io = open(joinpath(OUT, "flow_symp_L.csv"), "w"); println(io, "order,L,symp_abs,symp_rel,time_s")
for no in ORDERS, L in (0.1, 0.25, 0.5, 1.0, 2.0)
    e = Quadrupole("Q", L, (0.5 / 0.05) / CHIM)
    M = lie_flow_map(e, REF, no); t = @elapsed lie_flow_map(e, REF, no)
    s = da_symplectic_error(M)
    @printf(io, "%d,%g,%.3e,%.3e,%.4e\n", no, L, s, s / maximum(coeffnorm, M)^2, t)
end
close(io)

# (d) ring: Dragt-Finn one-turn factor (order 3) vs flow route
lat0 = load_structure(joinpath(ROOT, "..", "..", "chromaticity", "BELL", "magnetic_2_md.bell"))
lat = Lattice(lat0.name, AbstractElement[e for e in lat0.elements if !(e isa RFCavity)], lat0.ref, lat0.harmonic)
sp3 = flow_space(3)
toda(p) = DA(sp3, [getcoeff(p, collect(monomial_exponents(sp3, m))) for m in 1:nmonomials(sp3)])
t_df = @elapsed (MDF = factor_to_map(ring_lie_factor(lat)))
t_cl = @elapsed (MLc = (let M = collect(Z); for e in lat.elements; M = compose_maps(element_lie_map(e, lat.ref), M); end; M end))
DF = toda.(MDF); LC = toda.(MLc)
rows = String[]
t3 = @elapsed (F3 = one_turn_map_flow(lat, 3; cache = false))
push!(rows, @sprintf("flow N=3,%.4e,%.3e,%.3e", t3, maximum(coeffnorm(F3[i] - DF[i]) for i in 1:6), da_symplectic_error(F3)))
for N in (5, 7)
    t = @elapsed (FN = one_turn_map_flow(lat, N; cache = false))
    F3N = [change_order(f, sp3) for f in FN]
    push!(rows, @sprintf("flow N=%d (truncated to 3),%.4e,%.3e,%.3e", N, t, maximum(coeffnorm(F3N[i] - DF[i]) for i in 1:6), da_symplectic_error(FN)))
end
push!(rows, @sprintf("Lie maps composed (compose_maps),%.4e,%.3e,%.3e", t_cl, maximum(coeffnorm(LC[i] - DF[i]) for i in 1:6), da_symplectic_error(LC)))
push!(rows, @sprintf("Dragt-Finn factor_to_map(ring_lie_factor),%.4e,0,%.3e", t_df, da_symplectic_error(DF)))
open(joinpath(OUT, "ring_df_vs_flow.csv"), "w") do io
    println(io, "route,time_s,max_abs_diff_vs_DF_order3,symp_err,max_abs_coef_DF"); foreach(r -> println(io, r * @sprintf(",%.4e", maximum(coeffnorm, DF))), rows)
end
println(join(rows, "\n")); println("max |coef| of DF map: ", maximum(coeffnorm, DF))
println("flow N=3 vs composed Lie maps: ", maximum(coeffnorm(F3[i] - LC[i]) for i in 1:6))