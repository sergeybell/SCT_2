# hh_rk4_same.jl - identical algorithm check: Bell RK4 vs COSY FOX RK4 (Henon-Heiles, nv = 4, T = 1)
#   same step counts n = 10..320, compares maps coefficient-wise and time per map.
# Run: julia --project=julia julia/hh_rk4_same.jl  -> output/hh_rk4_same.csv
using Bell, Bell.DiffAlgebra, Bell.Integrators, Printf
include(joinpath(@__DIR__, "cosy_io.jl"))
const ROOT = dirname(@__DIR__); const CO = joinpath(ROOT, "cosy", "output")
f(z, s) = [z[3], z[4], -z[1] - 2z[1] * z[2], -z[2] - z[1] * z[1] + z[2] * z[2]]
io = open(joinpath(ROOT, "output", "hh_rk4_same.csv"), "w"); println(io, "order,nsteps,bell_s,cosy_s,ratio,max_abs_diff,max_abs_coef")
for no in (3, 5, 7, 10)
    sp = DASpace(4, no); z0 = [c + x for (c, x) in zip([0.1, -0.2, 0.25, 0.1], da_vars(sp))]
    ct = Dict(parse(Int, split(l)[1]) => parse(Float64, split(l)[3]) / parse(Float64, split(l)[2]) for l in eachline(joinpath(CO, "hh_time_o$no.txt")))
    for n in (10, 20, 40, 80, 160, 320)
        r = integrate(f, z0, 0.0, 1.0; method = RK4(), nsteps = n).z
        t = minimum(@elapsed(integrate(f, z0, 0.0, 1.0; method = RK4(), nsteps = n)) for _ in 1:max(1, round(Int, 0.02 / (ct[n] + 1e-9))) + 2)
        cz = read_cosy_da(joinpath(CO, "hh_o$(no)_n$n.txt"); nv = 4)
        d = 0.0; mc = 0.0
        for i in 1:4, m in 1:nmonomials(sp)
            e = monomial_exponents(sp, m); c = get(cz[i], e, 0.0); d = max(d, abs(r[i].c[m] - c)); mc = max(mc, abs(c))
        end
        @printf(io, "%d,%d,%.4e,%.4e,%.3f,%.3e,%.3e\n", no, n, t, ct[n], t / ct[n], d, mc)
    end
end
close(io); print(read(joinpath(ROOT, "output", "hh_rk4_same.csv"), String))