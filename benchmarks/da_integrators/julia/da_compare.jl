# da_compare.jl - (1) coefficient agreement Bell.DiffAlgebra vs COSY DA, (2) timing of DA operations.
# Needs cosy/output/acc_o<N>.txt, time_o<N>.txt (cosy/run_cosy_da.ps1).
# Run: julia --project=julia julia/da_compare.jl      -> output/da_accuracy.csv, da_timing.csv, da_compile.csv
t_start = time()
using Bell, Bell.DiffAlgebra
t_load = time() - t_start
using BenchmarkTools, Printf
include(joinpath(@__DIR__, "cosy_io.jl"))
const ROOT = dirname(@__DIR__); const OUT = joinpath(ROOT, "output"); mkpath(OUT)
const ORDERS = [3, 5, 7, 10]
const NAMES = ["sqrt(1+v)", "exp(v)", "sin(0.3+v)", "1/(1+x1+x2)", "log(2+v)", "atan(0.2+v)",
               "cosh(v)tan(0.1+v)", "sqrt(1+v)(1+v)^2", "exp(v)sin(0.3+v)", "exp∘M", "sin(0.3+·)∘M",
               "M^-1_1", "M^-1_2", "M^-1_3", "M^-1_4", "M^-1_5", "M^-1_6"]

function problems(no)                # identical definitions to cosy/da_ops_template.fox
    sp = DASpace(6, no); X = da_vars(sp)
    v = X[1] + 0.5X[2] - 0.3X[3] + 0.2X[4] + 0.1X[5] - 0.4X[6]
    F = [sqrt(1 + v), exp(v), sin(0.3 + v), 1 / (1 + X[1] + X[2]), log(2 + v), atan(0.2 + v),
         cosh(v) * tan(0.1 + v), sqrt(1 + v) * (1 + v) * (1 + v)]
    M = [X[i] + 0.5X[mod(i, 6) + 1] + 0.2X[i] * X[mod(i + 2, 6) + 1] + 0.1sin(v) * X[mod(i, 6) + 1]^2 for i in 1:6]
    C = compose([F[2], F[3]], M)
    Mi = inverse_map(M)
    return sp, X, v, F, M, vcat(F, [F[2] * F[3]], C, Mi)
end

# compile / first-call times (fresh process: the first problems() call compiles everything)
tc = @elapsed problems(3); tr = @elapsed problems(3)
open(joinpath(OUT, "da_compile.csv"), "w") do io
    println(io, "item,seconds"); println(io, "using Bell (load),$t_load")
    println(io, "first problems(3) call incl. compilation,$tc"); println(io, "second call,$tr")
end

acc = open(joinpath(OUT, "da_accuracy.csv"), "w"); println(acc, "order,name,ncoef_cosy,ncoef_bell,max_abs_diff,max_abs_coef,max_rel_diff,max_rel_per_coef")
tim = open(joinpath(OUT, "da_timing.csv"), "w"); println(tim, "order,op,nmonomials,bell_s,cosy_s,cosy_overhead_s,ratio_bell_over_cosy,bell_allocs,bell_bytes")
for no in ORDERS
    sp, X, v, F, M, R = problems(no)
    blocks = read_cosy_da(joinpath(ROOT, "cosy", "output", "acc_o$no.txt"))
    @assert length(blocks) == length(R)
    for (k, (b, a)) in enumerate(zip(blocks, R))
        maxd = 0.0; maxc = maximum(abs, a.c); maxr = 0.0; seen = 0
        for (e, c) in b
            cb = coefficient(a, e); d = abs(cb - c); maxd = max(maxd, d)
            abs(c) > 1e-12 * maxc && (maxr = max(maxr, d / abs(c)))
            seen += 1
        end
        nb = count(!iszero, a.c)
        # Bell coefficients that COSY did not print (COSY drops |c| < eps): their size
        cset = Set(keys(b))
        for m in 1:nmonomials(sp)
            e = monomial_exponents(sp, m); a.c[m] != 0 && !(e in cset) && (maxd = max(maxd, abs(a.c[m])))
        end
        @printf(acc, "%d,%s,%d,%d,%.3e,%.3e,%.3e,%.3e\n", no, NAMES[k], length(b), nb, maxd, maxc, maxd / maxc, maxr)
    end
    # timing (same dense arguments as COSY: A = F3, B = F2)
    A = F[3]; B = F[2]
    ops = [("mul", () -> A * B), ("add", () -> A + B), ("exp", () -> exp(A)), ("sin", () -> sin(A)),
           ("sqrt", () -> sqrt(B)), ("inv", () -> 1 / B), ("compose6", () -> compose(M, M)),
           ("invmap6", () -> inverse_map(M))]
    ct = read_cosy_times(joinpath(ROOT, "cosy", "output", "time_o$no.txt"))
    for (name, f) in ops
        b = @benchmark $f() seconds = 3 evals = 1
        tb = minimum(b).time * 1e-9; tm = median(b).time * 1e-9
        tco = ct[name] - ct["assign"]
        @printf(tim, "%d,%s,%d,%.4e,%.4e,%.4e,%.3f,%d,%d\n", no, name, nmonomials(sp), tm, tco, ct["assign"], tm / tco, b.allocs, b.memory)
        @printf("o=%2d %-9s Bell %.3e s  COSY %.3e s  ratio %.2f  allocs %d\n", no, name, tm, tco, tm / tco, b.allocs)
    end
end
close(acc); close(tim)