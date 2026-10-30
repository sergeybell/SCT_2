# da_fast_bench.jl - check and time the prototype optimisations (da_fast.jl) against Bell.DiffAlgebra
# Run: julia --project=julia julia/da_fast_bench.jl   -> output/da_fast.csv
using Bell, Bell.DiffAlgebra, BenchmarkTools, Printf
include(joinpath(@__DIR__, "da_fast.jl"))
const OUT = joinpath(dirname(@__DIR__), "output")
io = open(joinpath(OUT, "da_fast.csv"), "w"); println(io, "order,op,bell_s,proto_s,speedup,max_abs_diff")
for no in (3, 5, 7, 10)
    sp = DASpace(6, no); X = da_vars(sp)
    v = X[1] + 0.5X[2] - 0.3X[3] + 0.2X[4] + 0.1X[5] - 0.4X[6]
    M = [X[i] + 0.5X[mod(i, 6) + 1] + 0.2X[i] * X[mod(i + 2, 6) + 1] + 0.1sin(v) * X[mod(i, 6) + 1]^2 for i in 1:6]
    A = sin(0.3 + v); B = exp(v); c = zeros(nmonomials(sp))
    d1 = coeffnorm(DA(sp, mul_to!(c, A, B)) - A * B)
    t1 = minimum(@benchmark($A * $B, seconds = 2)).time * 1e-9
    t2 = minimum(@benchmark(mul_to!($c, $A, $B), seconds = 2)).time * 1e-9
    @printf(io, "%d,mul,%.4e,%.4e,%.2f,%.2e\n", no, t1, t2, t1 / t2, d1)
    C1 = compose(M, M); C2 = compose_fast(M, M); d2 = maximum(coeffnorm(C1[i] - C2[i]) for i in 1:6)
    t1 = minimum(@benchmark(compose($M, $M), seconds = 3, evals = 1)).time * 1e-9
    t2 = minimum(@benchmark(compose_fast($M, $M), seconds = 3, evals = 1)).time * 1e-9
    @printf(io, "%d,compose6,%.4e,%.4e,%.2f,%.2e\n", no, t1, t2, t1 / t2, d2)
    I1 = inverse_map(M); I2 = inverse_fast(M); d3 = maximum(coeffnorm(I1[i] - I2[i]) for i in 1:6)
    t1 = no == 10 ? (@elapsed inverse_map(M)) : minimum(@benchmark(inverse_map($M), seconds = 3, evals = 1)).time * 1e-9
    t2 = minimum(@benchmark(inverse_fast($M), seconds = 3, evals = 1)).time * 1e-9
    @printf(io, "%d,invmap6,%.4e,%.4e,%.2f,%.2e\n", no, t1, t2, t1 / t2, d3)
    flush(io); println("order $no done")
end
close(io); print(read(joinpath(OUT, "da_fast.csv"), String))