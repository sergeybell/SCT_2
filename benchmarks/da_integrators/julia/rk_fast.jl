# rk_fast.jl - PROTOTYPE of an allocation-light DOP853 step for DA states (SCT_2 only).
#   Bell's _step(::DormandPrince853) builds every stage as z + h*sum(a_ij*k_j) with new DA objects
#   (and an Any-typed stage vector).  Here the stage sums are accumulated in place on the
#   coefficient arrays (one buffer per component), the stage vector is concretely typed.
#   Same coefficients (Bell.Integrators.D8_*), same result to roundoff.
# Run: julia --project=julia julia/rk_fast.jl  -> output/rk_fast.csv
using Bell, Bell.DiffAlgebra, Bell.Integrators, BenchmarkTools, Printf
const BI = Bell.Integrators
function lin_into!(out::Vector{DA{Float64}}, z, h, cs, ks)
    for i in eachindex(z)
        o = out[i].c; copyto!(o, z[i].c)
        for j in eachindex(cs)
            c = cs[j]; iszero(c) && continue
            kc = ks[j][i].c; hc = h * c
            @inbounds @simd for m in eachindex(o); o[m] += hc * kc[m]; end
        end
    end
    out
end
function step853_fast(f, z::Vector{DA{Float64}}, s, h)
    sp = z[1].sp; nz = length(z); newbuf() = [DA(sp, zeros(nmonomials(sp))) for _ in 1:nz]
    ks = Vector{Vector{DA{Float64}}}(undef, 13); ks[1] = f(z, s); buf = newbuf()
    for st in 2:12
        ks[st] = f(lin_into!(buf, z, h, view(BI.D8_A, st, 1:st-1), ks), s + BI.D8_C[st] * h)
    end
    zn = lin_into!(newbuf(), z, h, BI.D8_B, ks); ks[13] = f(zn, s + h)
    e5 = lin_into!(newbuf(), [DA(sp, zeros(nmonomials(sp))) for _ in 1:nz], h, BI.D8_E5, ks)
    e3 = lin_into!(newbuf(), [DA(sp, zeros(nmonomials(sp))) for _ in 1:nz], h, BI.D8_E3, ks)
    zn, (e5, e3)
end
const ETA = 0.2131578052277212; const KQ = 10 / 2.149635688032703
io = open(joinpath(dirname(@__DIR__), "output", "rk_fast.csv"), "w")
println(io, "order,bell_step_s,fast_step_s,speedup,rhs13_s,max_abs_diff")
for no in (3, 5, 7, 10)
    sp = DASpace(6, no); X = da_vars(sp); dK = X[6]
    P2 = (1 + dK) * ((1 + dK) * ETA + 2) / (ETA + 2)
    f(z, s) = (ip = inv(sqrt(P2 - z[2] * z[2] - z[4] * z[4])); [z[2] * ip, -KQ * z[1], z[4] * ip, KQ * z[3]])
    z = X[1:4]; h = 0.02
    a = BI._step(DormandPrince853(), f, z, 0.0, h); b = step853_fast(f, z, 0.0, h)
    d = maximum(coeffnorm(a[1][i] - b[1][i]) for i in 1:4)
    tb = minimum(@benchmark(BI._step(DormandPrince853(), $f, $z, 0.0, $h), seconds = 3, evals = 1)).time * 1e-9
    tf = minimum(@benchmark(step853_fast($f, $z, 0.0, $h), seconds = 3, evals = 1)).time * 1e-9
    tr = 13 * minimum(@benchmark($f($z, 0.0), seconds = 1)).time * 1e-9
    @printf(io, "%d,%.4e,%.4e,%.2f,%.4e,%.2e\n", no, tb, tf, tb / tf, tr, d)
end
close(io); print(read(joinpath(dirname(@__DIR__), "output", "rk_fast.csv"), String))