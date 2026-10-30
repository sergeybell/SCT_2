# b2_strength.jl - integrated "sextupole-like" term of Bell's soft edge: K2L = int h0 F'' theta^2 ds per edge,
# and the gradient-like term int h0 F' theta ds (= h0^2 I), for Enge(2*GAP).   julia --project=../BELL b2_strength.jl
using Bell, Printf
h0 = 0.3926990753258494 / 5.329414799999999
open(joinpath(@__DIR__, "output", "b2_strength.csv"), "w") do io
    println(io, "gap,I,K1L_b1,K2L_b2,W,H1_hb1")
    for gap in (0.05, 0.025, 0.0125)
        F, W, _ = profile_function(Enge(gap = 2gap)); n = 20000; s = range(-W, W, length = n + 1); ds = step(s)
        th = cumsum([0.0; [h0 * (F(s[i]) + F(s[i+1])) / 2 * ds for i in 1:n]])
        e = 1e-6
        Fp(z) = (F(z + e) - F(z - e)) / 2e; Fpp(z) = (F(z + e) - 2F(z) + F(z - e)) / e^2
        k1 = sum(h0 * Fp(s[i]) * th[i] for i in 2:n) * ds
        k2 = sum(h0 * Fpp(s[i]) * th[i]^2 for i in 2:n) * ds
        # metric term missing in Bell: (1 + h x) B_y with h = h0 F gives h*b1 in the x^2 force: H1 = int h0 F * h0 F' theta ds
        h1 = sum(h0 * F(s[i]) * h0 * Fp(s[i]) * th[i] for i in 2:n) * ds
        @printf(io, "%g,%.10e,%.10e,%.10e,%.6g,%.10e\n", gap, fringe_integral(Enge(gap = 2gap)), k1, k2, W, h1)
    end
end
println(read(joinpath(@__DIR__, "output", "b2_strength.csv"), String))
