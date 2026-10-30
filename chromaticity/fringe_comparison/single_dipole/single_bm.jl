# single_bm.jl - Bell maps of a lone BM (drift 0.5 / BM / drift 0.5), coefficients to order 3
#   julia --project=../BELL single_bm.jl
using Bell, Printf
lat = load_structure(joinpath(@__DIR__, "single_bm.bell"))
mono = [e for e in Iterators.product(ntuple(_ -> 0:3, 6)...) if sum(e) in 1:3 && e[5] == 0]
models = [("NoEdge", 0.0, FringeModel()), ("HardEdge2", 0.0, FringeModel(dipole = HardEdge(order = 2))),
          ("HardEdge3", 0.0, FringeModel(dipole = HardEdge(order = 3)))]
for g in (0.05, 0.025, 0.0125)
    push!(models, ("SoftEnge", g, FringeModel(dipole = SoftEdge(profile = Enge(gap = 2g), method = TaylorIntegrate()))))
    push!(models, ("HardEdge2_FINT", g, FringeModel(dipole = HardEdge(order = 2, gap = 2g, fint = fringe_integral(Enge(gap = 2g)) / (2g)))))
end
io = open(joinpath(@__DIR__, "output", "bell_single_bm.csv"), "w")
println(io, "model,gap,out,ex,epx,ey,epy,el,ed,coef")
for (name, g, m) in models
    M = factor_to_map(ring_lie_factor(lat, m))
    for i in 1:6, e in mono
        c = coeff(M[i], collect(e))
        abs(c) > 1e-14 && println(io, name, ",", g, ",", i, ",", join(e, ","), ",", @sprintf("%.15e", c))
    end
end
close(io); println("ok")
