# ring_bm_optics.jl - for every BM of magnetic_2_md: one-turn linear matrix at its entrance and the
# dispersion vector D1 (dz/ddelta of the closed orbit) there; ring with Bell SoftEdge Enge(2*GAP) on dipoles.
# Also Bell maps of the bare BM.   julia --project=../BELL ring_bm_optics.jl
using Bell, Printf
src = read(joinpath(@__DIR__, "..", "..", "BELL", "magnetic_2_md.bell"), String)
head, seq = split(src, r"\nsequence\s*\n")
toks = split(seq)
idx = findall(==("BM"), toks)
io = open(joinpath(@__DIR__, "output", "ring_bm_optics.csv"), "w")
println(io, "gap,k,pos,", join(["M$(i)$(j)" for i in 1:4 for j in 1:4], ","), ",D1,D2,D3,D4,Qx,Qy")
for gap in (0.05, 0.025, 0.0125)
    model = FringeModel(dipole = SoftEdge(profile = Enge(gap = 2gap), method = TaylorIntegrate()))
    for (k, p) in enumerate(idx)
        rot = vcat(toks[p:end], toks[1:p-1])
        f = tempname() * ".bell"; write(f, head * "\nsequence\n  " * join(rot, " ") * "\n")
        lat = load_structure(f)
        M = factor_to_map(ring_lie_factor(lat, model))
        d = dispersion_orders(M); R = d.R
        ch = chromaticity_orders(M)
        println(io, gap, ",", k, ",", p, ",", join([@sprintf("%.15e", R[i, j]) for i in 1:4 for j in 1:4], ","), ",",
                join([@sprintf("%.15e", v) for v in d.D1], ","), ",", ch.x[1], ",", ch.y[1])
    end
end
close(io)
lat = load_structure(joinpath(@__DIR__, "single_bm_bare.bell"))
mono = [e for e in Iterators.product(ntuple(_ -> 0:3, 6)...) if sum(e) in 1:3 && e[5] == 0]
io = open(joinpath(@__DIR__, "output", "bell_bare_bm.csv"), "w")
println(io, "model,gap,out,ex,epx,ey,epy,el,ed,coef")
for (name, g, m) in [("NoEdge", 0.0, FringeModel()), ("HardEdge3", 0.0, FringeModel(dipole = HardEdge(order = 3))),
                     [("SoftEnge", g, FringeModel(dipole = SoftEdge(profile = Enge(gap = 2g), method = TaylorIntegrate()))) for g in (0.05, 0.025, 0.0125)]...]
    M = factor_to_map(ring_lie_factor(lat, m))
    for i in 1:6, e in mono
        c = coeff(M[i], collect(e)); abs(c) > 1e-14 && println(io, name, ",", g, ",", i, ",", join(e, ","), ",", @sprintf("%.15e", c))
    end
end
close(io); println("ok")
