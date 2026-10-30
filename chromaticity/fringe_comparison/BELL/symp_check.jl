using Bell, Printf
lat = load_structure(joinpath(@__DIR__, "..", "..", "BELL", "magnetic_2_md.bell"))
for g in (0.05, 0.00625), (mn, m) in [("Splitting", Splitting(n = 40, order = 4)), ("Magnus", Magnus(n = 40, order = 4)), ("TaylorIntegrate", TaylorIntegrate())]
    o = fringe_observables(lat, FringeModel(dipole = SoftEdge(profile = Enge(gap = 2g), method = m)))
    @printf("gap=%g %-16s symp_edge_lin=%.2e symp_edge_pt=%.2e symp_ring=%.2e\n", g, mn, o.symp_edge_lin, o.symp_edge_pt, o.symp_ring)
end
