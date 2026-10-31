# flow_chrom.jl -- exact high-order chromaticity of magnetic_2_md (RF removed) from the
# Bell flow route: chromaticity_series(one_turn_map_flow(lat, N)).  No finite differences.
#   julia --project=. flow_chrom.jl      -> output/flow_xi.csv
# Lattice: ../../BELL/magnetic_2_md.bell (= NoEdge of run.jl), BM angle exactly pi/8 (convert_optim close_ring);
# COSY now uses the same deuteron mass 1875.61294257 MeV, so no "matched" copy is needed any more.
# Fringes are not part of the flow route; the ring has no SBend edge angles (e1 = e2 = 0), no RF.
using Bell, Bell.DiffAlgebra, Printf, LinearAlgebra
const OUT = joinpath(@__DIR__, "output")
lat0 = load_structure(joinpath(@__DIR__, "..", "..", "BELL", "magnetic_2_md.bell"))
latF = Lattice(lat0.name, AbstractElement[e for e in lat0.elements if !(e isa RFCavity)], lat0.ref, lat0.harmonic)
println("elements: ", length(latF.elements), "  SBend: ", count(e -> e isa SBend, latF.elements), "  nonzero edges: ",
        count(e -> e isa SBend && (e.e1 != 0 || e.e2 != 0), latF.elements), "  RF removed: ", length(lat0.elements) - length(latF.elements))
println("BM angle - pi/8 = ", first(e for e in latF.elements if e isa SBend).angle - pi / 8, ",  mass = ", latF.ref.mass_MeV)

rows = ["lattice,source,N,plane,q0,xi1,xi2,xi3,xi4,xi5,xi6,xi7,seconds"]
rowx(lt, src, N, pl, q0, xi, t) = push!(rows, string(lt, ",", src, ",", N, ",", pl, ",", @sprintf("%.15g", q0), ",",
    join([k <= length(xi) ? @sprintf("%.15g", xi[k]) : "" for k in 1:7], ","), ",", @sprintf("%.4f", t)))
chromaticity_series(one_turn_map_flow(latF, 3; cache = false))   # compile
for (lt, lat, Ns) in (("file", latF, 3:9),)
    for N in Ns
        tm = @elapsed (M = one_turn_map_flow(lat, N; cache = false))
        tc = @elapsed (c = chromaticity_series(M))
        @printf("%s N=%d  map %.2f s  chrom %.2f s  symp %.1e\n", lt, N, tm, tc, da_symplectic_error(M))
        @printf("   x: Q0=%.13f xi=%s\n   y: Q0=%.13f xi=%s\n", c.x[1], c.xi_x, c.y[1], c.xi_y); flush(stdout)
        rowx(lt, "flow", N, "x", c.x[1], c.xi_x, tm + tc); rowx(lt, "flow", N, "y", c.y[1], c.xi_y, tm + tc)
        N == 7 && println("   coupling |R[1:2,3:4]| = ", norm(linear_part(M)[1:2, 3:4]))
    end
    # old Dragt-Finn / Lie route (order-3 map): finite differences (run.jl uses dd = 1e-4) and the exact series of that map
    t = @elapsed (MDF = factor_to_map(ring_lie_factor(lat)))
    for dd in (1e-5, 1e-4, 1e-3)
        h = chromaticity_orders(MDF; dd = dd)
        rowx(lt, "DF_findiff_dd$(dd)", 3, "x", h.x[1], [h.x[2], h.x[3]], t); rowx(lt, "DF_findiff_dd$(dd)", 3, "y", h.y[1], [h.y[2], h.y[3]], t)
    end
    sp = flow_space(3)
    tc(p, ex) = Bell.LieAlgebra.TaylorSeries.getcoeff(p, ex)
    MD = [DA(sp, [tc(MDF[i], monomial_exponents(sp, m)) for m in 1:nmonomials(sp)]) for i in 1:6]
    c = chromaticity_series(MD)
    rowx(lt, "DF_exact_series", 3, "x", c.x[1], c.xi_x, t); rowx(lt, "DF_exact_series", 3, "y", c.y[1], c.xi_y, t)
    @printf("%s DF: map %.2f s, exact-series xi1 = %.12f %.12f\n", lt, t, c.xi_x[1], c.xi_y[1])
end
write(joinpath(OUT, "flow_xi.csv"), join(rows, "\n") * "\n")
println("done")