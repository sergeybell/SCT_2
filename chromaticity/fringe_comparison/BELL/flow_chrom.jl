# flow_chrom.jl -- exact high-order chromaticity of magnetic_2_md (RF removed) from the
# Bell flow route: chromaticity_series(one_turn_map_flow(lat, N)).  No finite differences.
#   julia --project=. flow_chrom.jl      -> output/flow_xi.csv
# Lattices:  "file"    = ../../BELL/magnetic_2_md.bell as is (same as run.jl / NoEdge);
#            "matched" = the same with COSY's exact inputs: dipole angle pi/8 (file: 8 digits),
#                        deuteron mass 1875.612945 MeV (COSY PART), k1 = G/Brho with that Brho.
# Fringes are not part of the flow route; the ring has no SBend edge angles (e1 = e2 = 0), no RF.
using Bell, Bell.DiffAlgebra, Printf, LinearAlgebra
const OUT = joinpath(@__DIR__, "output")
lat0 = load_structure(joinpath(@__DIR__, "..", "..", "BELL", "magnetic_2_md.bell"))
latF = Lattice(lat0.name, AbstractElement[e for e in lat0.elements if !(e isa RFCavity)], lat0.ref, lat0.harmonic)
println("elements: ", length(latF.elements), "  SBend: ", count(e -> e isa SBend, latF.elements), "  nonzero edges: ",
        count(e -> e isa SBend && (e.e1 != 0 || e.e2 != 0), latF.elements), "  RF removed: ", length(lat0.elements) - length(latF.elements))
refC = reference_from_kinetic(Particle("deuteron (COSY mass)", 1875.612945, 1.0, lat0.ref.particle.G), 270.0)
fixe(e) = e isa SBend ? SBend(e.name, e.L, pi / 8, e.k1) :
          e isa Quadrupole ? Quadrupole(e.name, e.L, e.k1 * lat0.ref.brho_Tm / refC.brho_Tm) : e
latM = Lattice(lat0.name, AbstractElement[fixe(e) for e in latF.elements], refC, lat0.harmonic)

rows = ["lattice,source,N,plane,q0,xi1,xi2,xi3,xi4,xi5,xi6,xi7,seconds"]
rowx(lt, src, N, pl, q0, xi, t) = push!(rows, string(lt, ",", src, ",", N, ",", pl, ",", @sprintf("%.15g", q0), ",",
    join([k <= length(xi) ? @sprintf("%.15g", xi[k]) : "" for k in 1:7], ","), ",", @sprintf("%.4f", t)))
chromaticity_series(one_turn_map_flow(latF, 3; cache = false))   # compile
for (lt, lat, Ns) in (("file", latF, 3:9), ("matched", latM, 3:8))
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
    # flow map of order 3 with the old finite-difference reader (DA map -> same tunes_at_delta? no: exact series only)
    @printf("%s DF: map %.2f s, exact-series xi1 = %.12f %.12f\n", lt, t, c.xi_x[1], c.xi_y[1])
end
write(joinpath(OUT, "flow_xi.csv"), join(rows, "\n") * "\n")
println("done")