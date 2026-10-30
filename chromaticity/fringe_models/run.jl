# fringe_models/run.jl -- Bell fringe-field models vs COSY / MAD-X / OptiM on magnetic_2_md
#
#   run_fringe.bat   (= julia --project=. run.jl)
#
# Bell is used as a library (Project.toml: Bell dev'd from C:\GIT_REPOS\BELL).
# Reference numbers are copied from the other folders of SCT_2\chromaticity:
#   COSY  ..\COSY\output\cosy_ov321_fr.txt (OV 3 2 1 + TP, RINGFR FD FQ GAP)
#   MAD-X ..\analysis\exact_hamiltonian_check.out / ..\MADX (hard edge, E1=E2=0, FINT=0)
#   OptiM xi_x = -1.4615
# COSY GAP is the half-aperture: its Enge length scale is the full gap 2*GAP,
# which is Bell's Enge(gap = 2*GAP).
using Bell, Printf, Plots
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

const LAT = joinpath(@__DIR__, "..", "BELL", "magnetic_2_md.bell")
const OUT = joinpath(@__DIR__, "output"); mkpath(OUT)
lat = load_structure(LAT)

# ----------------------------------------------------------------- references
const COSY_FR0 = (Qx = 0.4617201483757378, Qy = 0.4613373212771295, xix = -1.461535714233974, xiy = -2.225484930000162)
# FD = 3 (dipole Enge fringe), FQ = 0; GAP [m] => (Qx, Qy, xi_x, xi_y)
const COSY_FR3 = [(0.05,    0.4617204423777027, 0.4557293388816404, -1.461498245468243, -1.807670104100061),
                  (0.025,   0.4617201851255041, 0.4585248206842834, -1.461525973542352, -1.812335929222331),
                  (0.0125,  0.4617201529694136, 0.4599289112349980, -1.461533163365978, -1.814691819197023),
                  (0.00625, 0.4617201489499543, 0.4606325722640402, -1.461535027478227, -1.815875595028895)]
const COSY_FR3_LIMIT = -1.817037      # extrapolation GAP -> 0 (COSY README)
const MADX_HARD = (xix = -1.46153562, xiy = -1.81706334)
const OPTIM_XIX = -1.4615

obs(m) = fringe_observables(lat, m)
io = open(joinpath(OUT, "fringe_vs_codes.txt"), "w")
out(args...) = (print(stdout, args...); print(io, args...))
outf(fmt, args...) = out(Printf.format(Printf.Format(fmt), args...))

out("Bell fringe models vs COSY / MAD-X / OptiM, lattice magnetic_2_md\n\n")

# ----------------------------------------------------------------- 1. hard edge
o0 = obs(FringeModel())
oh1 = obs(FringeModel(dipole = HardEdge(order = 1)))
oh2 = obs(FringeModel(dipole = HardEdge(order = 2)))
out("1. No fringe / hard edge\n")
outf("  %-34s %10s %10s %11s %11s\n", "model", "Qx", "Qy", "xi_x", "xi_y")
outf("  %-34s %10.6f %10.6f %11.6f %11.6f\n", "Bell NoEdge", o0.Qx, o0.Qy, o0.xix, o0.xiy)
outf("  %-34s %10.6f %10.6f %11.6f %11.6f\n", "COSY FR 0", COSY_FR0.Qx, COSY_FR0.Qy, COSY_FR0.xix, COSY_FR0.xiy)
outf("  %-34s %10.6f %10.6f %11.6f %11.6f\n", "Bell HardEdge(order=1)", oh1.Qx, oh1.Qy, oh1.xix, oh1.xiy)
outf("  %-34s %10.6f %10.6f %11.6f %11.6f\n", "Bell HardEdge(order=2)", oh2.Qx, oh2.Qy, oh2.xix, oh2.xiy)
outf("  %-34s %10s %10s %11.6f %11.6f\n", "MAD-X hard edge (E=0, FINT=0)", "", "", MADX_HARD.xix, MADX_HARD.xiy)
outf("  %-34s %10s %10s %11.4f %11s\n", "OptiM", "", "", OPTIM_XIX, "")
outf("  Bell NoEdge   - COSY FR0 : dxi_x = %.2e  dxi_y = %.2e\n", o0.xix - COSY_FR0.xix, o0.xiy - COSY_FR0.xiy)
outf("  Bell Hard(2)  - MAD-X    : dxi_x = %.2e  dxi_y = %.2e  (rel %.2e)\n",
     oh2.xix - MADX_HARD.xix, oh2.xiy - MADX_HARD.xiy, (oh2.xiy - MADX_HARD.xiy) / MADX_HARD.xiy)
outf("  Bell Hard(2)  - COSY FR3 limit (%.6f): dxi_y = %.2e\n\n", COSY_FR3_LIMIT, oh2.xiy - COSY_FR3_LIMIT)

# ----------------------------------------------------------------- 2. soft edge vs COSY FR 3
out("2. Enge soft edge (COSY default coefficients) vs COSY FR 3, dipoles only\n")
out("   Bell: SoftEdge(Enge(gap = 2*GAP), Splitting(n=40, order=4)) and Magnus(n=40, order=4)\n")
outf("  %8s | %10s %10s | %11s %11s %9s | %11s %9s\n", "GAP [m]", "Qy Bell", "Qy COSY", "xi_y Bell", "xi_y COSY", "diff", "xi_y Magnus", "diff")
bq, bx = Float64[], Float64[]
for (gap, cqx, cqy, cxx, cxy) in COSY_FR3
    os = obs(FringeModel(dipole = SoftEdge(profile = Enge(gap = 2gap), method = Splitting(n = 40, order = 4))))
    om = obs(FringeModel(dipole = SoftEdge(profile = Enge(gap = 2gap), method = Magnus(n = 40, order = 4))))
    push!(bq, os.Qy); push!(bx, os.xiy)
    outf("  %8.5f | %10.6f %10.6f | %11.6f %11.6f %9.2e | %11.6f %9.2e\n", gap, os.Qy, cqy, os.xiy, cxy, os.xiy - cxy, om.xiy, om.xiy - cxy)
    outf("  %8s | Qx %10.6f %10.6f | xi_x %9.6f %11.6f %9.2e\n", "", os.Qx, cqx, os.xix, cxx, os.xix - cxx)
end
# linear extrapolation in GAP of the two smallest gaps
g = [c[1] for c in COSY_FR3]
lin(y) = y[end] - (y[end-1] - y[end]) / (g[end-1] - g[end]) * g[end]
outf("\n  GAP -> 0 (linear from two smallest gaps): Bell %.6f  COSY %.6f  (COSY README %.6f); Bell HardEdge(2) %.6f\n\n",
     lin(bx), lin([c[5] for c in COSY_FR3]), COSY_FR3_LIMIT, oh2.xiy)

# ----------------------------------------------------------------- 3. linear FINT
I = fringe_integral(Enge(gap = 0.1))
ohf = obs(FringeModel(dipole = HardEdge(order = 2, gap = 0.1, fint = I / 0.1)))
out("3. Hard edge with FINT of the Enge profile (gap 0.1 m = COSY GAP 0.05)\n")
outf("  FINT = %.5f   Bell Hard(2)+FINT: Qy = %.6f xi_y = %.6f   COSY FR3: Qy = %.6f xi_y = %.6f\n\n",
     I / 0.1, ohf.Qy, ohf.xiy, COSY_FR3[1][3], COSY_FR3[1][5])
close(io)

p = scatter([c[1] for c in COSY_FR3], [c[5] for c in COSY_FR3], label = "COSY FR 3", marker = :square,
            xlabel = "COSY GAP (half aperture) [m]", ylabel = "xi_y", title = "magnetic_2_md: vertical chromaticity with dipole fringe")
plot!(p, [c[1] for c in COSY_FR3], bx, marker = :o, label = "Bell SoftEdge(Enge(2 GAP)), Splitting(40,4)")
hline!(p, [oh2.xiy], ls = :dash, label = "Bell HardEdge(order=2)")
hline!(p, [MADX_HARD.xiy], ls = :dot, label = "MAD-X hard edge")
hline!(p, [COSY_FR3_LIMIT], ls = :dashdot, label = "COSY FR 3, GAP -> 0")
savefig(p, joinpath(OUT, "xi_y_vs_gap.png"))
println("written: ", joinpath(OUT, "fringe_vs_codes.txt"), ", xi_y_vs_gap.png")
