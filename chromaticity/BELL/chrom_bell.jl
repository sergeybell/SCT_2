# chrom_bell.jl -- tunes, alpha_c, eta, chromaticity of magnetic_2_md via BELL (as-is)
#   julia --project=C:\GIT_REPOS\BELL chrom_bell.jl
using Bell, Printf, LinearAlgebra, TaylorSeries

const OPT  = raw"C:\GIT_REPOS\SCT_2\lattice\magnetic_2_md.opt"
const STR  = joinpath(@__DIR__, "magnetic_2_md.bell")
const OUT  = joinpath(@__DIR__, "output", "bell_results.txt")
mkpath(dirname(OUT))

of  = convert_optim(OPT, STR)          # BELL OptiM reader -> .bell structure
lat = load_structure(STR)
ref = lat.ref
C   = circumference(lat)
foreach(w -> println("WARNING: ", w), of.warnings)

# (A) linear optics: one-turn matrix
R  = ring_matrix(lat)
po = periodic_optics(R, ref, C)
opt = optics_along(lat, po)

# (B) one-turn Taylor map (RK4 in truncated polynomials), chromaticity_orders
t0 = time(); M = ring_taylor_map(lat); tM = time() - t0
pe = path_expansion(M, ref, C)
ch = chromaticity_orders(M)

# (C) one-turn Lie factor (Lie map) -> polynomial -> same readers
obs = observables(lat)

io = open(OUT, "w")
for s in (stdout, io)
  @printf(s, "# BELL (Julia) results for %s\n", OPT)
  @printf(s, "# reference: mass=%.8f MeV Ekin=%.4f MeV gamma=%.9f beta=%.9f p=%.6f MeV/c Brho=%.8f T*m\n",
          ref.mass_MeV, ref.Ekin_MeV, ref.gamma, ref.beta, ref.p0_MeV, ref.brho_Tm)
  @printf(s, "# delta = dp/p (canonical momentum deviation); xi = dQ/ddelta (central difference of Jacobian tunes on off-momentum closed orbit)\n")
  @printf(s, "C_m = %.7f\n", C)
  @printf(s, "linear_Qx = %.8f\nlinear_Qy = %.8f\n", opt.Qx, opt.Qy)
  @printf(s, "linear_alpha_c = %.8f\nlinear_eta = %.8f\n", po.alpha_c, po.eta_slip)
  @printf(s, "taylor_alpha0 = %.8f\ntaylor_alpha1 = %.8f\ntaylor_eta0 = %.8f\ntaylor_eta1 = %.8f\n", pe.alpha0, pe.alpha1, pe.eta0, pe.eta1)
  @printf(s, "taylor_Qx_frac = %.8f\ntaylor_xi_x = %.6f\ntaylor_xi2_x = %.4f\n", ch.x...)
  @printf(s, "taylor_Qy_frac = %.8f\ntaylor_xi_y = %.6f\ntaylor_xi2_y = %.4f\n", ch.y...)
  @printf(s, "lie_Qx = %.8f\nlie_Qy = %.8f\nlie_alpha0 = %.8f\nlie_eta0 = %.8f\n", obs.po.x.mu/2pi, obs.po.y.mu/2pi, obs.pe.alpha0, obs.pe.eta0)
  @printf(s, "lie_xi_x = %.6f\nlie_xi_y = %.6f\nlie_xi2_x = %.4f\nlie_xi2_y = %.4f\n", obs.xi.x[2], obs.xi.y[2], obs.xi.x[3], obs.xi.y[3])
  @printf(s, "# taylor map built in %.1f s\n", tM)
end
close(io)
open(joinpath(@__DIR__, "output", "bell_twiss.csv"), "w") do f
  println(f, "s_m,betx_m,bety_m,Dx_m")           # D per dp/p, at element boundaries (optics_along)
  for k in eachindex(opt.data.s)
    @printf(f, "%.9f,%.9f,%.9f,%.9e\n", opt.data.s[k], opt.data.betx[k], opt.data.bety[k], opt.data.D[k])
  end
end
println("written ", OUT)
