# =============================================================================
# run_v08.jl - nine steps of the thin-RF equilibrium-energy study.
# Ring: magnetic_2 + RFCavity in lat.elements. Models:
#   M1 formula, M2 Lie NF, M3 degree-3 map, M4 exact integrator.
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "steps", "01_longitudinal.jl"))
include(joinpath(@__DIR__, "steps", "02_running_mean.jl"))
include(joinpath(@__DIR__, "steps", "03_add_phase.jl"))
include(joinpath(@__DIR__, "steps", "04_betatron_motion.jl"))
include(joinpath(@__DIR__, "steps", "05_path_lengthening.jl"))
include(joinpath(@__DIR__, "steps", "06_add_betatron.jl"))
include(joinpath(@__DIR__, "steps", "07_lie_and_formula.jl"))
include(joinpath(@__DIR__, "steps", "08_torus.jl"))
include(joinpath(@__DIR__, "steps", "09_truncated_map.jl"))

section("0. Кольцо Bell с тонкой полостью в списке элементов")
S = setup_thin_rf_ring(OPTIM)
ctx = Ctx()
mac!(ctx, "Ekin", "270")
mac!(ctx, "hrf", "5")
mac!(ctx, "Vrf", "0.2")
mac!(ctx, "Circ", @sprintf("%.6f", S.C))
mac!(ctx, "gammas", @sprintf("%.8f", S.ref.gamma))
mac!(ctx, "betas", @sprintf("%.8f", S.ref.beta))
mac!(ctx, "massd", @sprintf("%.8f", S.ref.mass_MeV))
mac!(ctx, "massfile", @sprintf("%.3f", S.of.mass_MeV))
mac!(ctx, "pzero", @sprintf("%.5f", S.ref.p0_MeV))
mac!(ctx, "brho", @sprintf("%.6f", S.ref.brho_Tm))
mac!(ctx, "etazero", @sprintf("%.6f", S.pe.eta0))
mac!(ctx, "etaone", @sprintf("%.6f", S.pe.eta1))
mac!(ctx, "xix", @sprintf("%.4f", S.xi.x[2]))
mac!(ctx, "xiy", @sprintf("%.4f", S.xi.y[2]))
mac!(ctx, "Qx", @sprintf("%.6f", S.nf.mu[1] / 2pi))
mac!(ctx, "Qy", @sprintf("%.6f", S.nf.mu[2] / 2pi))
mac!(ctx, "Qs", @sprintf("%.6f", S.Qs))
mac!(ctx, "betax", @sprintf("%.4f", S.po.x.beta))
mac!(ctx, "betay", @sprintf("%.4f", S.po.y.beta))
mac!(ctx, "alphax", @sprintf("%.4f", S.po.x.alpha))
mac!(ctx, "alphay", @sprintf("%.4f", S.po.y.alpha))
mac!(ctx, "Dx", @sprintf("%.4f", S.po.D4[1]))
mac!(ctx, "Dpx", @sprintf("%.4f", S.po.D4[2]))
nrf = count(e -> e isa RFCavity, S.lat.elements)
@printf("  C=%.6f m  eta0=%.6f eta1=%.6f xi=(%.4f, %.4f) Qs=%.6f\n",
        S.C, S.pe.eta0, S.pe.eta1, S.xi.x[2], S.xi.y[2], S.Qs)
@printf("  elements=%d  RFCavity in list: %d  V=%.3f MV  h=%d\n",
        length(S.lat.elements), nrf, S.V, S.h)
nrf == 1 || error("expected one RFCavity in lat.elements")

step_01!(S, ctx)
step_02!(S, ctx)
step_03!(S, ctx)
step_04!(S, ctx)
step_05!(S, ctx)
step_06!(S, ctx)
step_07!(S, ctx)
step_08!(S, ctx)
step_09!(S, ctx)
write_numbers(ctx)

open(joinpath(OUT, "v08_summary.txt"), "w") do io
    println(io, "thin RF in lat.elements; M1 formula, M2 Lie NF, M3 deg-3 map, M4 exact")
    @printf(io, "C=%.6f eta0=%.6f eta1=%.6f xi=(%.4f,%.4f) Qs=%.6f\n",
            S.C, S.pe.eta0, S.pe.eta1, S.xi.x[2], S.xi.y[2], S.Qs)
    for (k, v) in sort(collect(ctx.macros))
        println(io, k, " = ", v)
    end
end
println("Done.")
