# =============================================================================
# run_v08.jl - delta_eq with THIN RF: four models M1..M4.
#              One-turn object is the Lie product; NF on generators;
#              torus average of delta.
# Outputs: output/v08_*.png, output/v08_summary.txt, note/figs/*.png,
#          note/gen/*.tex (numbers and tables \input by the note).
# =============================================================================
using EquilibriumEnergy
using Bell, Printf, Plots, Statistics, LinearAlgebra
gr(); default(fontfamily = "sans-serif", framestyle = :box, grid = true, dpi = 150)

const OPTIM = raw"C:\GIT_REPOS\SCT_2\OptiM\magnetic\magnetic_2.opt"
const OUT = joinpath(@__DIR__, "..", "output")
const NOTE = joinpath(@__DIR__, "..", "note")
const GEN = joinpath(NOTE, "gen"); const FIG = joinpath(NOTE, "figs")
mkpath(OUT); mkpath(GEN); mkpath(FIG)
section(t) = println("\n", "="^72, "\n ", t, "\n", "="^72)
savefig2(p, name) = (savefig(p, joinpath(OUT, "v08_" * name * ".png")); savefig(p, joinpath(FIG, name * ".png")))
tex(x; d = 4) = (x == 0 ? "0" : (e = floor(Int, log10(abs(x))); m = x / 10.0^e;
                 e == 0 ? @sprintf("%.*f", d, m) : @sprintf("%.*f\\cdot10^{%d}", d, m, e)))
macros = Dict{String,String}()
mac(k, v) = (macros[k] = v)

# ---------------------------------------------------------------- setup
section("0. Setup (thin RF, magnetic lattice)")
const EKIN = 270.0; const HARM = 5; const VRF = 0.2; const PHIS = Float64(pi)
lat, ref, of = load_magnetic(OPTIM; Ekin = EKIN, h = HARM)
C = circumference(lat)
Fmag = ring_lie_factor(lat)
Mmag = factor_to_map(Fmag)
pe = path_expansion(Mmag, ref, C)
xi = chromaticity_orders(Mmag)
po = periodic_optics(ring_matrix(lat), ref, C)
cav = RFCavity("RF", VRF, HARM, PHIS)
Qs_smooth = synchrotron_tune(ref, pe.eta0, VRF, HARM, PHIS; T = 1.0)
F = one_turn_factor_rf(lat, cav, ref, C)
tnf = @elapsed nf = dragt_finn_nf6d(F)
M = factor_to_map(F)
Qlin = nf.mu ./ (2pi)
Qs = abs(Qlin[3])
@printf("  C=%.6f m  eta0=%.6f eta1=%.6f xi=(%.4f, %.4f) Qs(lin)=%.6f Qs(smooth)=%.6f\n",
        C, pe.eta0, pe.eta1, xi.x[2], xi.y[2], Qs, Qs_smooth)
@printf("  mass used %.8f MeV (OptiM file %.3f) gamma=%.9f beta=%.9f\n", ref.mass_MeV, of.mass_MeV, ref.gamma, ref.beta)
@printf("  NF: z*=%s\n", string(nf.zstar))
@printf("  NF: Q=(%.6f, %.6f, %.6f) |A^-1RA-rot|=%.1e |A'SA-S|=%.1e f3 %.2e -> %.1e, f4nr %.1e (%.1f s)\n",
        Qlin..., nf.Rn_residual, nf.symp_residual, nf.f3_before, nf.f3_after, nf.f4_nonres_after, tnf)
d_note = [2pi / C * xi.x[2] / pe.eta0, 2pi / C * xi.y[2] / pe.eta0, -pe.eta1 * nf.sigma_delta_s^2 / pe.eta0]
@printf("  d (NF)   = %.8e %.8e %.8e\n  d (grid) = %.8e %.8e %.8e\n  d (M1)   = %.8e %.8e %.8e\n",
        nf.d..., nf.d_grid..., d_note...)
println("  h4 c = ", nf.c)

mac("Ekin", "270"); mac("hrf", "5"); mac("Vrf", "0.2")
mac("Circ", @sprintf("%.6f", C)); mac("gammas", @sprintf("%.8f", ref.gamma)); mac("betas", @sprintf("%.8f", ref.beta))
mac("massd", @sprintf("%.8f", ref.mass_MeV)); mac("massfile", @sprintf("%.3f", of.mass_MeV))
mac("pzero", @sprintf("%.5f", ref.p0_MeV)); mac("brho", @sprintf("%.6f", ref.brho_Tm))
mac("etazero", @sprintf("%.6f", pe.eta0)); mac("etaone", @sprintf("%.6f", pe.eta1))
mac("xix", @sprintf("%.4f", xi.x[2])); mac("xiy", @sprintf("%.4f", xi.y[2]))
mac("Qx", @sprintf("%.6f", Qlin[1])); mac("Qy", @sprintf("%.6f", Qlin[2])); mac("Qs", @sprintf("%.6f", Qs))
mac("Qssmooth", @sprintf("%.6f", Qs_smooth))
mac("betax", @sprintf("%.4f", po.x.beta)); mac("betay", @sprintf("%.4f", po.y.beta))
mac("alphax", @sprintf("%.4f", po.x.alpha)); mac("alphay", @sprintf("%.4f", po.y.alpha))
mac("Dx", @sprintf("%.4f", po.D4[1])); mac("Dpx", @sprintf("%.4f", po.D4[2]))
mac("zstarsix", tex(nf.zstar[6]; d = 2)); mac("zstarfive", tex(nf.zstar[5]; d = 2))
mac("RnRes", tex(nf.Rn_residual; d = 1)); mac("SympRes", tex(nf.symp_residual; d = 1))
mac("fthreeBefore", tex(nf.f3_before; d = 2)); mac("fthreeAfter", tex(nf.f3_after; d = 1))
mac("ffourNR", tex(nf.f4_nonres_after; d = 1))
for (k, n) in enumerate(("x", "y", "s"))
    mac("dNF" * n, tex(nf.d[k]; d = 6)); mac("dNote" * n, tex(d_note[k]; d = 6))
    mac("dRel" * n, tex(nf.d[k] / d_note[k] - 1; d = 1))
end
mac("sigds", @sprintf("%.6f", nf.sigma_delta_s)); mac("sigls", @sprintf("%.4f", nf.sigma_l_s))
for (i, a) in enumerate(("x", "y", "s")), (j, b) in enumerate(("x", "y", "s"))
    j >= i && mac("c" * a * b, tex(nf.c[i, j]; d = 4))
end

# ---------------------------------------------------------------- models
step_M3(z; ds = 0.1) = rf_kick(exact_track_turn(z, lat, ref, C; ds = ds), cav, ref, C)
step_M4(z) = evaluate_map(M, z)
function run(step, z0, N)
    Zs = Matrix{Float64}(undef, 6, N); z = copy(z0)
    for n in 1:N
        z = step(z); Zs[:, n] = z
    end
    return Zs
end
deltas(Zs) = Zs[6, :]
phis_note(Zs) = [note_phi(canonical_to_user(Zs[:, n], ref, C, HARM)[5]) for n in 1:size(Zs, 2)]
csJ(A, p) = (1 + p.alpha^2) * A^2 / (2 * p.beta)          # CS action of (x=A, x'=0)
M1(dm, Jx, Jy; phat = 0.0) = -(pe.eta1 * dm^2 / 2 * (1 - phat^2 / 32) +
                              formula_B_Lambda(C, Jx, xi.x[2], Jy, xi.y[2])) / pe.eta0
init(dm, Ax, Ay) = user_to_canonical([po.D4[1] * dm + Ax, po.D4[2] * dm, Ay, 0.0, 0.0, dm], ref, C, HARM)
function period_mean(d, Q, nd)
    K = floor(Int, (length(d) - nd) * Q); n = round(Int, K / Q)
    return mean(d[nd+1:nd+n])
end
function analyse(Zs; Ax = 0.0, Ay = 0.0, z0)
    d = deltas(Zs); p = phis_note(Zs)
    m = measure_synchrotron(d, p; Qs = Qs)
    nd = length(d) - m.n_used
    J0 = nf6d_actions(nf, z0); J0lin = nf6d_actions(nf, z0; linear = true)
    Jx = csJ(Ax, po.x); Jy = csJ(Ay, po.y)
    return (m = m, pm = period_mean(d, m.Qs_meas, nd), J0 = J0, J0lin = J0lin, Jx = Jx, Jy = Jy,
            m1 = M1(m.delta_m_fit, Jx, Jy), m1p = M1(m.delta_m_fit, Jx, Jy; phat = m.phi_hat),
            m2 = nf6d_delta_eq(nf, J0), dm2 = sqrt(2 * J0[3]) * nf.sigma_delta_s)
end

# ---------------------------------------------------------------- Test A
section("1. Test A: delta_eq vs delta_m (no betatron)")
NA = 12000
dmsA = [0.5e-3, 0.75e-3, 1e-3, 1.5e-3, 2e-3]
rowsA = []; seriesA = Dict{Float64,Vector{Float64}}()
for dm in dmsA
    z0 = init(dm, 0.0, 0.0)
    t3 = @elapsed Z3 = run(step_M3, z0, NA)
    Z4 = run(step_M4, z0, NA)
    r = analyse(Z3; z0 = z0)
    m4 = measure_synchrotron(deltas(Z4), phis_note(Z4); Qs = Qs)
    seriesA[dm] = deltas(Z3)
    Qnf = nf6d_tunes(nf, r.J0).Q[3]; Qnf = min(Qnf, 1 - Qnf)
    push!(rowsA, (dm0 = dm, r = r, m4 = m4.delta_eq, m4err = m4.err_split, Qnf = Qnf,
                  Qtr = r.m.Qs_meas, Qnote = Qs * (1 - r.m.phi_hat^2 / 16)))
    @printf("  dm0=%.2e dm_fit=%.5e phat=%.4f | M1=%+.5e M1p=%+.5e M2=%+.5e M3=%+.5e(+-%.1e, pm %+.5e) M4=%+.5e | Qs NF=%.6f tr=%.6f note=%.6f (%.0fs)\n",
            dm, r.m.delta_m_fit, r.m.phi_hat, r.m1, r.m1p, r.m2, r.m.delta_eq, r.m.err_split, r.pm, m4.delta_eq,
            Qnf, r.m.Qs_meas, Qs * (1 - r.m.phi_hat^2 / 16), t3)
end
# integration-step check
zc = init(2e-3, 0.0, 0.0)
Zc = run(z -> step_M3(z; ds = 0.05), zc, 4000)
mc = measure_synchrotron(deltas(Zc), phis_note(Zc); Qs = Qs)
Zc1 = run(step_M3, zc, 4000)
mc1 = measure_synchrotron(deltas(Zc1), phis_note(Zc1); Qs = Qs)
@printf("  ds check (dm=2e-3, 4000 turns): ds=0.1 %.6e  ds=0.05 %.6e  diff %.1e\n", mc1.delta_eq, mc.delta_eq, mc.delta_eq - mc1.delta_eq)
mac("dsCheck", tex(abs(mc.delta_eq - mc1.delta_eq); d = 1)); mac("NA", string(NA))

e6(x) = @sprintf("%.4f", x * 1e6)
open(joinpath(GEN, "tabA.tex"), "w") do io
    for row in rowsA
        r = row.r
        @printf(io, "%.3f & %.5f & %.4f & %s & %s & %s & %s & %.1e & %s & %+.2f & %+.2f & %+.2f\\\\\n",
                row.dm0 * 1e3, r.m.delta_m_fit * 1e3, r.m.phi_hat, e6(r.m1), e6(r.m1p), e6(r.m2), e6(r.m.delta_eq),
                r.m.err_split * 1e6, e6(row.m4),
                100 * (r.m1 / r.m.delta_eq - 1), 100 * (r.m2 / r.m.delta_eq - 1), 100 * (row.m4 / r.m.delta_eq - 1))
    end
end
open(joinpath(GEN, "tabQs.tex"), "w") do io
    for row in rowsA
        r = row.r
        @printf(io, "%.3f & %.4f & %.6f & %.6f & %.6f & %.6f & %.1e & %.1e\\\\\n", row.dm0 * 1e3, r.m.phi_hat,
                r.J0[3] * 1e6, row.Qtr, row.Qnf, row.Qnote, row.Qnf - row.Qtr, row.Qnote - row.Qtr)
    end
end
open(joinpath(GEN, "tabAmean.tex"), "w") do io
    for row in rowsA
        r = row.r
        @printf(io, "%.3f & %s & %s & %s & %.2e & %.5f & %.5f\\\\\n", row.dm0 * 1e3, e6(r.m.delta_eq), e6(r.pm),
                e6(r.m.delta_eq_h1 - r.m.delta_eq_h2), r.m.rms, r.dm2 * 1e3, r.m.delta_m_fit * 1e3)
    end
end
xA = [r.r.m.delta_m_fit^2 * 1e6 for r in rowsA]
pA = plot(xlabel = "δ_m² · 10⁶", ylabel = "δ_eq · 10⁶", legend = :topleft, size = (720, 480))
plot!(pA, [0; xA], [0; [M1(sqrt(x * 1e-6), 0, 0) for x in xA]] .* 1e6, lw = 2, label = "M1: −η₁δ_m²/(2η₀)")
plot!(pA, xA, [r.r.m2 for r in rowsA] .* 1e6, lw = 2, ls = :dash, label = "M2: нормальная форма")
scatter!(pA, xA, [r.r.m.delta_eq for r in rowsA] .* 1e6, ms = 7, label = "M3: трекинг")
scatter!(pA, xA, [r.m4 for r in rowsA] .* 1e6, ms = 4, shape = :diamond, label = "M4: усечённая карта")
savefig2(pA, "fig_A_delta_eq")

# ---------------------------------------------------------------- Test B
section("2. Test B: betatron contribution (Lambda)")
NB = 6000; dmB = 3e-4
ampsB = [5e-3, 10e-3, 15e-3, 20e-3]
rowsB = []
for pl in (:X, :Y), A in ampsB
    Ax = pl == :X ? A : 0.0; Ay = pl == :Y ? A : 0.0
    z0 = init(dmB, Ax, Ay)
    Z3 = run(step_M3, z0, NB)
    r = analyse(Z3; Ax = Ax, Ay = Ay, z0 = z0)
    Lnote = formula_B_Lambda(C, r.Jx, xi.x[2], r.Jy, xi.y[2])
    Lfit = -pe.eta0 * r.m.delta_eq - pe.eta1 * r.m.delta_m_fit^2 / 2
    LNF = -pe.eta0 * (nf.d[1] * r.J0[1] + nf.d[2] * r.J0[2])
    row = pl == :X ? 1 : 3
    Qb = measure_tune(Z3[row, :]); Qb = min(Qb, 1 - Qb)
    Qbnf = nf6d_tunes(nf, r.J0).Q[pl == :X ? 1 : 2]; Qbnf = min(Qbnf, 1 - Qbnf)
    push!(rowsB, (pl = pl, A = A, r = r, Lnote = Lnote, Lfit = Lfit, LNF = LNF, Qb = Qb, Qbnf = Qbnf))
    @printf("  %s A=%2.0f mm J_CS=%.4e J_NF=%.4e | Lam note=%+.4e NF=%+.4e fit=%+.4e | deq M1=%+.5e M2=%+.5e M3=%+.5e(+-%.1e) | Q tr=%.6f NF=%.6f\n",
            pl, A * 1e3, pl == :X ? r.Jx : r.Jy, r.J0[pl == :X ? 1 : 2], Lnote, LNF, Lfit, r.m1, r.m2, r.m.delta_eq, r.m.err_split, Qb, Qbnf)
end
mac("NB", string(NB)); mac("dmB", "3\\cdot10^{-4}")
open(joinpath(GEN, "tabB.tex"), "w") do io
    for b in rowsB
        r = b.r; J = b.pl == :X ? r.Jx : r.Jy; Jn = r.J0[b.pl == :X ? 1 : 2]
        @printf(io, "%s & %.0f & %.4f & %.4f & %s & %s & %s & %s & %s & %s & %+.2f & %+.2f\\\\\n",
                b.pl == :X ? "x" : "y", b.A * 1e3, J * 1e6, Jn * 1e6, e6(b.Lnote), e6(b.LNF), e6(b.Lfit),
                e6(r.m1), e6(r.m2), e6(r.m.delta_eq), 100 * (r.m1 / r.m.delta_eq - 1), 100 * (r.m2 / r.m.delta_eq - 1))
    end
end
open(joinpath(GEN, "tabBQ.tex"), "w") do io
    for b in rowsB
        Jn = b.r.J0[b.pl == :X ? 1 : 2]
        @printf(io, "%s & %.0f & %.4f & %.6f & %.6f & %.1e\\\\\n", b.pl == :X ? "x" : "y", b.A * 1e3, Jn * 1e6,
                b.Qb, b.Qbnf, b.Qbnf - b.Qb)
    end
end
pB = plot(xlabel = "J, мкм (= 10⁻⁶ м)", ylabel = "Λ · 10⁶", legend = :topleft, size = (720, 480))
for (pl, col) in ((:X, 1), (:Y, 2))
    rs = filter(b -> b.pl == pl, rowsB); nm = pl == :X ? "x" : "y"
    Js = [pl == :X ? b.r.Jx : b.r.Jy for b in rs]
    plot!(pB, [0; Js] .* 1e6, [0; [formula_B_Lambda(C, pl == :X ? J : 0, xi.x[2], pl == :Y ? J : 0, xi.y[2]) for J in Js]] .* 1e6,
          lw = 2, color = col, label = "M1: −(2π/C)J_$(nm)ξ_$(nm)")
    scatter!(pB, Js .* 1e6, [b.Lfit for b in rs] .* 1e6, color = col, ms = 6, label = "M3: −η₀δ_eq−η₁δ_m²/2 ($nm)")
    scatter!(pB, [b.r.J0[pl == :X ? 1 : 2] for b in rs] .* 1e6, [b.LNF for b in rs] .* 1e6, color = col, shape = :x, ms = 7, label = "M2: −η₀ d_$(nm) J_$(nm)")
end
savefig2(pB, "fig_B_lambda")

# ---------------------------------------------------------------- Fast betatron averaging (magnetic only)
section("3. Fast betatron averaging of the path length (RF off, delta=0)")
Afast = 20e-3
function within_turn(Ax, Ay, nturns)
    z = user_to_canonical([Ax, 0.0, Ay, 0.0, 0.0, 0.0], ref, C, HARM)
    s = [0.0]; L = [0.0]; st = 0.0
    for t in 1:nturns, e in lat.elements
        z = exact_track_element(z, e, ref, C); st += (hasproperty(e, :L) ? e.L : 0.0)
        push!(s, st); push!(L, -z[5])     # path excess ΔL = -l at v=v0
    end
    return s, L
end
sX, LX = within_turn(Afast, 0.0, 3)
JxF = csJ(Afast, po.x); LamX = formula_B_Lambda(C, JxF, xi.x[2], 0, 0)
sY, LY = within_turn(0.0, Afast, 3)
JyF = csJ(Afast, po.y); LamY = formula_B_Lambda(C, 0, 0, JyF, xi.y[2])
pF1 = plot(xlabel = "s, м", ylabel = "ΔL(s), мм", legend = :topleft, size = (760, 440))
plot!(pF1, sX, LX .* 1e3, lw = 1.5, label = "трекинг M3, A_x = 20 мм")
plot!(pF1, [0, 3C], [0, 3C * LamX] .* 1e3, ls = :dash, lw = 2, label = "Λ_x s (M1)")
plot!(pF1, sY, LY .* 1e3, lw = 1.5, label = "трекинг M3, A_y = 20 мм")
plot!(pF1, [0, 3C], [0, 3C * LamY] .* 1e3, ls = :dash, lw = 2, label = "Λ_y s (M1)")
vline!(pF1, [C, 2C], color = :gray, ls = :dot, label = "")
savefig2(pF1, "fig_fast_within_turn")
# turn-by-turn path excess and its running mean
Nf = 400
function per_turn_L(Ax, Ay, N)
    z = user_to_canonical([Ax, 0.0, Ay, 0.0, 0.0, 0.0], ref, C, HARM); out = zeros(N)
    for n in 1:N
        l0 = z[5]; z = exact_track_turn(z, lat, ref, C); out[n] = -(z[5] - l0) / C
    end
    return out
end
LnX = per_turn_L(Afast, 0.0, Nf); LnY = per_turn_L(0.0, Afast, Nf)
runmean(v) = cumsum(v) ./ (1:length(v))
pF2 = plot(xlabel = "номер оборота n", ylabel = "Λ · 10⁶", legend = :topright, size = (760, 440))
scatter!(pF2, 1:Nf, LnX .* 1e6, ms = 2, msw = 0, alpha = 0.5, label = "Λ_n (x), за один оборот")
plot!(pF2, 1:Nf, runmean(LnX) .* 1e6, lw = 2, label = "скользящее среднее (x)")
hline!(pF2, [LamX * 1e6], ls = :dash, lw = 2, label = "M1 (x)")
scatter!(pF2, 1:Nf, LnY .* 1e6, ms = 2, msw = 0, alpha = 0.5, label = "Λ_n (y)")
plot!(pF2, 1:Nf, runmean(LnY) .* 1e6, lw = 2, label = "скользящее среднее (y)")
hline!(pF2, [LamY * 1e6], ls = :dash, lw = 2, label = "M1 (y)")
savefig2(pF2, "fig_fast_turns")
@printf("  Ax=20mm: Lambda note=%.5e, <Lambda_n>_%d=%.5e, single-turn spread=%.2e\n", LamX, Nf, mean(LnX), std(LnX))
@printf("  Ay=20mm: Lambda note=%.5e, <Lambda_n>_%d=%.5e, single-turn spread=%.2e\n", LamY, Nf, mean(LnY), std(LnY))
mac("LamXfast", e6(LamX)); mac("LamXfastTr", e6(mean(LnX))); mac("LamXspread", e6(std(LnX)))
mac("LamYfast", e6(LamY)); mac("LamYfastTr", e6(mean(LnY))); mac("LamYspread", e6(std(LnY)))
mac("JxFast", @sprintf("%.4f", JxF * 1e6)); mac("JyFast", @sprintf("%.4f", JyF * 1e6)); mac("Nf", string(Nf))

# ---------------------------------------------------------------- Slow synchrotron averaging
section("4. Slow synchrotron averaging: running mean of delta")
NT = 12000
cases = [("δ_m=1·10⁻³", 1e-3, 0.0, 0.0), ("δ_m=1·10⁻³, A_x=10 мм", 1e-3, 10e-3, 0.0),
         ("δ_m=1·10⁻³, A_x=20 мм", 1e-3, 20e-3, 0.0), ("δ_m=1·10⁻³, A_y=10 мм", 1e-3, 0.0, 10e-3),
         ("δ_m=1·10⁻³, A_y=20 мм", 1e-3, 0.0, 20e-3), ("δ_m=1·10⁻³, A_x=A_y=10 мм", 1e-3, 10e-3, 10e-3)]
res = []
for (lab, dm, Ax, Ay) in cases
    z0 = init(dm, Ax, Ay)
    Z3 = (Ax == 0 && Ay == 0 && haskey(seriesA, dm)) ? nothing : run(step_M3, z0, NT)
    d = Z3 === nothing ? seriesA[dm] : deltas(Z3)
    Zx = Z3 === nothing ? nothing : Z3
    r = Z3 === nothing ? analyse(run(step_M3, z0, 200); z0 = z0) : analyse(Z3; Ax = Ax, Ay = Ay, z0 = z0)
    m = measure_synchrotron(d, d; Qs = Qs)
    push!(res, (lab = lab, dm = dm, Ax = Ax, Ay = Ay, d = d, J0 = r.J0, m2 = nf6d_delta_eq(nf, r.J0),
                m1 = M1(m.delta_m_fit, csJ(Ax, po.x), csJ(Ay, po.y)), m3 = m.delta_eq, err = m.err_split, Qm = m.Qs_meas, Z = Zx))
    @printf("  %-28s M1=%+.5e M2=%+.5e M3=%+.5e  run-mean(N)=%+.5e\n", lab, res[end].m1, res[end].m2, m.delta_eq, mean(d))
end
open(joinpath(GEN, "tabConv.tex"), "w") do io
    for c in res
        d = c.d; rmv = runmean(d)
        Ts = 1 / c.Qm
        errs = [abs(rmv[N] - c.m2) for N in (100, 1000, 10000)]
        @printf(io, "%s & %.4f & %.4f & %.5f & %s & %s & %s & %.2e & %.2e & %.2e\\\\\n",
                replace(c.lab, "δ_m=1·10⁻³" => "\$\\delta_m{=}10^{-3}\$", "A_x" => "\$A_x\$", "A_y" => "\$A_y\$", "мм" => "мм"),
                c.J0[1] * 1e6, c.J0[2] * 1e6, c.J0[3] * 1e6, e6(c.m1), e6(c.m2), e6(c.m3), errs...)
    end
end
mac("NT", string(NT))
# (a) time series and running mean, base case
c0 = res[1]; nshow = 400
pT1 = plot(xlabel = "номер оборота n", ylabel = "δ_n · 10³", legend = :topright, size = (760, 380))
plot!(pT1, 1:nshow, c0.d[1:nshow] .* 1e3, lw = 1.2, marker = :circle, ms = 1.5, label = "δ_n, трекинг M3 (δ_m = 10⁻³, J_x=J_y=0)")
hline!(pT1, [0], color = :black, label = "")
savefig2(pT1, "fig_T_series")
pT2 = plot(xlabel = "номер оборота N", ylabel = "⟨δ⟩_N · 10⁶", legend = :outerright, size = (900, 440), xscale = :log10, ylims = (-8, 18))
for (ic, c) in enumerate(res)
    plot!(pT2, 1:length(c.d), runmean(c.d) .* 1e6, lw = 0.6, alpha = 0.5, color = ic, label = "")
    Tsc = 1 / c.Qm; Kc = floor(Int, length(c.d) / Tsc); Nkc = [round(Int, k * Tsc) for k in 1:Kc]
    plot!(pT2, Nkc, [mean(c.d[1:n]) for n in Nkc] .* 1e6, lw = 2.5, color = ic, label = c.lab)
    hline!(pT2, [c.m2 * 1e6], ls = :dash, color = ic, lw = 1, label = "")
end
savefig2(pT2, "fig_T_runmean")
# (b) dm scan running mean, normalised to theory
pT3 = plot(xlabel = "номер оборота N", ylabel = "⟨δ⟩_N / δ_eq(M2)", legend = :topright, size = (760, 420), xscale = :log10, ylims = (-2, 4))
for dm in dmsA
    d = seriesA[dm]; z0 = init(dm, 0, 0); th = nf6d_delta_eq(nf, nf6d_actions(nf, z0))
    plot!(pT3, 1:length(d), runmean(d) ./ th, lw = 1.5, label = @sprintf("δ_m = %.2f·10⁻³", dm * 1e3))
end
hline!(pT3, [1.0], ls = :dash, color = :black, lw = 2, label = "M2 (= M1)")
savefig2(pT3, "fig_T_runmean_dm")
# (c) convergence |<δ>_N - δ_eq| log-log with envelope
pT4 = plot(xlabel = "номер оборота N", ylabel = "|⟨δ⟩_N − δ_eq(M2)|", legend = :bottomleft, size = (760, 460),
           xscale = :log10, yscale = :log10)
c0 = res[1]; d = c0.d; th = c0.m2; Q = c0.Qm
Ns = 1:length(d)
plot!(pT4, Ns, max.(abs.(runmean(d) .- th), 1e-14), lw = 1, label = "простое среднее, δ_m = 10⁻³")
plot!(pT4, Ns, c0.dm ./ (Ns .* sin(pi * Q)), lw = 2, ls = :dash, label = "огибающая δ_m /(N sin πQ_s)")
Ts = 1 / Q
Kmax = floor(Int, length(d) / Ts)
Nk = [round(Int, k * Ts) for k in 1:Kmax]
plot!(pT4, Nk, max.(abs.([mean(d[1:n]) for n in Nk] .- th), 1e-14), lw = 2, marker = :circle, ms = 2, label = "среднее по целому числу периодов T_s")
Nfit = unique(round.(Int, exp10.(range(log10(300), log10(length(d)); length = 25))))
fitv = [abs(measure_synchrotron(d[1:n], d[1:n]; Qs = Qs, n_drop = 0).delta_eq - th) for n in Nfit]
plot!(pT4, Nfit, max.(fitv, 1e-14), lw = 2, marker = :square, ms = 3, label = "гармонический МНК-фит")
hline!(pT4, [abs(th) * 1e-2], ls = :dot, color = :black, label = "1 % от δ_eq")
savefig2(pT4, "fig_T_convergence")
N1 = findfirst(n -> all(abs.(runmean(d)[n:end] .- th) .< 0.01 * abs(th)), 1:length(d))
mac("NonePct", N1 === nothing ? ">$(length(d))" : string(N1))
mac("NsimpleOnePct", tex(c0.dm / (sin(pi * Q) * 0.01 * abs(th)); d = 1))
ifit = findfirst(i -> all(fitv[i:end] .< 0.01 * abs(th)), eachindex(fitv))
mac("NfitOnePct", ifit === nothing ? "—" : string(Nfit[ifit]))
pmerr = abs.([mean(d[1:n]) for n in Nk] .- th)
mac("PmErrEnd", tex(pmerr[end]; d = 1)); mac("SimpleErrEnd", tex(abs(runmean(d)[end] - th); d = 1))
mac("FitErrEnd", tex(fitv[end]; d = 1))
mac("Ts", @sprintf("%.1f", Ts))
# (d) spectrum
cS = res[end]; dS = cS.d .- mean(cS.d); NS = length(dS)
win = 0.5 .- 0.5 .* cos.(2pi .* (0:NS-1) ./ (NS - 1))
F = abs.([sum(dS .* win .* cis.(-2pi * k * (0:NS-1) / NS)) for k in 0:NS÷2]) .* (2 / sum(win))
fr = (0:NS÷2) ./ NS
pT5 = plot(xlabel = "частота, доля частоты обращения", ylabel = "амплитуда гармоники δ", yscale = :log10,
           legend = :topright, size = (760, 440), ylims = (1e-12, 1e-2))
plot!(pT5, fr[2:end], max.(F[2:end], 1e-13), lw = 1, label = "спектр δ_n, " * cS.lab)
for (q, nm) in ((Qs, "Q_s"), (2Qs, "2Q_s"), (3Qs, "3Q_s"), (min(Qlin[1], 1 - Qlin[1]), "Q_x"), (min(Qlin[2], 1 - Qlin[2]), "Q_y"))
    vline!(pT5, [abs(q)], ls = :dot, label = nm)
end
savefig2(pT5, "fig_T_spectrum")
println("  spectrum peak list (top 6):")
idx = sortperm(F[2:end], rev = true)[1:6] .+ 1
for k in idx
    @printf("    f=%.5f amp=%.3e\n", fr[k], F[k])
end
# zoom of the spectrum near the betatron tunes
qx = min(Qlin[1], 1 - Qlin[1]); qy = min(Qlin[2], 1 - Qlin[2])
pT6 = plot(xlabel = "частота, доля частоты обращения", ylabel = "амплитуда гармоники δ", yscale = :log10,
           legend = :topleft, size = (760, 420), xlims = (0.40, 0.50), ylims = (1e-13, 1e-8))
sel = findall(f -> 0.40 <= f <= 0.50, fr)
plot!(pT6, fr[sel], max.(F[sel], 1e-14), lw = 1, label = "спектр δ_n, " * cS.lab)
vline!(pT6, [qx, qy], color = :red, ls = :dash, label = "Q_x, Q_y")
vline!(pT6, vcat([qx + k * Qs for k in (-3, -2, -1, 1, 2, 3)]), color = :gray, ls = :dot, label = "Q_x ± kQ_s")
savefig2(pT6, "fig_T_spectrum_beta")
betamax = maximum(F[sel]); mac("betaLineMax", tex(betamax; d = 1))
@printf("  max betatron-line amplitude in delta: %.2e\n", betamax)
# ---------------------------------------------------------------- M4 long term vs M3
section("5. Long-term mean: M4 truncated map vs M3")
let dm = 1e-3, z0 = init(1e-3, 0, 0)
    Z4 = run(step_M4, z0, NA)
    d4 = deltas(Z4); d3 = seriesA[dm]
    th = nf6d_delta_eq(nf, nf6d_actions(nf, z0))
    pS = plot(xlabel = "номер оборота N", ylabel = "⟨δ⟩_N · 10⁶", legend = :topright, size = (760, 420), xscale = :log10, ylims = (-5, 8))
    plot!(pS, 1:NA, runmean(d3) .* 1e6, lw = 1.5, label = "M3: поэлементный трекинг")
    plot!(pS, 1:NA, runmean(d4) .* 1e6, lw = 1.5, ls = :dash, label = "M4: усечённая карта 3-го порядка")
    hline!(pS, [th * 1e6], color = :black, lw = 2, ls = :dot, label = "M2")
    savefig2(pS, "fig_M3_M4")
    @printf("  N=%d: <d>M3=%.6e <d>M4=%.6e M2=%.6e; amp drift M4: first/last 1000 rms %.4e %.4e\n", NA, mean(d3), mean(d4), th,
            std(d4[1:1000]), std(d4[end-999:end]))
    Tsl = 1 / Qs; Kl = floor(Int, NA / Tsl); nl = round(Int, Kl * Tsl)
    mac("MthreeMean", e6(mean(d3[1:nl]))); mac("MfourMean", e6(mean(d4[1:nl])))
    mac("MthreeFourDiff", tex(mean(d4[1:nl]) - mean(d3[1:nl]); d = 1)); mac("MthreeFourK", string(Kl))
    m3f = measure_synchrotron(d3, d3; Qs = Qs); m4f = measure_synchrotron(d4, d4; Qs = Qs)
    mac("MthreeFit", e6(m3f.delta_eq)); mac("MfourFit", e6(m4f.delta_eq))
end

# ---------------------------------------------------------------- write macros + summary
open(joinpath(GEN, "numbers.tex"), "w") do io
    for k in sort(collect(keys(macros)))
        println(io, "\\newcommand{\\v", k, "}{", macros[k], "}")
    end
end
open(joinpath(OUT, "v08_summary.txt"), "w") do io
    println(io, "v0.8 summary (thin RF everywhere)")
    @printf(io, "C=%.6f eta0=%.6f eta1=%.6f xi=(%.4f,%.4f) Q=(%.6f,%.6f,%.6f)\n", C, pe.eta0, pe.eta1, xi.x[2], xi.y[2], Qlin...)
    @printf(io, "NF d=(%.8e,%.8e,%.8e) note d=(%.8e,%.8e,%.8e) z*6=%.2e\n", nf.d..., d_note..., nf.zstar[6])
    println(io, "Test A: dm0 dmfit phat M1 M1p M2 M3 err M4")
    for row in rowsA
        r = row.r
        @printf(io, "%.2e %.5e %.4f %+.6e %+.6e %+.6e %+.6e %.1e %+.6e\n", row.dm0, r.m.delta_m_fit, r.m.phi_hat, r.m1, r.m1p, r.m2, r.m.delta_eq, r.m.err_split, row.m4)
    end
    println(io, "Test B: plane A Jcs Jnf Lnote LNF Lfit M1 M2 M3")
    for b in rowsB
        r = b.r
        @printf(io, "%s %.0e %.4e %.4e %+.5e %+.5e %+.5e %+.6e %+.6e %+.6e\n", b.pl, b.A, b.pl == :X ? r.Jx : r.Jy, r.J0[b.pl == :X ? 1 : 2], b.Lnote, b.LNF, b.Lfit, r.m1, r.m2, r.m.delta_eq)
    end
end
println("Done.")
flush(stdout); exit(0)
