# fringe_comparison/BELL/run.jl -- all Bell fringe models on magnetic_2_md, machine-readable output
#   run_bell.bat  (= julia --project=. run.jl)      Bell is a library: Project.toml -> C:\GIT_REPOS\BELL
#
# GAP = half-aperture as in COSY; dipole full gap = 2*GAP.  Sextupoles: NoEdge everywhere.
# Widths of the non-Enge profiles are chosen so that the fringe integral I = int F(1-F) dz
# equals that of Enge(gap = 2*GAP) (same FINT => same leading-order vertical edge focusing).
# Tabulated = Enge(gap = 2*GAP) sampled on 403..409 points (checks the interpolation path).
#
# Map-order limits: Bell's ring map is a Lie factor with H to degree JET_ORDER = 4, i.e. the
# polynomial map is exact to order MAP_ORDER = 3 in (x,px,y,py,l,delta).  The Jacobian on the
# closed orbit is therefore exact to delta^2, so Q(delta) is exact to O(delta^2): xi1 and xi2 are
# reliable, the fitted delta^3, delta^4 coefficients are truncation artefacts (written, flagged).
using Bell, Printf, LinearAlgebra
const LAT = joinpath(@__DIR__, "..", "..", "BELL", "magnetic_2_md.bell")
const OUT = joinpath(@__DIR__, "output"); mkpath(OUT)
lat = load_structure(LAT)
GAPS = [0.05, 0.025, 0.0125, 0.00625]

# ---- profiles with the same fringe integral as Enge(2*GAP)
function match(mk, I0, a, b)          # bisection on the width parameter
    for _ in 1:60
        m = (a + b) / 2
        fringe_integral(mk(m)) < I0 ? (a = m) : (b = m)
    end
    (a + b) / 2
end
function profiles(gap)
    E = Enge(gap = 2gap); I0 = fringe_integral(E)
    lam = match(Tanh, I0, 1e-5, 1.0)
    w = match(LinearRamp, I0, 1e-5, 2.0)
    F, W, _ = profile_function(E)
    # NB Bell caches edge maps by repr(model) and repr(Tabulated) is only 'Tabulated(N points)':
    # tables of equal length for different gaps would hit the same cache entry (Bell bug, see README).
    # Workaround: a different number of points for every gap.
    z = collect(range(-W, W, length = 401 + 2 * findfirst(==(gap), GAPS)))
    T = Tabulated(z, F.(z))
    return [("Enge", E), ("Tanh", Tanh(lam)), ("LinearRamp", LinearRamp(w)), ("Tabulated", T)], I0, lam, w
end

# ---- Q(delta) on closed orbits, polynomial fit
const DGRID = collect(range(-4e-3, 4e-3, length = 17))
function qseries(M; deg = 4)
    Q = [tunes_at_delta(M, d) for d in DGRID]
    A = [d^k for d in DGRID, k in 0:deg]
    (A \ [q[1] for q in Q], A \ [q[2] for q in Q])
end

rows = String[]; series = String[]
push!(rows, "variant,dipole,quad,profile,method,gap,width_param,I_fringe,Qx,Qy,xi1x,xi1y,xi2x,xi2y,seconds")
push!(series, "variant,gap,k,qx_k,qy_k,reliable")
function run!(variant, dip, quad, prof, meth, gap, wpar, I0, model)
    t = time()
    F = ring_lie_factor(lat, model); M = factor_to_map(F)
    ch = chromaticity_orders(M; dd = 1e-4)
    cx, cy = qseries(M)
    dt = time() - t
    push!(rows, @sprintf("%s,%s,%s,%s,%s,%.6g,%.8g,%.8g,%.12f,%.12f,%.10f,%.10f,%.8f,%.8f,%.2f",
        variant, dip, quad, prof, meth, gap, wpar, I0, ch.x[1], ch.y[1], ch.x[2], ch.y[2], ch.x[3], ch.y[3], dt))
    for k in 0:length(cx)-1
        push!(series, @sprintf("%s,%.6g,%d,%.12e,%.12e,%d", variant, gap, k, cx[k+1], cy[k+1], k <= 2))
    end
    @printf("%-48s gap=%-8g Qx=%.6f Qy=%.6f xi_x=%.6f xi_y=%.6f  (%.1f s)\n", variant, gap, ch.x[1], ch.y[1], ch.x[2], ch.y[2], dt)
    flush(stdout)
end

# gap-independent models
run!("NoEdge", "NoEdge", "NoEdge", "-", "-", NaN, NaN, NaN, FringeModel())
for o in 1:3
    run!("HardEdge$o", "HardEdge($o)", "NoEdge", "-", "-", NaN, NaN, NaN, FringeModel(dipole = HardEdge(order = o)))
end
run!("Quad_HardEdge2", "NoEdge", "HardEdge(2)", "-", "-", NaN, NaN, NaN, FringeModel(quad = HardEdge(order = 2)))
run!("HardEdge2+Quad_HardEdge2", "HardEdge(2)", "HardEdge(2)", "-", "-", NaN, NaN, NaN,
     FringeModel(dipole = HardEdge(order = 2), quad = HardEdge(order = 2)))

methods_ = [("Splitting", () -> Splitting(n = 40, order = 4)), ("Magnus", () -> Magnus(n = 40, order = 4)),
            ("TaylorIntegrate", () -> TaylorIntegrate(tol = 1e-10, nmax = 512))]
for gap in GAPS
    profs, I0, lam, w = profiles(gap)
    @printf("# GAP = %g: I(Enge 2GAP) = %.6g m, Tanh lambda = %.6g m, LinearRamp width = %.6g m\n", gap, I0, lam, w)
    # hard edge with FINT of the Enge profile
    run!("HardEdge2_FINT", "HardEdge(2;FINT)", "NoEdge", "Enge", "-", gap, I0 / (2gap), I0,
         FringeModel(dipole = HardEdge(order = 2, gap = 2gap, fint = I0 / (2gap))))
    for (pn, p) in profs, (mn, mk) in methods_
        wpar = pn == "Enge" ? 2gap : pn == "Tanh" ? lam : pn == "LinearRamp" ? w : 2gap
        run!("Soft_$(pn)_$(mn)", "SoftEdge", "NoEdge", pn, mn, gap, wpar, I0,
             FringeModel(dipole = SoftEdge(profile = p, method = mk())))
    end
    run!("Soft_Enge_Splitting+Quad_HardEdge2", "SoftEdge", "HardEdge(2)", "Enge", "Splitting", gap, 2gap, I0,
         FringeModel(dipole = SoftEdge(profile = profs[1][2], method = Splitting(n = 40, order = 4)), quad = HardEdge(order = 2)))
end
write(joinpath(OUT, "bell_models.csv"), join(rows, "\n") * "\n")
write(joinpath(OUT, "bell_q_series.csv"), join(series, "\n") * "\n")
println("written: output/bell_models.csv, output/bell_q_series.csv")
