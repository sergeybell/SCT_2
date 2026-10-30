# ode_compare.jl - (3) ODE integration with DA initial conditions: Bell.Integrators vs COSY
#   problems: hh     Henon-Heiles, nv = 4, z0 = (0.1,-0.2,0.25,0.1) + DA, T = 1
#             quad   MQ 0.5 m, 0.5 T pole / 0.05 m, proton 200 MeV, exact sqrt Hamiltonian, nv = 6
#             dipole MS R = 2 m, 30 deg, homogeneous sector, exact Hamiltonian, nv = 6
#   Bell: RK4, DP5, DOP853, GBS(6), Yoshida 4/6/8 (separable problems), DA-flow prototype (elements)
#   COSY: FOX RK4 (hh, identical algorithm), element maps MQ / MS (COSY's own DA flow).
# Needs cosy/output from cosy/run_cosy_ode.ps1.  Run: julia --project=julia julia/ode_compare.jl [orders...]
using Bell, Bell.DiffAlgebra, Bell.Integrators, Printf
include(joinpath(@__DIR__, "cosy_io.jl")); include(joinpath(@__DIR__, "da_fast.jl"))
const ROOT = dirname(@__DIR__); const OUT = joinpath(ROOT, "output"); const CO = joinpath(ROOT, "cosy", "output")
const ORDERS = isempty(ARGS) ? [3, 5, 7, 10] : parse.(Int, ARGS)
const CHIM, ETA = parse.(Float64, split(read(joinpath(CO, "refpart.txt"), String)))   # B rho [T m], K/mc^2

# ---------------------------------------------------------------- problems
struct Problem
    name::String; nv::Int; L::Float64
    z0::Function           # sp -> initial DA state (4 components)
    rhs::Function          # sp -> f(z, s)
    split::Union{Nothing,Function}   # sp -> SplitFlows
    flowF::Union{Nothing,Function}   # sp -> (F on identity, X) for daflow
end
hh = Problem("hh", 4, 1.0, sp -> [c + x for (c, x) in zip([0.1, -0.2, 0.25, 0.1], da_vars(sp))],
    sp -> (z, s) -> [z[3], z[4], -z[1] - 2z[1] * z[2], -z[2] - z[1] * z[1] + z[2] * z[2]],
    sp -> SplitFlows((z, h) -> [z[1] + h * z[3], z[2] + h * z[4], z[3], z[4]],
                     (z, h) -> [z[1], z[2], z[3] - h * (z[1] + 2z[1] * z[2]), z[4] - h * (z[2] + z[1] * z[1] - z[2] * z[2])]),
    nothing)
onepd(sp) = (dK = da_var(sp, 6); sqrt((1 + dK) * ((1 + dK) * ETA + 2) / (ETA + 2)))   # p/p0 from COSY dK
const KQ = (0.5 / 0.05) / CHIM
function quad_rhs(sp)
    P2 = onepd(sp)^2
    (z, s) -> (ip = inv(sqrt(P2 - z[2] * z[2] - z[4] * z[4])); [z[2] * ip, -KQ * z[1], z[4] * ip, KQ * z[3]])
end
function quad_split(sp)
    P2 = onepd(sp)^2
    SplitFlows((z, h) -> (ip = h * inv(sqrt(P2 - z[2] * z[2] - z[4] * z[4])); [z[1] + z[2] * ip, z[2], z[3] + z[4] * ip, z[4]]),
               (z, h) -> [z[1], z[2] - (h * KQ) * z[1], z[3], z[4] + (h * KQ) * z[3]])
end
quad = Problem("quad", 6, 0.5, sp -> da_vars(sp)[1:4], quad_rhs, quad_split,
    sp -> (X = da_vars(sp); (quad_rhs(sp)(X[1:4], 0.0), X)))
const HD = 0.5
function dip_rhs(sp)
    P = onepd(sp); P2 = P^2
    (z, s) -> (e = 1 + HD * z[1]; ps = sqrt(P2 - z[2] * z[2] - z[4] * z[4]); ip = e * inv(ps);
               [z[2] * ip, HD * (ps - e), z[4] * ip, 0 * z[4]])
end
dipole = Problem("dipole", 6, 2 * pi / 6, sp -> da_vars(sp)[1:4], dip_rhs, nothing,
    sp -> (X = da_vars(sp); (dip_rhs(sp)(X[1:4], 0.0), X)))

# ---------------------------------------------------------------- runs
function timed(f)                     # warm, then min over repeats for >= ~0.3 s
    r = f(); t = @elapsed f(); n = t < 0.3 ? clamp(ceil(Int, 0.3 / t), 1, 200) : 1
    for _ in 2:n; t = min(t, @elapsed f()); end
    r, t
end
cfg(p) = vcat([("RK4", RK4(), :n, n) for n in (4, 8, 16, 32, 64, 128)],
              [("DP5", DormandPrince5(), :tol, t) for t in (1e-6, 1e-8, 1e-10, 1e-12, 1e-14)],
              [("DOP853", DormandPrince853(), :tol, t) for t in (1e-8, 1e-10, 1e-12, 1e-14, 1e-15)],
              [("GBS(6)", GraggBulirschStoer(6), :tol, t) for t in (1e-8, 1e-10, 1e-12, 1e-14, 1e-15)],
              p.split === nothing ? [] :
              [(nm, m, :split, n) for (nm, m) in (("Yoshida4", Yoshida(4)), ("Yoshida6", Yoshida(6)), ("Yoshida8", Yoshida(8)))
                                  for n in (2, 4, 8, 16, 32, 64)])
function runone(p, sp, (nm, m, kind, par))
    z0 = p.z0(sp)
    kind == :n     && return integrate(p.rhs(sp), z0, 0.0, p.L; method = m, nsteps = par)
    kind == :tol   && return integrate(p.rhs(sp), z0, 0.0, p.L; method = m, tol = par)
    kind == :split && return integrate(p.split(sp), z0, 0.0, p.L; method = m, nsteps = par)
end
cosymap(file, sp) = [(d = Dict(b); DA(sp, [get(d, monomial_exponents(sp, m), 0.0) for m in 1:nmonomials(sp)])) for b in read_cosy_da(file; nv = sp.nv)][1:4]
relerr(z, ref) = maximum(coeffnorm(z[i] - ref[i]) for i in 1:4) / maximum(coeffnorm(r) for r in ref)

wp = open(joinpath(OUT, "ode_workprecision.csv"), "w")
println(wp, "problem,order,method,param,time_s,nsteps,nfev,err_vs_ref,err_vs_cosy")
cm = open(joinpath(OUT, "ode_cosy.csv"), "w")
println(cm, "problem,order,what,time_s,err_vs_ref")
for no in ORDERS, p in (hh, quad, dipole)
    sp = DASpace(p.nv, no)
    # reference: hh - GBS(8) at tol 1e-14 (agrees with DOP853 1e-15 to ~1e-14)
    # elements: the DA-flow series (exact to roundoff; agrees with COSY to ~1e-16) is the reference
    ref = p.flowF === nothing ?
        integrate(p.rhs(sp), p.z0(sp), 0.0, p.L; method = GraggBulirschStoer(8), tol = 1e-14, maxsteps = 20000).z :
        (FX = p.flowF(sp); daflow(FX[1], FX[2], p.L; nstate = 4))
    if p.name == "hh"
        cosy = cosymap(joinpath(CO, "hh_o$(no)_n320.txt"), sp)
        ctimes = Dict(parse(Int, split(l)[1]) => parse(Float64, split(l)[3]) / parse(Float64, split(l)[2]) for l in eachline(joinpath(CO, "hh_time_o$no.txt")))
        for (n, t) in sort(collect(ctimes))
            cz = cosymap(joinpath(CO, "hh_o$(no)_n$n.txt"), sp)
            @printf(cm, "hh,%d,COSY FOX RK4 n=%d,%.4e,%.3e\n", no, n, t, relerr(cz, ref))
        end
    else
        k = p.name == "quad" ? 1 : 2
        cosy = cosymap(joinpath(CO, "elem_o$(no)_$k.txt"), sp)
        l = [split(x) for x in eachline(joinpath(CO, "elem_time_o$no.txt"))][k]
        @printf(cm, "%s,%d,COSY element map,%.4e,%.3e\n", p.name, no, parse(Float64, l[3]) / parse(Float64, l[2]), relerr(cosy, ref))
        (F, X) = p.flowF(sp)
        zf, tf = timed(() -> daflow(p.flowF(sp)[1], X, p.L; nstate = 4))
        @printf(cm, "%s,%d,Bell DA-flow prototype,%.4e,%.3e\n", p.name, no, tf, relerr(zf, ref))
        @printf(wp, "%s,%d,DA-flow (prototype),0,%.4e,0,0,%.3e,%.3e\n", p.name, no, tf, relerr(zf, ref), relerr(zf, cosy))
    end
    for c in cfg(p)
        r, t = timed(() -> runone(p, sp, c))
        ns = hasproperty(r, :nsteps) ? r.nsteps : c[4]; nf = hasproperty(r, :nfev) ? r.nfev : 0
        @printf(wp, "%s,%d,%s,%g,%.4e,%d,%d,%.3e,%.3e\n", p.name, no, c[1], c[4], t, ns, nf, relerr(r.z, ref), relerr(r.z, cosy))
        flush(wp)
    end
    @printf("%-6s order %2d done\n", p.name, no); flush(cm)
end
close(wp); close(cm)