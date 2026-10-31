# Shared paths, plot style, and helpers for scripts/steps/*.jl
using EquilibriumEnergy
using Bell, Printf, Plots, Statistics, LinearAlgebra, TaylorSeries

gr()
default(fontfamily = "sans-serif", framestyle = :box, grid = true, dpi = 150,
        legendfontsize = 8, tickfontsize = 9, guidefontsize = 11)

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OPTIM = raw"C:\GIT_REPOS\SCT_2\OptiM\magnetic\magnetic_2.opt"
const OUT = joinpath(ROOT, "output")
const NOTE = joinpath(ROOT, "note")
const GEN = joinpath(NOTE, "gen")
const FIG = joinpath(NOTE, "figs")
mkpath(OUT); mkpath(GEN); mkpath(FIG)

section(t) = println("\n", "="^72, "\n ", t, "\n", "="^72)

function savefig_step(p, name)
    savefig(p, joinpath(OUT, name * ".png"))
    savefig(p, joinpath(FIG, name * ".png"))
    println("  wrote ", name, ".png")
end

texnum(x; d = 4) = begin
    x == 0 && return "0"
    e = floor(Int, log10(abs(x)))
    m = x / 10.0^e
    e == 0 ? @sprintf("%.*f", d, m) : @sprintf("%.*f\\cdot10^{%d}", d, m, e)
end

e6(x) = @sprintf("%.4f", x * 1e6)

mutable struct Ctx
    macros::Dict{String,String}
    series::Dict{Symbol,Any}
end
Ctx() = Ctx(Dict{String,String}(), Dict{Symbol,Any}())
mac!(ctx::Ctx, k, v) = (ctx.macros[k] = string(v))

function write_numbers(ctx::Ctx)
    open(joinpath(GEN, "numbers.tex"), "w") do io
        for k in sort(collect(keys(ctx.macros)))
            println(io, "\\newcommand{\\v", k, "}{", ctx.macros[k], "}")
        end
    end
end

function write_step_tex(name, body)
    open(joinpath(GEN, name * ".tex"), "w") do io
        write(io, body)
    end
end

function write_conclusion(name, text)
    open(joinpath(OUT, name * "_conclusion.txt"), "w") do io
        println(io, text)
    end
    println("  conclusion: ", text)
end
