# cosy_io.jl - readers for COSY INFINITY output
"""
    read_cosy_da(file; nv = 6) -> Vector{Dict{Vector{Int},Float64}}

Reads every DA written by `WRITE unit X` in COSY ("I COEFFICIENT ORDER EXPONENTS" blocks).
Exponents are parsed by fixed width (2 characters per exponent, pairs separated by a blank),
so that two-digit exponents of order 10 are read correctly.
"""
function read_cosy_da(file; nv = 6)
    blocks = Dict{Vector{Int},Float64}[]; cur = nothing
    for ln in eachline(file)
        if occursin("COEFFICIENT", ln)
            cur = Dict{Vector{Int},Float64}(); push!(blocks, cur)
        elseif startswith(strip(ln), "---")
            cur = nothing
        elseif occursin("ALL COMPONENTS ZERO", ln)
            push!(blocks, Dict{Vector{Int},Float64}()); cur = nothing
        elseif cur !== nothing && !isempty(strip(ln))
            f = split(ln); coef = parse(Float64, replace(f[2], 'D' => 'E'))
            ordr = parse(Int, f[3])
            ex = _cosy_exps(ln, ordr, nv ÷ 2)
            cur[ex] = coef
        end
    end
    blocks
end
# exponent field: the last 15 characters, three groups (1X,I2,1X?,I2) = 5 chars " a  b" right-aligned
function _cosy_exps(ln, ordr, ng = 3)
    t = rstrip(ln); t = t[end-5ng+1:end]
    vals = Int[]
    for g in 0:ng-1
        s = t[5g+1:5g+5]
        push!(vals, parse(Int, s[1:3]), parse(Int, s[4:5]))
    end
    sum(vals) == ordr || error("exponent parse failed: $ln")
    vals
end
read_cosy_times(file) = Dict(split(l)[1] => parse(Float64, split(l)[3]) / parse(Float64, split(l)[2]) for l in eachline(file) if !isempty(strip(l)))