# da_fast.jl - PROTOTYPE optimisations of Bell.DiffAlgebra, kept in SCT_2 (BELL core untouched).
#   mul_to!(c, a, b, k)  in-place truncated product to order k, loops only over nonzero coefficients
#                        of b and over degree blocks that can contribute
#   compose_fast(M, N)   composition with in-place monomial table and in-place accumulation
#   inverse_fast(M)      fixed-point inversion, iteration k computed at order k only
#   daflow(F, X, L)      flow of an autonomous ODE z' = F(z) by the Lie series
#                        z(L) = sum_k L^k/k! (F.grad)^k id   (requires F(0) = 0)
using Bell.DiffAlgebra
const DAM = Bell.DiffAlgebra

nzidx(a::DA) = findall(!iszero, a.c)

"c = a*b truncated to order k (k <= no); c is overwritten."
function mul_to!(c::Vector{Float64}, a::DA, b::DA, k::Int = a.sp.no; bnz = nzidx(b))
    sp = a.sp; ds = sp.degstart; key = sp.key; deg = sp.deg
    fill!(c, 0.0)
    @inbounds for i in 1:ds[k+2]-1
        ai = a.c[i]; iszero(ai) && continue
        jmax = ds[k-deg[i]+2] - 1; ki = key[i]
        for j in bnz
            j > jmax && break
            c[DAM._index(sp, ki + key[j])] += ai * b.c[j]
        end
    end
    c
end

function compose_fast(M::AbstractVector{<:DA}, N::AbstractVector{<:DA}; k::Int = N[1].sp.no)
    spM = M[1].sp; spN = N[1].sp; nm = nmonomials(spM); len = nmonomials(spN)
    V = Vector{Vector{Float64}}(undef, nm); V[1] = zeros(len); V[1][1] = 1.0
    Nnz = [nzidx(n) for n in N]
    out = [zeros(len) for _ in M]
    # degree of each monomial of spM: if deg > k and N(0) = 0 it contributes nothing
    n0 = all(iszero(n.c[1]) for n in N)
    for m in 2:nm
        if n0 && spM.deg[m] > k
            V[m] = Float64[]; continue
        end
        V[m] = mul_to!(zeros(len), DA(spN, V[spM.parent[m]]), N[spM.pvar[m]], k; bnz = Nnz[spM.pvar[m]])
    end
    for (o, Mi) in zip(out, M), m in 1:nm
        c = Mi.c[m]; (iszero(c) || isempty(V[m])) && continue
        v = V[m]; @inbounds @simd for i in eachindex(v); o[i] += c * v[i]; end
    end
    [DA(spN, o) for o in out]
end

function inverse_fast(M::AbstractVector{<:DA})
    sp = M[1].sp; nv = sp.nv
    L = linear_part(M); Li = inv(L); I_ = identity_map(sp)
    NL = [Mi - sum(L[i, j] * I_[j] for j in 1:nv) for (i, Mi) in enumerate(M)]
    X = [sum(Li[i, j] * I_[j] for j in 1:nv) for i in 1:nv]
    for k in 2:sp.no              # after iteration k, X is exact through order k
        C = compose_fast(NL, X; k = k)
        R = [I_[i] - C[i] for i in 1:nv]
        X = [truncate_da(sum(Li[i, j] * R[j] for j in 1:nv), k) for i in 1:nv]
    end
    X
end

"""
    daflow(F, X, L; nstate) -> Vector{DA}

Taylor (Lie) series of the flow of the autonomous system z' = F(z) over length L,
F evaluated on the identity X (DA vars 1..nstate are the state; others are parameters).
Exact to roundoff in nDv if F(0) = 0 (no loss of order by the derivative).
"""
function daflow(F::AbstractVector{<:DA}, X::AbstractVector{<:DA}, L::Real; nstate = length(F), kmax = 200)
    out = similar(X, nstate)
    for i in 1:nstate
        g = X[i]; acc = g; small = 0
        for k in 1:kmax
            g = (L / k) * sum(F[j] * deriv(g, j) for j in 1:nstate)
            acc += g
            small = coeffnorm(g) < 1e-17 * coeffnorm(acc) ? small + 1 : 0
            small >= 2 && break
        end
        out[i] = acc
    end
    out
end