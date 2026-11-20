# Chromaticity of the magnetic one-turn Lie factor.
#
# The factor is ring_lie_factor of the RF-off lattice. Its linear part is a
# betatron rotation plus a longitudinal shear, so the normalising matrix is
# built here: a symplectic dispersion shear (including the l row), then the
# Courant–Snyder matrix of each transverse plane at the ring start.
# Actions in those coordinates are J = (q² + p²) / 2.
# The first chromaticity is the coefficient of J δ in the pulled cubic.
# The second is the coefficient of J δ² in the quartic after F3 removes the
# non-resonant cubic. δ is a parameter of the transverse homological
# operator. The kernel of that operator contains δ³ and J δ.
# Both bracket signs are returned. On this magnetic factor the combination
# that reproduces d²Q/dδ² is f4 − ½[f3, F3].

module GeneratorChromaticity

using LinearAlgebra
using Bell

export dispersion_shear, twiss_normalisation, normalisation_matrix,
       max_l_coefficient, action_coefficient, transverse_rotation,
       remove_nonresonant_cubic, generator_chromaticity,
       with_one_sextupole, sextupole_rows, feeddown_xi, feeddown_xi_integral

"""
    dispersion_shear(D4) -> Matrix{Float64}

Map from dispersion-free coordinates to canonical coordinates at one location.
`D4 = (D_x, D_x', D_y, D_y')`.

```
x  = x_β + D_x δ
px = px_β + D_x' δ
y  = y_β + D_y δ
py = py_β + D_y' δ
l  = l_β - D_x' x_β + D_x px_β - D_y' y_β + D_y py_β
δ  = δ
```

The `l` row is the one that keeps `dx∧dpx + dy∧dpy + dl∧dδ`.
"""
function dispersion_shear(D4)
    length(D4) == 4 || error("D4 must have four components, got $(length(D4))")
    Dx, Dpx, Dy, Dpy = Float64.(D4)
    A = Matrix{Float64}(I, 6, 6)
    A[index_x, index_delta] = Dx
    A[index_px, index_delta] = Dpx
    A[index_y, index_delta] = Dy
    A[index_py, index_delta] = Dpy
    A[index_l, index_x] = -Dpx
    A[index_l, index_px] = Dx
    A[index_l, index_y] = -Dpy
    A[index_l, index_py] = Dy
    return A
end

function _cs_block(beta, alpha)
    beta > 0 || error("beta must be positive, got $beta")
    sβ = sqrt(float(beta))
    return [sβ 0.0; -float(alpha) / sβ  1 / sβ]
end

"""
    twiss_normalisation(βx, αx, βy, αy) -> Matrix{Float64}

Courant–Snyder matrix at one location. Columns are normal coordinates.
Each plane is `[√β  0; −α/√β  1/√β]`. The `(l, δ)` block is the identity.
`J = (q² + p²) / 2` is the Courant–Snyder action.
"""
function twiss_normalisation(beta_x, alpha_x, beta_y, alpha_y)
    A = Matrix{Float64}(I, 6, 6)
    A[1:2, 1:2] = _cs_block(beta_x, alpha_x)
    A[3:4, 3:4] = _cs_block(beta_y, alpha_y)
    return A
end

"""
    normalisation_matrix(D4, βx, αx, βy, αy) -> Matrix{Float64}

`A = A_shear A_Twiss`. `z_canonical = A z_normal`, so
`pull_generator(f, A)` is `f` in normal coordinates.
"""
function normalisation_matrix(D4, beta_x, alpha_x, beta_y, alpha_y)
    return dispersion_shear(D4) * twiss_normalisation(beta_x, alpha_x, beta_y, alpha_y)
end

"""Largest |coefficient| of monomials that contain `l`."""
function max_l_coefficient(f)
    m = 0.0
    for k in 0:JET_ORDER
        for e in exponents_of_order(k)
            e[index_l] == 0 && continue
            m = max(m, abs(coeff(f, e)))
        end
    end
    return m
end

function _J(plane::Integer)
    i = 2 * plane - 1
    return (Z[i]^2 + Z[i + 1]^2) / 2
end

"""
    action_coefficient(f, plane, δ_power) -> NamedTuple

Coefficient `c` of `J δ^n` in `f`, with `J = (q² + p²) / 2`.
`c = coeff(q² δ^n) + coeff(p² δ^n)`. Also returns the two quadratic
coefficients and the `q p δ^n` coefficient.
"""
function action_coefficient(f, plane::Integer, delta_power::Integer)
    plane == 1 || plane == 2 || error("plane must be 1 or 2, got $plane")
    delta_power >= 0 || error("delta power must be non-negative")
    i = 2 * plane - 1
    eq = zeros(Int, 6)
    ep = zeros(Int, 6)
    ec = zeros(Int, 6)
    eq[i] = 2
    ep[i + 1] = 2
    ec[i] = 1
    ec[i + 1] = 1
    eq[index_delta] = delta_power
    ep[index_delta] = delta_power
    ec[index_delta] = delta_power
    q2 = coeff(f, eq)
    p2 = coeff(f, ep)
    cross = coeff(f, ec)
    return (c = q2 + p2, q2 = q2, p2 = p2, cross = cross)
end

"""Transverse rotation by `(μx, μy)`. The `(l, δ)` block stays the identity, so `δ` is a parameter."""
function transverse_rotation(mu_x::Real, mu_y::Real)
    R = Matrix{Float64}(I, 6, 6)
    for (k, mu) in ((1, float(mu_x)), (3, float(mu_y)))
        c, s = cos(mu), sin(mu)
        R[k:k + 1, k:k + 1] = [c s; -s c]
    end
    return R
end

function _monomial(e)
    p = Z[1]^0
    for i in 1:6
        e[i] == 0 && continue
        p *= Z[i]^e[i]
    end
    return p
end

function _coeffs(f, exps)
    return [coeff(f, e) for e in exps]
end

function _from_coeffs(x, exps)
    p = Z[1]^0 * 0
    for j in eachindex(exps)
        x[j] == 0 && continue
        p += x[j] * _monomial(exps[j])
    end
    return p
end

function _homological_matrix(exps, Rrot)
    n = length(exps)
    L = zeros(n, n)
    for j in 1:n
        mono = _monomial(exps[j])
        diffp = mono - pull_generator(mono, Rrot)
        for i in 1:n
            L[i, j] = coeff(diffp, exps[i])
        end
    end
    return L
end

"""
    remove_nonresonant_cubic(f3, μx, μy; rtol = 1e-9)

Solve `F3 − F3 ∘ r = −f3` on cubic polynomials, with `r` the transverse
rotation. Singular values below `rtol` times the largest are the kernel
(`δ³`, `J δ`, and the same monomials multiplied by powers of `l`).
The minimum-norm solution has no kernel component.
"""
function remove_nonresonant_cubic(f3, mu_x::Real, mu_y::Real; rtol::Real = 1.0e-9)
    Rrot = transverse_rotation(mu_x, mu_y)
    exps = exponents_of_order(3)
    L = _homological_matrix(exps, Rrot)
    x1 = action_coefficient(f3, 1, 1)
    y1 = action_coefficient(f3, 2, 1)
    cd = coeff(f3, [0, 0, 0, 0, 0, 3])
    f_res = x1.c * _J(1) * Z[index_delta] + y1.c * _J(2) * Z[index_delta] + cd * Z[index_delta]^3
    f_nr = hom(f3 - f_res, 3)
    rhs = -_coeffs(f_nr, exps)
    F = svd(L)
    smax = maximum(F.S)
    cutoff = rtol * smax
    keep = F.S .> cutoff
    any(keep) || error("homological operator has no invertible singular value")
    sinv = [keep[i] ? 1 / F.S[i] : 0.0 for i in eachindex(F.S)]
    x = F.V * (sinv .* (F.U' * rhs))
    F3 = hom(_from_coeffs(x, exps), 3)
    # Drop a round-off piece along J δ and δ³. It is in the kernel, so the
    # residual is unchanged, and it does not enter the Poisson bracket.
    sx = action_coefficient(F3, 1, 1)
    sy = action_coefficient(F3, 2, 1)
    sd = coeff(F3, [0, 0, 0, 0, 0, 3])
    F3 = hom(F3 - sx.c * _J(1) * Z[index_delta] - sy.c * _J(2) * Z[index_delta] - sd * Z[index_delta]^3, 3)
    residual_poly = hom(f_nr + F3 - pull_generator(F3, Rrot), 3)
    kept = F.S[keep]
    return (
        F3 = F3,
        residual = maxabs(residual_poly),
        F3_max = maxabs(F3),
        kernel_dimension = count(!, keep),
        smallest_kept = minimum(kept),
        largest_singular = smax,
        cutoff = cutoff,
    )
end

function _rotation_residual(Rn, mu_x, mu_y)
    rx = transverse_rotation(mu_x, mu_y)
    block = max(
        maximum(abs, Rn[1:2, 1:2] - rx[1:2, 1:2]),
        maximum(abs, Rn[3:4, 3:4] - rx[3:4, 3:4]),
    )
    coupling = maximum(abs, vcat(vec(Rn[1:2, 3:4]), vec(Rn[3:4, 1:2])))
    return (block = block, coupling = coupling)
end

"""
    generator_chromaticity(F, po) -> NamedTuple

Pull `F.f3` and `F.f4` through the normalising matrix of `po`
(`periodic_optics` at the ring start). Read `J δ` from the cubic.

`F3` solves `F3 − F3 ∘ r = −f3` on the non-resonant cubic. The quartic
is formed both ways:

- `x2`, `y2` from `f4 − ½[f3, F3]`
- `x2_plus`, `y2_plus` from `f4 + ½[f3, F3]`

`po` is the periodic optics of `F.R`. The numbers are coefficients, not
yet converted to `dQ/dδ`. The conversion factor is fixed by the caller
against `chromaticity_orders`.
"""
function generator_chromaticity(F::LieFactor, po)
    po.x.stable && po.y.stable || error("periodic optics is unstable")
    D4 = Float64.(collect(po.D4))
    A_shear = dispersion_shear(D4)
    A_twiss = twiss_normalisation(po.x.beta, po.x.alpha, po.y.beta, po.y.alpha)
    A = A_shear * A_twiss
    Rn = inv(A) * F.R * A
    Rs = inv(A_shear) * F.R * A_shear
    f3n = hom(pull_generator(F.f3, A), 3)
    f4n = hom(pull_generator(F.f4, A), 4)
    removal = remove_nonresonant_cubic(f3n, po.x.mu, po.y.mu)
    half_bracket = poisson(f3n, removal.F3) / 2
    f4_minus = hom(f4n - half_bracket, 4)
    f4_plus = hom(f4n + half_bracket, 4)
    rot = _rotation_residual(Rn, po.x.mu, po.y.mu)
    return (
        symplectic_error_A = symplectic_error(A),
        symplectic_error_shear = symplectic_error(A_shear),
        symplectic_error_twiss = symplectic_error(A_twiss),
        dispersion_column = maximum(abs, Rn[1:4, index_delta]),
        shear_dispersion_column = maximum(abs, Rs[1:4, index_delta]),
        longitudinal_coupling = maximum(abs, Rn[index_l, 1:4]),
        delta_row = maximum(abs, Rn[index_delta, 1:5]),
        delta_delta = Rn[index_delta, index_delta],
        l_delta = Rn[index_l, index_delta],
        rotation_residual = rot.block,
        transverse_coupling = rot.coupling,
        max_l_f3 = max_l_coefficient(f3n),
        max_l_f4 = max_l_coefficient(f4n),
        x1 = action_coefficient(f3n, 1, 1),
        y1 = action_coefficient(f3n, 2, 1),
        x2 = action_coefficient(f4_minus, 1, 2),
        y2 = action_coefficient(f4_minus, 2, 2),
        x2_plus = action_coefficient(f4_plus, 1, 2),
        y2_plus = action_coefficient(f4_plus, 2, 2),
        x2_bare = action_coefficient(f4n, 1, 2),
        y2_bare = action_coefficient(f4n, 2, 2),
        delta3 = coeff(f3n, [0, 0, 0, 0, 0, 3]),
        homological_residual = removal.residual,
        F3_max = removal.F3_max,
        kernel_dimension = removal.kernel_dimension,
        smallest_kept = removal.smallest_kept,
        largest_singular = removal.largest_singular,
        optics_coupling = po.coupling,
        mu_x = po.x.mu,
        mu_y = po.y.mu,
        beta_x = po.x.beta,
        alpha_x = po.x.alpha,
        beta_y = po.y.beta,
        alpha_y = po.y.alpha,
        D4 = D4,
    )
end

"""Copy `lattice` with the sextupole at `index` set to `k2`."""
function with_one_sextupole(lattice::Lattice, index::Integer, k2::Real)
    1 <= index <= length(lattice.elements) || error("element index $index is outside the ring")
    elems = AbstractElement[e for e in lattice.elements]
    e = elems[index]
    e isa Sextupole || error("element $index ($(e.name)) is not a sextupole")
    elems[index] = Sextupole(e.name, e.L, float(k2))
    return Lattice(lattice.name, elems, lattice.ref, lattice.harmonic)
end

function _sample_index(opt, s)
    k = findfirst(x -> abs(x - s) <= 1.0e-8, opt.data.s)
    k === nothing && error("s = $s m is not on the optics table")
    return k
end

"""β(s) and D_x(s) of a drift, integrated from the entrance Twiss."""
function _drift_feed(beta, alpha, D, Dp, L)
    beta > 0 || error("beta must be positive, got $beta")
    g = (1 + alpha^2) / beta
    b0 = beta * D
    b1 = beta * Dp - 2 * alpha * D
    b2 = -2 * alpha * Dp + g * D
    b3 = g * Dp
    integral = b0 * L + b1 * L^2 / 2 + b2 * L^3 / 3 + b3 * L^4 / 4
    s = L / 2
    return (
        integral = integral,
        beta_mid = beta - 2 * alpha * s + g * s^2,
        D_mid = D + Dp * s,
    )
end

"""
    sextupole_rows(lattice, opt) -> Vector{NamedTuple}

One row per sextupole. `beta` and `D` are the drift transport from the
element entrance to its centre. `D` is horizontal dispersion. The integrals
are `∫ β D_x ds` along the element.
"""
function sextupole_rows(lattice::Lattice, opt)
    rows = NamedTuple[]
    d = opt.data
    for (i, e) in enumerate(lattice.elements)
        e isa Sextupole || continue
        s0 = opt.bounds[i][1]
        k = _sample_index(opt, s0)
        hx = _drift_feed(d.betx[k], d.alfx[k], d.D[k], d.Dp[k], e.L)
        hy = _drift_feed(d.bety[k], d.alfy[k], d.D[k], d.Dp[k], e.L)
        push!(rows, (
            index = i,
            name = e.name,
            L = e.L,
            s0 = s0,
            k2 = e.k2,
            beta_x = hx.beta_mid,
            beta_y = hy.beta_mid,
            D = hx.D_mid,
            beta_x_D = hx.beta_mid * hx.D_mid,
            integral_beta_x_D = hx.integral,
            integral_beta_y_D = hy.integral,
        ))
    end
    return rows
end

"""`k2 L β D / (4π)` with the sign of the horizontal feed-down."""
feeddown_xi(k2::Real, L::Real, beta::Real, D::Real) = k2 * L * beta * D / (4π)

"""`k2 ∫ β D ds / (4π)`."""
feeddown_xi_integral(k2::Real, integral_beta_D::Real) = k2 * integral_beta_D / (4π)

end # module
