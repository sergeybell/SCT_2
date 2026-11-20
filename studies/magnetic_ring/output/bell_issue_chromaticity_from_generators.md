# Magnetic chromaticity is an invariant of the RF-off Lie factor

## Summary

The magnetic one-turn factor has two elliptic planes and a longitudinal shear. `symplectic_eigen_A` asks for three elliptic pairs. `resonant_cubic` treats a cubic generator as having no pure invariant; that statement fits three oscillators. On the magnetic factor the cubic invariants include `J_x δ`, `J_y δ`, and `δ³`. Those `J δ` coefficients are the first chromaticity.

## What the magnetic ring shows

Lattice: `magnetic_2` with every sector bend named BM set to π/8, sextupoles at `k2 = 0`, no cavity. The closed file used here is `studies/magnetic_ring/lattices/magnetic_2_closed.bell`. The factor is `ring_lie_factor(magnetic_lattice(lat))`.

Elliptic pairs with Im λ > 1e-8: 2.
Eigenvalues: -9.712144e-01-2.382071e-01im, -9.712144e-01+2.382071e-01im, -9.706386e-01-2.405425e-01im, -9.706386e-01+2.405425e-01im, +1.000000e+00+0.000000e+00im, +1.000000e+00+0.000000e+00im.
symplectic_eigen_A: threw: one-turn matrix not stable

A normalising matrix built as a symplectic dispersion shear (with the `l` row `l → l − D' x + D p_x`) times the Courant–Snyder matrix at the ring start has symplectic error 6.67123917e-21. After the pull, `max |R_n[1:4, δ]|` is 8.67361738e-19 and the transverse blocks match the rotations to 1.87627691e-14.

Monomials that contain `l` have max |coefficient| 0.00000000e+00 in `f3` and 0.00000000e+00 in `f4`.

Coefficients of `J δ`: horizontal 9.18309974e+00, vertical 1.39831342e+01. Against `chromaticity_orders` on this ring the conversion is ξ1 = -1/(2π) × c, relative residual 1.304e-10, accepted: true.

The transverse homological operator, with `δ` held fixed, has kernel dimension 8 at the cutoff used (analytic count is 8: `δ³`, `J δ`, and the same families times powers of `l`). `F3` solves `F3 − F3 ∘ r = −f3` on the non-resonant cubic, residual 3.50015475e-15.

Coefficients of `J δ²` in `f4 − ½[f3, F3]`: horizontal -2.68750545e+01, vertical -2.54588110e+01. In `f4 + ½[f3, F3]`: horizontal -2.65654545e+01, vertical -2.50448183e+01. With the same factor −1/(2π), the minus combination reproduces `d²Q/dδ²`. The plus combination does not: once a sextupole makes the bracket the leading term, the plus combination flips the sign of ξ2.

Second chromaticity: f4 − ½[f3, F3] reproduces d²Q/dδ² with the ξ1 factor, on the bare ring and with one sextupole. f4 + ½[f3, F3] is stored as xi2_plus and is not the reported value.

## Recommendation

Do not pass the RF-off factor to `symplectic_eigen_A`. Document that the magnetic `(l, δ)` block is a shear. For that factor, read the first chromaticity from the invariant `J δ` in the cubic generator after a dispersion shear and a transverse Twiss normalisation. Read the second from `J δ²` in `f4 − ½[f3, F3]`, where `F3` solves `F3 − F3 ∘ r = −f3` on the transverse rotation and `δ` is a parameter. The plus sign in front of that bracket, as written for the three-oscillator normal form, flips the sextupole chromaticity on this ring. State that the kernel of the transverse homological operator includes `δ³` and `J δ`. `resonant_cubic`, as written, is the three-oscillator statement and drops this invariant.

## Reproducer

```julia
using Bell
using Bell.NormalForm: symplectic_eigen_A
lat = load_structure(raw"C:/GIT_REPOS/SCT_2/studies/magnetic_ring/lattices/magnetic_2_closed.bell")
F = ring_lie_factor(magnetic_lattice(lat))
symplectic_eigen_A(F.R)   # three elliptic pairs are required; this R has two
```

The same lattice is `examples/lattices/magnetic_2.bell` with the BM angle set to `π/8` and with no cavity. SCT_2 study step `studies/magnetic_ring/scripts/02b_generator_chromaticity.jl` records the generator coefficients. This note recommends an API reading. It does not change BELL.
