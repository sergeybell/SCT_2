# Chromaticity from the magnetic Lie generators

Lattice file: `C:\GIT_REPOS\SCT_2\studies\magnetic_ring\lattices\magnetic_2_closed.bell`. Reference: deuteron, kinetic energy 270.0 MeV, harmonic 1. No RF. The bare case has every sextupole at `k2 = 0`.

The one-turn object is `ring_lie_factor` of `magnetic_lattice`. `chromaticity_orders` reads `factor_to_map` of that factor. `chromaticity_series` builds `one_turn_map_da` at order 4. The normalising matrix is not `symplectic_eigen_A`.

## Normalisation

`A` maps normal coordinates to canonical coordinates: dispersion shear, then the Courant–Snyder matrix of each plane at the ring start. The shear is `x → x + D δ`, `p_x → p_x + D' δ`, and `l → l − D' x + D p_x`, with the same pattern in `y`. Actions are `J = (q² + p²) / 2`.

| quantity | k2 = 0 | one sextupole |
| --- | --- | --- |
| symplectic error of A | 6.67123917e-21 | 6.67123917e-21 |
| symplectic error of the shear | 2.03078623e-23 | 2.03078623e-23 |
| max |R_n[1:4, δ]| | 8.67361738e-19 | 8.67361738e-19 |
| max |R_shear[1:4, δ]| | 1.73472348e-18 | 1.73472348e-18 |
| max |R_n[l, 1:4]| | 3.80121282e-15 | 3.80121282e-15 |
| rotation residual | 1.87627691e-14 | 1.87627691e-14 |
| transverse coupling of R_n | 0.00000000e+00 | 0.00000000e+00 |
| max |R_n[δ, 1:5]| | 0.00000000e+00 | 0.00000000e+00 |
| (l|δ) in normal coordinates | 6.26487986e+01 | 6.26487986e+01 |
| max |coef| of f3 monomials with l | 0.00000000e+00 | 0.00000000e+00 |
| max |coef| of f4 monomials with l | 0.00000000e+00 | 0.00000000e+00 |
| |q² − p²| coefficient, x, Jδ | 2.75731501e-01 | 1.34913220e+01 |
| |q p| coefficient, x, Jδ | 1.21341809e-02 | 3.21357783e+00 |
| |q² − p²| coefficient, y, Jδ | 2.56772944e-01 | 3.02078843e+00 |
| homological residual | 3.50015475e-15 | 3.50067197e-14 |
| max |coef| of F3 | 2.28476779e+00 | 2.91100959e+01 |
| kernel dimension | 8 | 8 |
| smallest kept singular value | 2.40537126e-03 | 2.40537126e-03 |
| |Lie R − ring_matrix| | 2.70006240e-13 | 2.70006240e-13 |
| |periodic D − dispersion_orders| | 0.00000000e+00 | 0.00000000e+00 |

Normalisation checks are within 1.0e-8: true.

## First chromaticity

The cubic generator in normal coordinates contains `c_x J_x δ + c_y J_y δ`. The map is `exp(:f:)`, and the angle that matches `rotation` by `+μ` is the generator `−μ J`. The overall factor was required to be one entry of a fixed list, the same number on both planes, and was compared with `chromaticity_orders` at step `1e-6` on the `k2 = 0` ring.

| quantity | value |
| --- | --- |
| c_x of J_x δ | 9.18309974e+00 |
| c_y of J_y δ | 1.39831342e+01 |
| ξ1_x / c_x | -1.59154943e-01 |
| ξ1_y / c_y | -1.59154943e-01 |
| candidate | -1/(2π) |
| factor | -1.59154943e-01 |
| max relative residual | 1.304e-10 |
| accepted | yes |

Accepted factor: ξ1 = -1/(2π) × (coefficient of J δ). On this factor, `d²Q/dδ²` predicted from a term `c₂ J δ²` is `2 × -1/(2π) × c₂`.

`chromaticity_orders` is a central difference. Steps `1e-4` (the default), `1e-5`, and `1e-6` are all stored. The comparison with `chromaticity_series` uses the step `1e-6` once the three steps have settled. ξ1 is `dQ/dδ`. ξ2 is `d²Q/dδ²`, the third component of `chromaticity_orders` and the second derivative from `chromaticity_series`.

## Bare ring, k2 = 0

| plane | quantity | generator | orders 1e-4 | orders 1e-6 | series |
| --- | --- | --- | --- | --- | --- |
| x | ξ1 | -1.46153572e+00 | -1.46153752e+00 | -1.46153572e+00 | -1.46153572e+00 |
| x | ξ2 | 8.55459555e+00 | 8.55473511e+00 | 8.55454596e+00 | 8.55459555e+00 |
| y | ξ1 | -2.22548493e+00 | -2.22548794e+00 | -2.22548493e+00 | -2.22548493e+00 |
| y | ξ2 | 8.10379124e+00 | 8.10412260e+00 | 8.10368439e+00 | 8.10379124e+00 |

Generator ξ2 in that table uses `c₂` from `f4 − ½[f3, F3]` and the direct prediction `2 × factor × c₂`. It is not rescaled. `F3` is the solution of `F3 − F3 ∘ r = −f3` on the non-resonant cubic. The row `ξ2 plus` below is the same prediction with the bracket added instead of subtracted.

| plane | ξ2 from f4 − ½[f3,F3] | ξ2 from f4 + ½[f3,F3] | orders 1e-6 | series |
| --- | --- | --- | --- | --- |
| x | 8.55459555e+00 | 8.45604679e+00 | 8.55454596e+00 | 8.55459555e+00 |
| y | 8.10379124e+00 | 7.97201327e+00 | 8.10368439e+00 | 8.10379124e+00 |

| plane | c₂ of f4 − ½[f3,F3] | c₂ of f4 + ½[f3,F3] | c₂ in f4 alone | |q² − p²| of the minus quartic |
| --- | --- | --- | --- | --- |
| x | -2.68750545e+01 | -2.65654545e+01 | -2.67202545e+01 | 4.47604967e+00 |
| y | -2.54588110e+01 | -2.50448183e+01 | -2.52518147e+01 | 5.64485767e+00 |

Second chromaticity: f4 − ½[f3, F3] reproduces d²Q/dδ² with the ξ1 factor, on the bare ring and with one sextupole. f4 + ½[f3, F3] is stored as xi2_plus and is not the reported value. Extra factor on top of the ξ1 conversion: 1.

## One sextupole

The slot is the existing sextupole with the largest `|β_x D_x|` at its centre. Centre values are the drift transport of the entrance Twiss across half the element. The horizontal formula is `+k2 L β_x D_x / (4π)`. The vertical formula is `−k2 L β_y D_x / (4π)`. `D_x` is the horizontal dispersion in both planes. The integral column uses `k2 ∫ β D_x ds / (4π)` with the same signs.

| quantity | value |
| --- | --- |
| element index | 59 |
| name | SF1 |
| s at entrance [m] | 8.21118072e+01 |
| L [m] | 1.20000000e-01 |
| k2 [1/m³] | 1.00000000e+00 |
| βx at centre [m] | 1.93675294e+01 |
| βy at centre [m] | 4.05447112e+00 |
| Dx at centre [m] | 1.17016103e+01 |
| formula Δξx centre | 2.16416932e+00 |
| formula Δξy centre | -4.53055310e-01 |
| formula Δξx integral | 2.16417625e+00 |
| formula Δξy integral | -4.53088382e-01 |

Other sextupoles stay at `k2 = 0`. Linear matrices of the bare and excited rings differ by 0.00000000e+00.

| plane | quantity | generator | orders 1e-4 | orders 1e-6 | series |
| --- | --- | --- | --- | --- | --- |
| x | ξ1 | 7.02640529e-01 | 7.02649751e-01 | 7.02640530e-01 | 7.02640529e-01 |
| x | ξ2 | 1.22897420e+02 | 1.22898754e+02 | 1.22897470e+02 | 1.22897420e+02 |
| y | ξ1 | -2.67857332e+00 | -2.67857656e+00 | -2.67857332e+00 | -2.67857332e+00 |
| y | ξ2 | 7.34244467e+00 | 7.34284980e+00 | 7.34229344e+00 | 7.34244467e+00 |

The same bracket test on the sextupole. The plus combination changes sign relative to `chromaticity_series`. The minus combination stays with it.

| plane | ξ2 from f4 − ½[f3,F3] | ξ2 from f4 + ½[f3,F3] | orders 1e-6 | series |
| --- | --- | --- | --- | --- |
| x | 1.22897420e+02 | -1.20506701e+02 | 1.22897470e+02 | 1.22897420e+02 |
| y | 7.34244467e+00 | -6.65847996e+00 | 7.34229344e+00 | 7.34244467e+00 |

| plane | Δξ1 generator | Δξ1 orders 1e-6 | Δξ1 series | formula centre | formula integral |
| --- | --- | --- | --- | --- | --- |
| x | 2.16417625e+00 | 2.16417625e+00 | 2.16417625e+00 | 2.16416932e+00 | 2.16417625e+00 |
| y | -4.53088382e-01 | -4.53088382e-01 | -4.53088382e-01 | -4.53055310e-01 | -4.53088382e-01 |

Sextupole survey, centre values, bare strengths:

| index | name | s0 [m] | L [m] | βx [m] | βy [m] | Dx [m] | role |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 9 | SF2 | 1.16788296e+01 | 1.20000000e-01 | 1.64584225e+01 | 3.32369705e+00 | 5.85183022e+00 |  |
| 13 | SD | 1.75482444e+01 | 1.20000000e-01 | 3.64896769e+00 | 1.87680353e+01 | 4.89102216e+00 |  |
| 17 | SF1 | 2.34176592e+01 | 1.20000000e-01 | 1.93675294e+01 | 4.05447112e+00 | 1.17016103e+01 |  |
| 21 | SD | 2.92870740e+01 | 1.20000000e-01 | 3.64896769e+00 | 1.87680353e+01 | 4.89102216e+00 |  |
| 25 | SF2 | 3.51564888e+01 | 1.20000000e-01 | 1.64584225e+01 | 3.32369705e+00 | 5.85183022e+00 |  |
| 51 | SF2 | 7.03729776e+01 | 1.20000000e-01 | 1.64584225e+01 | 3.32369705e+00 | 5.85183022e+00 |  |
| 55 | SD | 7.62423924e+01 | 1.20000000e-01 | 3.64896769e+00 | 1.87680353e+01 | 4.89102216e+00 |  |
| 59 | SF1 | 8.21118072e+01 | 1.20000000e-01 | 1.93675294e+01 | 4.05447112e+00 | 1.17016103e+01 | chosen |
| 63 | SD | 8.79812220e+01 | 1.20000000e-01 | 3.64896769e+00 | 1.87680353e+01 | 4.89102216e+00 |  |
| 67 | SF2 | 9.38506368e+01 | 1.20000000e-01 | 1.64584225e+01 | 3.32369705e+00 | 5.85183022e+00 |  |

## Agreement

On each case, `chromaticity_orders` at step `1e-6` and `chromaticity_series` agree within 0.0001 on ξ1 and within 0.001 on ξ2. `ring_matrix` and the Lie factor `R` agree within 1.0e-8, and so do `periodic_optics` dispersion and `dispersion_orders`. The finite-difference steps are printed so the truncation of the default step `1e-4` can be seen.

A generator value is this study's reading of the Lie factor. Where it disagrees with a BELL column, both numbers stay in the table.

## What symplectic_eigen_A does on this factor

Called on `F.R` of the bare ring and not used afterwards.

```
Elliptic pairs with Im λ > 1e-8: 2.
Eigenvalues: -9.712144e-01-2.382071e-01im, -9.712144e-01+2.382071e-01im, -9.706386e-01-2.405425e-01im, -9.706386e-01+2.405425e-01im, +1.000000e+00+0.000000e+00im, +1.000000e+00+0.000000e+00im.
symplectic_eigen_A: threw: one-turn matrix not stable
```
