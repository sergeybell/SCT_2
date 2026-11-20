# Lie invariants of the closed magnetic ring

Lattice file: `C:\GIT_REPOS\SCT_2\studies\magnetic_ring\lattices\magnetic_2_closed.bell`. Sextupoles are at `k2 = 0`. Reference: deuteron, kinetic energy 270.0 MeV, gamma 1.14395294e+00, harmonic 1.

The magnetic factor is `ring_lie_factor` of `magnetic_lattice`. Slip coefficients are `path_expansion` of `factor_to_map`. Chromaticity in the generator column is `generator_chromaticity` with the factor fixed in step 2b, `ξ1 = −c / (2π)` and `ξ2 = 2 ξ1_factor c₂`, placed next to `chromaticity_orders`.

## Slip

The note formula is `η1 = α1 − α0/γ² + (3γ² − 1)/(2γ⁴)`. `eta1 − eta1_formula` uses that expression. `path_expansion` also returns its own `eta1_formula`.

| quantity | value |
| --- | --- |
| C [m] | 1.17388296e+02 |
| gamma | 1.14395294e+00 |
| 1/gamma^2 | 7.64158910e-01 |
| alpha0 | 2.30470282e-01 |
| alpha1 | 3.39007701e-01 |
| eta0 | -5.33688628e-01 |
| eta1 | 1.01716073e+00 |
| eta1_formula | 1.01716073e+00 |
| eta1 − eta1_formula | 0.000000e+00 |
| path_expansion eta1_formula | 1.01716073e+00 |
| eta1 − path_expansion eta1_formula | 0.000000e+00 |
| note formula − path_expansion formula | 0.000000e+00 |
| (-c2 + c1^2) − (3γ²−1)/(2γ⁴) | 0.000000e+00 |
| |Lie R − ring_matrix| | 2.700062e-13 |

η1 and the note formula agree within 1.0e-8. The note formula and `path_expansion`'s `eta1_formula` agree within the same tolerance. The velocity identity `−c2 + c1² = (3γ² − 1)/(2γ⁴)` is the algebraic reason those two formulas meet.

## Chromaticity

Generator ξ1 is `−1/(2π)` times the coefficient of `J δ`. Generator ξ2 is `−1/π` times the coefficient of `J δ²` in `f4 − ½[f3, F3]`. `chromaticity_orders` is the central difference of the same magnetic map. The step `1e-6` is the settled difference. The step `1e-4` is the default of that function.

| plane | quantity | generator | orders 1e-4 | orders 1e-6 | generator − orders 1e-6 |
| --- | --- | --- | --- | --- | --- |
| x | xi1 | -1.46153572e+00 | -1.46153752e+00 | -1.46153572e+00 | 1.67234226e-10 |
| x | xi2 | 8.55459555e+00 | 8.55473511e+00 | 8.55454596e+00 | 4.95857707e-05 |
| y | xi1 | -2.22548493e+00 | -2.22548794e+00 | -2.22548493e+00 | 2.90293123e-10 |
| y | xi2 | 8.10379124e+00 | 8.10412260e+00 | 8.10368439e+00 | 1.06853598e-04 |

| plane | c of J δ | c of J δ² |
| --- | --- | --- |
| x | 9.18309974e+00 | -2.68750545e+01 |
| y | 1.39831342e+01 | -2.54588110e+01 |

Generator and `chromaticity_orders` at step `1e-6` agree within 0.0001 on ξ1 and within 0.001 on ξ2.

## Thin cavity

One thin `RFCavity` is appended to the element list. Voltage and harmonic are `regimes.linear` in `02_rf_regimes.json`. The cavity phase passed to BELL is `3.14159265e+00` rad. The regime file records the notes' stationary phase `0.00000000e+00` rad. `eta0` is negative, and `eta0 * cos(phi_s)` is positive at the BELL phase, so the smooth tune is real.

The one-turn factor used below is `ring_lie_factor` of this lattice. The cavity is already an element, so its kick is inside that product. `one_turn_factor_rf` is not called.

| quantity | value |
| --- | --- |
| h | 1 |
| V [MV] | 6.78732526e-03 |
| phi_s on the cavity [rad] | 3.14159265e+00 |
| phi_s in the regime file [rad] | 0.00000000e+00 |
| Q_s note | 1.06737726e-03 |
| Q_s smooth model at phi_s = π | 1.06737726e-03 |
| Q_s stored in the regime file | 1.06737726e-03 |
| Q_s note − regime file | 0.000000e+00 |
| Q_s note − smooth model | 0.000000e+00 |

## Two ways to attach one kick, and a second kick

`ring_lie_factor` of the lattice with the cavity is compared with `concatenate` of the magnetic factor and `dragt_finn_factor(rf_kick_taylor(...))`. A further `concatenate` of the cavity-lattice factor with that same kick is the map `one_turn_factor_rf` builds when the lattice already contains the cavity: `ring_lie_factor` includes the cavity, and the function appends `rf_kick_taylor` again. Synchrotron tunes in this section are `|μ_s| / 2π` from `symplectic_eigen_A` of those linear matrices.

| quantity | value |
| --- | --- |
| |R(cavity lattice) − R(magnetic + one kick)| | 2.842171e-14 |
| |f3| gap of that comparison | 3.907985e-14 |
| |f4| gap of that comparison | 6.661338e-14 |
| Q_s of the cavity lattice | 1.06737926e-03 |
| Q_s after a second kick | 1.50950505e-03 |
| Q_s second / Q_s one | 1.41421621e+00 |

## Normal form

`dragt_finn_nf6d` returned. `Q_s = |μ_s| / 2π`. The detuning matrix is stored as `h4 = Σ_{i≤j} c_ij J_i J_j`, with zeros below the diagonal. `μ_k(J) = μ_k − ∂h4/∂J_k`, and `∂h4/∂J_k = 2 c_kk J_k + Σ_{j≠k} c_{min,max} J_j`.

| quantity | value |
| --- | --- |
| mu_x [rad] | 2.90107326e+00 |
| mu_y [rad] | 2.89866789e+00 |
| mu_s [rad] | 6.70654166e-03 |
| Q_s | 1.06737926e-03 |
| Q_s − Q_s note | 2.000339e-09 |
| wrapped |Q_x − Q_x magnetic| | 6.452061e-13 |
| wrapped |Q_y − Q_y magnetic| | 1.110223e-15 |
| |A^-1 R A − rotation| | 4.519475e-14 |
| |A^T S A − S| | 3.591450e-15 |
| max |coef| of f3 | 7.366624e-02 |
| cubic residual after F3 | 1.776357e-15 |
| |f4 − h4| after F4 | 6.221412e-14 |
| sigma_delta per sqrt(2 Js) | 1.03465069e-02 |
| sigma_l per sqrt(2 Js) | 9.66515199e+01 |

Detuning matrix `c`, rows and columns `(x, y, s)`:

| c_i1 | c_i2 | c_i3 |
| --- | --- | --- |
| -2.32364474e+00 | -3.81134833e+00 | -4.71742246e-03 |
| 0.00000000e+00 | -3.48047667e+00 | -5.53396969e-03 |
| 0.00000000e+00 | 0.00000000e+00 | 1.12162519e-02 |

| plane | d_k |
| --- | --- |
| x | 1.46580620e-01 |
| y | 2.23198761e-01 |
| s | 2.04027704e-04 |

`d` is the normal-form coefficient in `δ_eq = z*_6 + d · J`. `z*` from this call is the origin.

## Tracking

Each magnet is `exact_track_element` with step 0.1 m. The thin cavity is `rf_kick`. The observation point is the start of the element list, which is the entrance of the one-turn map. The launch point is the periodic dispersion orbit at `δ = δ_m` plus a betatron offset. The horizontal offset is twice `|D_x δ_m|` at this location (`a_x = 2.21882208e-06` m, `a_y = 2.21882208e-06` m, `δ_m = 1.00000000e-04`). `l = 0`, so the particle starts at the momentum crest of the stationary bucket.

Normal-form actions use `u = A⁻¹ (z − z*)` and `w = Tinv(u)`, then `J_k = (w_{2k−1}² + w_{2k}²) / 2`. Courant–Snyder `J` uses the on-momentum Twiss at this location. The geometric Courant–Snyder action keeps the dispersive piece of `(x, p_x)`. The dispersion-free column subtracts `D δ` first.

1874 turns, including the launch point as turn 0. Elapsed time 46.9 s. Synchrotron periods at the normal-form Q_s: 2.000.

| action | mean | min | max | max − min | (max−min)/|mean| |
| --- | --- | --- | --- | --- | --- |
| Jx NF | 1.39216466e-13 | 1.38353798e-13 | 1.40069956e-13 | 1.71615771e-15 | 1.232726e-02 |
| Jy NF | 6.07161030e-13 | 6.07046053e-13 | 6.07341523e-13 | 2.95470436e-16 | 4.866426e-04 |
| Js NF | 4.67076657e-05 | 4.67076564e-05 | 4.67076751e-05 | 1.87153965e-11 | 4.006922e-07 |
| Jx Courant-Snyder | 1.67708748e-13 | 1.32330336e-14 | 3.99169943e-13 | 3.85936909e-13 | 2.301233e+00 |
| Jy Courant-Snyder | 6.07053494e-13 | 6.06972067e-13 | 6.07134915e-13 | 1.62847309e-16 | 2.682586e-04 |
| Jx Courant-Snyder, dispersion removed | 1.39379666e-13 | 1.23468653e-13 | 1.55794627e-13 | 3.23259738e-14 | 2.319275e-01 |
| Jy Courant-Snyder, dispersion removed | 6.07053494e-13 | 6.06972067e-13 | 6.07134915e-13 | 1.62847309e-16 | 2.682586e-04 |

Figure: `03_actions.svg`. Samples: `03_actions.csv`.
