# Sextupole families

Lattice file: `C:\GIT_REPOS\SCT_2\studies\magnetic_ring\lattices\magnetic_2_closed.bell`. Reference: deuteron, kinetic energy 270.0 MeV, harmonic 1. Placement uses `optics_along` at ds = 0.02 m and `place_sextupole_families`. Sextupoles in the closed file are at `k2 = 0`.

Chromaticity in the in-study matrix is `generator_chromaticity` with the step-2b factor `ξ = −c / (2π)`. `η1` and `α1` are `path_expansion` of `factor_to_map` of `ring_lie_factor` on `magnetic_lattice`. The neighbour cancels `ξx`, `ξy`, and `α1` with `response` and `solve_A`. The finite-difference step is dk = 1.0 for both matrices.

## Slots

Each row is one magnet. `beta_x`, `beta_y`, and `D` are the values at the optics maximum that claimed the slot. The slot is the nearest sextupole or short drift.

| family | quantity | s_max [m] | slot | s_slot [m] | beta_x [m] | beta_y [m] | D [m] | nearest quad |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SXD | D_x | 2.34176592e+01 | SF1 | 2.34776592e+01 | 1.93677152e+01 | 4.05535903e+00 | 1.17016103e+01 | QF (F) |
| SXD | D_x | 8.22318072e+01 | SF1 | 8.21718072e+01 | 1.93677152e+01 | 4.05535903e+00 | 1.17016103e+01 | QF (F) |
| SXB | beta_x | 1.16788296e+01 | SF2 | 1.17388296e+01 | 1.64586564e+01 | 3.32479147e+00 | 5.83218650e+00 | QF (F) |
| SXB | beta_x | 3.51564888e+01 | SF2 | 3.52164888e+01 | 1.64586262e+01 | 3.32476890e+00 | 5.87147395e+00 | QF (F) |
| SXB | beta_x | 7.03729776e+01 | SF2 | 7.04329776e+01 | 1.64586564e+01 | 3.32479147e+00 | 5.83218650e+00 | QF (F) |
| SXB | beta_x | 9.38506368e+01 | SF2 | 9.39106368e+01 | 1.64586262e+01 | 3.32476890e+00 | 5.87147395e+00 | QF (F) |
| SYB | beta_y | 5.80941480e+00 | OS1 | 5.86941480e+00 | 3.66065321e+00 | 1.87836872e+01 | 1.12776958e+00 | QD (D) |
| SYB | beta_y | 1.76682444e+01 | SD | 1.76082444e+01 | 3.65973724e+00 | 1.87801732e+01 | 4.92180878e+00 | QD (D) |
| SYB | beta_y | 2.92870740e+01 | SD | 2.93470740e+01 | 3.65973724e+00 | 1.87801732e+01 | 4.92180878e+00 | QD (D) |
| SYB | beta_y | 4.11459036e+01 | OS1 | 4.10859036e+01 | 3.66065321e+00 | 1.87836872e+01 | 1.12776958e+00 | QD (D) |
| SYB | beta_y | 6.45035628e+01 | OS1 | 6.45635628e+01 | 3.66065321e+00 | 1.87836872e+01 | 1.12776958e+00 | QD (D) |
| SYB | beta_y | 7.63623924e+01 | SD | 7.63023924e+01 | 3.65973724e+00 | 1.87801732e+01 | 4.92180878e+00 | QD (D) |
| SYB | beta_y | 8.79812220e+01 | SD | 8.80412220e+01 | 3.65973724e+00 | 1.87801732e+01 | 4.92180878e+00 | QD (D) |
| SYB | beta_y | 9.98400516e+01 | OS1 | 9.97800516e+01 | 3.66065321e+00 | 1.87836872e+01 | 1.12776958e+00 | QD (D) |

Notes from `place_sextupole_families`:

- beta_x max at s=23.418 m shares slot SF1 with family SXD
- beta_x max at s=82.112 m shares slot SF1 with family SXD

Family SXD: 2 magnets, L = 0.12 m.
Family SXB: 4 magnets, L = 0.12 m.
Family SYB: 8 magnets, L = 0.12 m.

## Picture

The placement picture is SXD at a dispersion maximum, SXB near a horizontal beta maximum, and SYB near a vertical beta maximum.

SXD has 2 magnets, at the kept dispersion maximum (placed 11.7016, highest kept 11.7016). Nearest quadrupoles are focusing (QF).

SXD keeps the slot that a horizontal beta peak also asked for.

SXB has 4 magnets, below the highest kept horizontal beta maximum (placed 16.4587, highest kept 19.3677). Nearest quadrupoles are focusing (QF). The higher horizontal peaks share their slot with SXD, so this family takes the remaining horizontal peaks.

SYB has 8 magnets, at the kept vertical beta maximum (placed 18.7837, highest kept 18.7837). Nearest quadrupoles are defocusing (QD).

The largest `beta_x` sample is at s = 5.87541480e+01 m, beta_x = 1.93677385e+01 m, D = 1.10637404e-02 m. `find_optics_maxima` keeps a horizontal peak within 0.3 m of it: false. The largest `beta_y` sample is at s = 1.11458881e+02 m, D = 6.77386256e-03 m, kept: false. The largest dispersion sample is at s = 8.22318072e+01 m, D = 1.17016103e+01 m.

`find_optics_maxima` keeps a local maximum only where D > 0.5 m, including the beta curves. A horizontal peak that falls in a slot already given to SXD is omitted and printed as a shared-slot note.

## Own response

Columns are family strengths in 1/m³. Rows are dξx/dk2, dξy/dk2, and dη1/dk2. Bare readings at k2 = 0: ξx = -1.46153572e+00, ξy = -2.22548493e+00, η1 = 1.01716073e+00, α1 = 3.39007701e-01, α0 = 2.30470282e-01.

| row | SXD | SXB | SYB |
| --- | --- | --- | --- |
| dξx/dk2 | 4.32835249e+00 | 3.67885967e+00 | 8.43302449e-01 |
| dξy/dk2 | -9.06176765e-01 | -7.43004929e-01 | -4.33699520e+00 |
| dη1/dk2 | -1.63792247e+00 | -4.09700496e-01 | -2.42403976e-01 |

Condition number cond(A) = 8.84452055e+00. Rank = 3. Square = true.

The α1 that sets η1 to zero is α0/γ² − (3γ² − 1)/(2γ⁴) = -6.78153026e-01. γ = 1.14395294e+00.

Achieved after applying the own solution: α1 = -6.78153026e-01, η1 = -7.77156117e-16, ξx = 3.00385720e-16, ξy = 8.01764826e-16. α1 − target = -8.88178420e-16. Residual max(|ξx|, |ξy|, |η1|) = 8.01764826e-16.

Strengths k2 [1/m³]:

| family | own linear | own used | solve_A |
| --- | --- | --- | --- |
| SXD | 8.17687255e-01 | 8.17687255e-01 | 2.30414107e-01 |
| SXB | -4.24654677e-01 | -4.24654677e-01 | 2.65267228e-01 |
| SYB | -6.11237385e-01 | -6.11237385e-01 | -6.06726257e-01 |

NelderMead did not run. All three families are present and the re-evaluated residual is at most 0.001.

## Neighbour, solve_A

`solve_A` cancels ξx, ξy, and α1 using `ring_metrics` (`chromaticity_orders` and `path_expansion`). It does not target η1.

| family | dξx/dk2 | dξy/dk2 | dα1/dk2 | dη1/dk2 |
| --- | --- | --- | --- | --- |
| SXD | 4.32833403e+00 | -9.06177393e-01 | -1.63792247e+00 | -1.63792247e+00 |
| SXB | 3.67886770e+00 | -7.43005237e-01 | -4.09700496e-01 | -4.09700496e-01 |
| SYB | 8.43304089e-01 | -4.33701180e+00 | -2.42403976e-01 | -2.42403976e-01 |

At the `solve_A` strengths, `ring_metrics` reads ξx = 4.91867769e-06, ξy = -6.83152590e-06, α1 = -2.77555756e-16, η1 = 6.78153026e-01, η0 = -5.33688628e-01. Leftover η1 is the slip that remains when α1 is removed and the kinematic piece (3γ² − 1)/(2γ⁴) − α0/γ² is left. Formula δ_eq from these BELL readings, with Λ from the same ξ and the regime actions, is 6.35345704e-09.

The study reading of that same lattice is in the table below. |η1_study − η1_bell| = 0.00000000e+00. |ξx_study − ξx_bell| = 4.45244219e-09.

## Linear-regime equilibrium

Regime file `output/02_rf_regimes.json`, object `linear`. Voltage 0.00678732526132278 MV, harmonic 1, δ_m = 0.0001, Q_s = 0.001067377255883277, φ̂ = 0.05000000000000002. The notes' stationary phase in that file is 0.0. The cavity in the track is at BELL phase π. Actions are the placeholder Jx = 2.581612713258228e-8, Jy = 2.5820256644667e-8. Λ = −(2π/C) (Jx ξx + Jy ξy) uses the ξ of that setting. The formula is δ_eq = −(η1 δ_m²/2 + Λ) / η0. γ_eff = γ + β² γ δ_eq + β² ⟨δ²⟩ / (2γ), with ⟨δ²⟩ = δ_m²/2 in the formula and the harmonic mean square in the track.

The exact column is `exact_track_element` on the magnets and `rf_kick` on one thin cavity, ds = 0.1 m, 937 turns. δ_eq is the constant term of a harmonic fit. The tune starts at the regime Q_s and is refined on the residual, so a fraction of a turn of phase error does not leak into the constant.

| setting | ξx | ξy | η1 | α1 | Λ | formula δ_eq | exact δ_eq | formula γ_eff | exact γ_eff |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| k2=0 | -1.46153572e+00 | -2.22548493e+00 | 1.01716073e+00 | 3.39007701e-01 | 5.09522845e-09 | 1.90767267e-08 | 1.90710469e-08 | 1.14395295e+00 | 1.14395295e+00 |
| own | 3.00385720e-16 | 8.01764826e-16 | -7.77156117e-16 | -6.78153026e-01 | -1.52313225e-24 | -1.01349599e-23 | -7.13955154e-13 | 1.14395294e+00 | 1.14395294e+00 |
| solve_A | 4.92313013e-06 | -6.83486651e-06 | 6.78153026e-01 | -2.77555756e-16 | 2.64315396e-15 | 6.35345704e-09 | 6.35055992e-09 | 1.14395294e+00 | 1.14395294e+00 |

No stop. Formula and exact δ_eq agree within the step-4 tolerance, and η1 agrees with the slip identity.

Strengths for the later steps are in `output/05_sext_strengths.json`. The primary vector is `own`: ξx = ξy = η1 = 0. Rebuild by calling `place_sextupole_families` on `optics_along` of the closed lattice at ds = 0.02 m, then `set_family_strengths`.
