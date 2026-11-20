# Dipole fringe

Lattice file: `C:\GIT_REPOS\SCT_2\studies\magnetic_ring\lattices\magnetic_2_closed.bell`. Reference: deuteron, kinetic energy 270.0 MeV, harmonic 1. Placement is `place_sextupole_families` on `optics_along` at ds = 0.02 m of the closed ring, then `set_family_strengths` with `k2_per_m3.own` from `output/05_sext_strengths.json`. Strengths are not solved again.

Own strengths k2 [1/m³]: SXD = 8.17687255e-01, SXB = -4.24654677e-01, SYB = -6.11237385e-01. Family indices match step 5.

Each row is `ring_lie_factor(lat, model)` with the fringe model on dipoles. Quadrupoles and sextupoles stay `NoEdge`. `Qx` and `Qy` are the body tunes from `optics_along` plus the shortest fringe shift of μ/2π from `periodic_optics` of that factor. `ξ` is `generator_chromaticity` with the step-2b factor `ξ = −c / (2π)`. The normalising matrix is built from `periodic_optics` of the same factor. `η1` is `path_expansion` of `factor_to_map` of that factor.

## Assumption

The soft-edge width is an assumption of this step: `LinearRamp(0.05)`, width 0.05 m. The closed lattice does not carry this width. BELL's half-width of that profile is W = 2.50000000e-02 m (profile shift 8.67361738e-19 m). The fringe integral ∫ F(1−F) dz = 8.33333333e-03 m. The default soft-edge method is `Splitting(n=20, order=4)`.

Loaded BM card: 16 magnets, L = 5.32941480e+00 m, angle = 3.92699082e-01 rad, e1 = 0.00000000e+00, e2 = 0.00000000e+00, gap = 0.00000000e+00 m, FINT = 0.00000000e+00, k1 = 0.00000000e+00. The ramp satisfies 2W < L.

Body tunes, with no edge model: Qx = 2.46172015e+00, Qy = 2.46133732e+00.

## Comparison

Regime file `output/02_rf_regimes.json`, object `linear`. Voltage 0.00678732526132278 MV, harmonic 1, δ_m = 0.0001, Q_s = 0.001067377255883277, φ̂ = 0.05000000000000002. The notes' stationary phase in that file is 0.0. The cavity in the track is at BELL phase π. Actions are the placeholder Jx = 2.581612713258228e-8, Jy = 2.5820256644667e-8, the same values as steps 4 and 5. Λ = −(2π/C) (Jx ξx + Jy ξy) uses the ξ of that fringe model. The linear-regime formula is δ_eq = −(η1 δ_m²/2 + Λ) / η0. η0 is the slip of that same factor. γ_eff = γ + β² γ δ_eq + β² ⟨δ²⟩ / (2γ), with ⟨δ²⟩ = δ_m²/2 in the formula and the harmonic mean square in the track.

| model | Qx | Qy | ξx | ξy | η0 | η1 | Λ | formula δ_eq | track δ_eq | formula γ_eff | track γ_eff |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NoEdge | 2.46172015e+00 | 2.46133732e+00 | 3.00385720e-16 | 8.01764826e-16 | -5.33688628e-01 | -7.77156117e-16 | -1.52313225e-24 | -1.01349599e-23 | -7.13955154e-13 | 1.14395294e+00 | 1.14395294e+00 |
| HardEdge() | 2.46172015e+00 | 2.46133732e+00 | 3.00385720e-16 | 4.08421638e-01 | -5.33688628e-01 | -7.77156117e-16 | -5.64448558e-10 | -1.05763647e-09 | -1.05835152e-09 | 1.14395294e+00 | 1.14395294e+00 |
| SoftEdge(profile=LinearRamp(0.05)) | 2.46171907e+00 | 2.46006817e+00 | 3.83856209e-04 | 4.11149163e-01 | -5.33687726e-01 | -1.03809949e-04 | -5.68748476e-10 | -1.06666783e-09 | -1.06729934e-09 | 1.14395294e+00 | 1.14395294e+00 |

`SoftEdge` uses the assumed width 0.05 m. `HardEdge()` is `HardEdge(order=2)` on dipoles, with the magnet gap and FINT (both zero on this card).

## Each model

### NoEdge

Dipole model: `NoEdge`. Fractional tunes μ/2π: Qx = 4.61720150e-01, Qy = 4.61337323e-01. Fringe shift: ΔQx = 1.34290087e-15, ΔQy = 4.80617152e-15. |R − body matrix| = 2.70006240e-13. Symplectic error of R = 1.67782455e-14.

The linear map matches the body within 1.0e-8. The normalising matrix was still built from this factor's `periodic_optics`. βx = 1.93675224e+01 m, βy = 4.05447047e+00 m, Dx = 1.10941104e-02 m.

Generator chromaticity ξx = 3.00385720e-16, ξy = 8.01764826e-16. `chromaticity_orders` reads ξx = -3.79909437e-07, ξy = 2.06814288e-07. Gaps: 3.79909438e-07, -2.06814287e-07. Rotation residual of the normalised linear map: 1.87627691e-14. Normalisation within the step-2b tolerance: true.

Slip: η0 = -5.33688628e-01, η1 = -7.77156117e-16, α0 = 2.30470282e-01, α1 = -6.78153026e-01. Identity η1 = α1 − α0/γ² + (3γ² − 1)/(2γ⁴) gives -7.77156117e-16.

Track: 0 edge items, 937 turns, ds = 0.1 m, seed Q_s = 0.001067377255883277, note Q_s from this η0 = 0.001067377255883277. Fit Q_s = 1.06721221e-03, rms = 3.89338930e-13, half-window split = 1.81501151e-14. Window left early: false.

### HardEdge

Dipole model: `HardEdge()`. Fractional tunes μ/2π: Qx = 4.61720150e-01, Qy = 4.61337323e-01. Fringe shift: ΔQx = 1.34290087e-15, ΔQy = 4.80617152e-15. |R − body matrix| = 2.70006240e-13. Symplectic error of R = 1.67782455e-14.

The linear map matches the body within 1.0e-8. The normalising matrix was still built from this factor's `periodic_optics`. βx = 1.93675224e+01 m, βy = 4.05447047e+00 m, Dx = 1.10941104e-02 m.

Generator chromaticity ξx = 3.00385720e-16, ξy = 4.08421638e-01. `chromaticity_orders` reads ξx = -3.79909437e-07, ξy = 4.08420630e-01. Gaps: 3.79909438e-07, 1.00818760e-06. Rotation residual of the normalised linear map: 1.87627691e-14. Normalisation within the step-2b tolerance: true.

Slip: η0 = -5.33688628e-01, η1 = -7.77156117e-16, α0 = 2.30470282e-01, α1 = -6.78153026e-01. Identity η1 = α1 − α0/γ² + (3γ² − 1)/(2γ⁴) gives -7.77156117e-16.

Track: 32 edge items, 937 turns, ds = 0.1 m, seed Q_s = 0.001067377255883277, note Q_s from this η0 = 0.001067377255883277. Fit Q_s = 1.06721221e-03, rms = 3.89505177e-13, half-window split = 6.63318214e-14. Window left early: false.

### SoftEdge

Dipole model: `SoftEdge(profile=LinearRamp(0.05))`. Fractional tunes μ/2π: Qx = 4.61719074e-01, Qy = 4.60068172e-01. Fringe shift: ΔQx = -1.07612720e-06, ΔQy = -1.26915133e-03. |R − body matrix| = 3.14000626e-02. Symplectic error of R = 1.19074889e-13.

The linear optics differ from the body, so the normalising matrix was built from this factor: βx = 1.93675231e+01 m, βy = 4.05466819e+00 m, Dx = 1.10820459e-02 m.

Generator chromaticity ξx = 3.83856209e-04, ξy = 4.11149163e-01. `chromaticity_orders` reads ξx = 3.83477369e-04, ξy = 4.11148145e-01. Gaps: 3.78840200e-07, 1.01805676e-06. Rotation residual of the normalised linear map: 2.08499884e-13. Normalisation within the step-2b tolerance: true.

Slip: η0 = -5.33687726e-01, η1 = -1.03809949e-04, α0 = 2.30471184e-01, α1 = -6.78256146e-01. Identity η1 = α1 − α0/γ² + (3γ² − 1)/(2γ⁴) gives -1.03809949e-04.

Track: 32 edge items, 937 turns, ds = 0.1 m, seed Q_s = 0.001067377255883277, note Q_s from this η0 = 0.0010673763537413106. Fit Q_s = 1.06721130e-03, rms = 3.84239380e-13, half-window split = 6.90410032e-14. Window left early: false.

## Tracking

A turn follows `fringe_sequence`. An `EdgeMap` is `apply_factor` of `dragt_finn_factor(edge_map(...))`, the same Lie factor the one-turn product uses, truncated at map order 3. A magnet is `exact_track_element` at ds = 0.1 m. The thin cavity is `rf_kick` at phase π. `exact_track_element` is not called on an `EdgeMap`.

Probe: `exact_track_element` on `BM_in` (BM, side in) raises `MethodError`. The message is recorded in `output/06_bell_finding.md`.

No stop. Formula and fringe-track δ_eq agree within the step-4 tolerance, generator chromaticity agrees with `chromaticity_orders` within 0.001, and η1 agrees with the slip identity.
