# Stationary RF regimes for the closed magnetic ring

Chosen after `eta0`, `eta1`, and the tunes of this step. The lattice harmonic is kept. `V = 0.2` MV and `h = 5` are not reused.

The bucket is stationary: `sin phi_s = 0` and `cos phi_s = -sign(eta0)`, so a positive voltage is stable. `phi_s = 0.0` rad.

Reference energy 2145.612943 MeV, beta 0.48563473, gamma 1.14395294, eta0 -5.33688628e-01, eta1 1.01716073e+00.
Revolution period T_s = 8.06295728e-07 s.

Selection cuts, used as a concrete stand-in for much-less-than-one:

| cut | value |
| --- | --- |
| abs(eta1*delta_m/eta0) | 2.00000000e-02 |
| abs(Lambda/(eta0*delta_m)) | 2.00000000e-02 |
| delta_m | 1.00000000e-02 |
| omega_s T_s | 2.00000000e-02 |
| omega_s T_s / betatron phase | 5.00000000e-02 |
| linear phi_hat target | 5.00000000e-02 |
| cubic phi_hat target | 3.00000000e-01 |

Shared harmonic h = 1, shared delta_m = 1.00000000e-04.

Placeholder action, not a measured amplitude. Peak displacement `a` at the beta maximum, `J = a^2 / (2 beta_max)`. Lambda uses the natural first chromaticity of the same Lie polynomial (`xi_x = -1.46153752e+00`, `xi_y = -2.22548794e+00`). That chromaticity is only a size for Lambda. The chromaticity table is a later step.

| quantity | value |
| --- | --- |
| amplitude_m | 1.00000000e-03 |
| beta_x_max_m | 1.93677385e+01 |
| beta_y_max_m | 1.93646410e+01 |
| J_x_m | 2.58161271e-08 |
| J_y_m | 2.58202566e-08 |
| Lambda | 5.09523510e-09 |

## A. Linear regime

Target phase half-range 0.05, so `sin phi ≈ phi` applies.

| quantity | value | role |
| --- | --- | --- |
| abs(eta1*delta_m/eta0) | 1.90590669e-04 | much less than 1 |
| abs(Lambda/(eta0*delta_m)) | 9.54720568e-05 | much less than 1 |
| phi_hat | 5.00000000e-02 | much less than 1 |
| phi_hat^2 | 2.50000000e-03 | follows from phi_hat |
| delta_m | 1.00000000e-04 | much less than 1 |
| omega_s T_s | 6.70652909e-03 | much less than 1 |
| 2 pi {nu_x} | 2.90107326e+00 | printed |
| 2 pi {nu_y} | 2.89866789e+00 | printed |
| omega_s T_s / (2 pi {nu_x}) | 2.31174069e-03 | much less than 1 |
| omega_s T_s / (2 pi {nu_y}) | 2.31365901e-03 | much less than 1 |
| 2 pi dist(nu_x, integer) | 2.90107326e+00 | phase to the nearest integer |
| 2 pi dist(nu_y, integer) | 2.89866789e+00 | phase to the nearest integer |
| omega_s T_s / that x distance | 2.31174069e-03 | much less than 1 |
| omega_s T_s / that y distance | 2.31365901e-03 | much less than 1 |
| Q_s | 1.06737726e-03 | omega_s T_s / 2 pi |
| V_MV | 6.78732526e-03 | stationary bucket |
| formula delta_eq | 1.90767391e-08 | notes, with placeholder Lambda |

Smooth-model Q_s at its own stable phase (3.141593 rad) is 1.06737726e-03. The absolute difference from the note tune was required to be at most 1e-9.

## B. Cubic-phase regime

Target phase half-range 0.30. `phi_hat` is not small. `phi_hat^2 = 0.09` is the condition that replaces `phi_hat ≪ 1`. The other conditions stay.

| quantity | value | role |
| --- | --- | --- |
| abs(eta1*delta_m/eta0) | 1.90590669e-04 | much less than 1 |
| abs(Lambda/(eta0*delta_m)) | 9.54720568e-05 | much less than 1 |
| phi_hat | 3.00000000e-01 | not the linear condition |
| phi_hat^2 | 9.00000000e-02 | much less than 1 |
| delta_m | 1.00000000e-04 | much less than 1 |
| omega_s T_s | 1.11775485e-03 | much less than 1 |
| 2 pi {nu_x} | 2.90107326e+00 | printed |
| 2 pi {nu_y} | 2.89866789e+00 | printed |
| omega_s T_s / (2 pi {nu_x}) | 3.85290114e-04 | much less than 1 |
| omega_s T_s / (2 pi {nu_y}) | 3.85609836e-04 | much less than 1 |
| 2 pi dist(nu_x, integer) | 2.90107326e+00 | phase to the nearest integer |
| 2 pi dist(nu_y, integer) | 2.89866789e+00 | phase to the nearest integer |
| omega_s T_s / that x distance | 3.85290114e-04 | much less than 1 |
| omega_s T_s / that y distance | 3.85609836e-04 | much less than 1 |
| Q_s | 1.77896209e-04 | omega_s T_s / 2 pi |
| V_MV | 1.88536813e-04 | stationary bucket |
| formula delta_eq | 1.90499373e-08 | notes, with placeholder Lambda |

Smooth-model Q_s at its own stable phase (3.141593 rad) is 1.77896209e-04. The absolute difference from the note tune was required to be at most 1e-9.

`formula delta_eq` evaluates `-(eta1 * delta_m^2 / 2 * f + Lambda) / eta0`, with `f = 1` in the linear regime and `f = 1 - phi_hat^2/32` in the cubic regime. It is not a tracking result.
