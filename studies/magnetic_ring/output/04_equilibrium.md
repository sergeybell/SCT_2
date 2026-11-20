# Equilibrium momentum of the closed magnetic ring

## STOP

- C cubic J=0 map δ_eq is closer to chapter 1-2 than to 1 - phi_hat^2/32. measured 9.544986450512834e-9, chapter 1-2 9.529533453586332e-9, chapter 3 9.50273164074812e-9. The tune and the third harmonic of this series do carry phi_hat^2; the mean momentum does not.

Sextupoles are at `k2 = 0`. One thin cavity is appended. BELL stores `phi_s = π` because `l > 0` is an early arrival. The notes' phase is late-positive, so `φ^N = −φ^B`. Below transition the notes' stationary phase is 0, and `η0 cos φ_s^B > 0` at `φ_s^B = π`.

The normal-form torus is the average of `U6` on the action torus. It is the quadratic piece `z*_6 + d · J` (`d_grid` is the same average at action `1e-8`). That polynomial stops before `φ̂²`, so the torus is compared with chapters 1–2 even in the cubic-phase regime.

`δ_eq`, `δ_m`, and the second harmonic of a tracked series come from a harmonic least-squares fit, harmonics 1–3, with the tune refined on a grid. `δ_m` in the table is the half-range `(max − min) / 2`. The signed second harmonic is the coefficient of `cos 2ψ` after `ψ` is measured from the crest of the fundamental. `⟨δ²⟩` is `(a1² + a2² + a3²) / 2`, the average of that waveform over a full period. `γ_eff` and `ΔE/E` use that column's `δ_eq` and `⟨δ²⟩`. Agreement is judged on `δ_eq`.

## Phase dictionary

Linear-regime voltage, `J = 0`, launch at `δ = δ_m` and `l = 0`. The control uses the same voltage with `phi_s = 0`. A run that reaches `|δ| = 0.01` stops early.

| phi_s BELL | turns asked | samples | aborted | min delta | max delta | max |delta| | oscillates | left bucket |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3.14159265e+00 | 937 | 938 | false | -9.99929508e-05 | 1.00000046e-04 | 1.00000046e-04 | true | false |
| 0.00000000e+00 | 937 | 938 | false | 1.00000000e-04 | 4.00735480e-03 | 4.00735480e-03 | false | true |

## Slip, chromaticity, inequalities

| quantity | value |
| --- | --- |
| C_m | 1.17388296e+02 |
| gamma | 1.14395294e+00 |
| beta | 4.85634729e-01 |
| eta0 | -5.33688628e-01 |
| eta1 | 1.01716073e+00 |
| xi_x_generator | -1.46153572e+00 |
| xi_y_generator | -2.22548493e+00 |
| Jx | 2.58161271e-08 |
| Jy | 2.58202566e-08 |
| Lambda | 5.09522845e-09 |
| Lambda_json_placeholder | 5.09523510e-09 |
| abs_eta1_delta_m_over_eta0 | 1.90590669e-04 |
| delta_m | 1.00000000e-04 |
| abs_Lambda_over_eta0_delta_m | 9.54719321e-05 |
| linear_phi_hat | 5.00000000e-02 |
| linear_phi_hat_squared | 2.50000000e-03 |
| linear_omega_s_Ts | 6.70652909e-03 |
| linear_omega_over_betatron_x | 2.31174069e-03 |
| linear_omega_over_betatron_y | 2.31365901e-03 |
| linear_Qs | 1.06737726e-03 |
| linear_Qs_smooth_phi_pi | 1.06737726e-03 |
| linear_voltage_MV | 6.78732526e-03 |
| cubic_phase_phi_hat | 3.00000000e-01 |
| cubic_phase_phi_hat_squared | 9.00000000e-02 |
| cubic_phase_omega_s_Ts | 1.11775485e-03 |
| cubic_phase_omega_over_betatron_x | 3.85290114e-04 |
| cubic_phase_omega_over_betatron_y | 3.85609836e-04 |
| cubic_phase_Qs | 1.77896209e-04 |
| cubic_phase_Qs_smooth_phi_pi | 1.77896209e-04 |
| cubic_phase_voltage_MV | 1.88536813e-04 |

Generator `ξ` is `−c / (2π)` for the coefficient of `J δ`. `Λ = −2π (J_x ξ_x + J_y ξ_y) / C` uses that `ξ` and the placeholder actions of `02_rf_regimes.json`.

`|η1 δ_m / η0|`, `δ_m`, `|Λ / (η0 δ_m)|`, the linear `φ̂`, the cubic `φ̂²`, and `ω_s T_s` (also divided by the betatron phase) are much less than 1. The cubic `φ̂` is not required to be much less than 1.

## A linear J=0

Linear voltage, `J_x = J_y = 0`. The formula is chapter 1: `⟨δ²⟩ = δ_m²/2`, `Λ = 0`.

| target | value |
| --- | --- |
| delta_eq chapter 1-2 | 9.52953345e-09 |
| delta_eq chapter 3 | 9.52878896e-09 |
| delta_eq full pendulum | 9.52878888e-09 |
| second harmonic cos 2psi | -5.72370752e-09 |
| Omega/omega_s chapter 3 | 0.99984375 |
| Omega/omega_s pendulum | 0.99984375 |
| Lambda | 0.00000000e+00 |

NF `d · J` = 9.52953345e-09, `d_grid · J` = 9.52953345e-09, torus mean = 9.52953345e-09. `J_s` = 4.67070563e-05 from `δ_m² / (2 σ_δ²)` with `σ_δ` = 1.03465069e-02.

| quantity | formula | NF torus | degree-3 map | exact track |
| --- | --- | --- | --- | --- |
| delta_eq | 9.52953345e-09 | 9.52953345e-09 | 9.52935139e-09 | 9.52817969e-09 |
| delta_m half-range | 1.00000000e-04 | 9.99955318e-05 | 9.99965405e-05 | 9.99965410e-05 |
| second harmonic cos 2psi | -5.72370752e-09 | -5.72376469e-09 | -5.72338773e-09 | -5.72180968e-09 |
| second harmonic amplitude | 5.72370752e-09 | 5.72376469e-09 | 5.72338773e-09 | 5.72180968e-09 |
| Q / Q_note | 1.00000000 | 1.00000187 | 0.99984558 | 0.99984539 |
| mean square of delta | 5.00000000e-09 | 4.99921881e-09 | 4.99928487e-09 | 4.99928493e-09 |
| phi_hat | 5.00000000e-02 | nan | 5.00035181e-02 | 5.00035074e-02 |
| gamma_eff | 1.14395294e+00 | 1.14395294e+00 | 1.14395294e+00 | 1.14395294e+00 |
| Delta E / E | 2.69800573e-09 | 2.69793534e-09 | 2.69789835e-09 | 2.69762202e-09 |

Map residual rms 3.24059766e-13, exact residual rms 6.08084243e-13. Map `δ_eq` split 6.37737359e-17, exact split 6.82055006e-16. Turns: map 2811, exact 2811. Degree-3 image of the origin, delta component 0.00000000e+00.

## B linear finite J

Same linear voltage. `J_x` and `J_y` are the placeholder actions. `Λ` uses generator chromaticity.

| target | value |
| --- | --- |
| delta_eq chapter 1-2 | 1.90767267e-08 |
| delta_eq chapter 3 | 1.90759822e-08 |
| delta_eq full pendulum | 1.90759821e-08 |
| second harmonic cos 2psi | -5.72370752e-09 |
| Omega/omega_s chapter 3 | 0.99984375 |
| Omega/omega_s pendulum | 0.99984375 |
| Lambda | 5.09522845e-09 |

NF `d · J` = 1.90767267e-08, `d_grid · J` = 1.90767267e-08, torus mean = 1.90767267e-08. `J_s` = 4.67070563e-05 from `δ_m² / (2 σ_δ²)` with `σ_δ` = 1.03465069e-02.

| quantity | formula | NF torus | degree-3 map | exact track |
| --- | --- | --- | --- | --- |
| delta_eq | 1.90767267e-08 | 1.90767267e-08 | 1.90752571e-08 | 1.90740626e-08 |
| delta_m half-range | 1.00000000e-04 | 9.99955313e-05 | 9.99869934e-05 | 9.99869940e-05 |
| second harmonic cos 2psi | -5.72370752e-09 | -5.72376469e-09 | -5.72229984e-09 | -5.72072182e-09 |
| second harmonic amplitude | 5.72370752e-09 | 5.72376469e-09 | 5.72229984e-09 | 5.72072182e-09 |
| Q / Q_note | 1.00000000 | 1.00000187 | 0.99984561 | 0.99984542 |
| mean square of delta | 5.00000000e-09 | 4.99921876e-09 | 4.99833040e-09 | 4.99833046e-09 |
| phi_hat | 5.00000000e-02 | nan | 4.99987657e-02 | 4.99987548e-02 |
| gamma_eff | 1.14395295e+00 | 1.14395295e+00 | 1.14395295e+00 | 1.14395295e+00 |
| Delta E / E | 4.94962618e-09 | 4.94955579e-09 | 4.94912916e-09 | 4.94884745e-09 |

Map residual rms 3.83111088e-13, exact residual rms 6.40989748e-13. Map `δ_eq` split 1.39377800e-15, exact split 2.35009032e-15. Turns: map 2811, exact 2811. Degree-3 image of the origin, delta component 0.00000000e+00.

## C cubic J=0

Cubic-phase voltage, `J = 0`. The formula column uses `⟨δ²⟩ = (δ_m²/2) (1 − φ̂²/32)` and `Ω/ω_s = 1 − φ̂²/16`. The torus is still the chapter 1–2 value. The pendulum column is the full sine, kept as the size of the dropped `φ^5` term.

| target | value |
| --- | --- |
| delta_eq chapter 1-2 | 9.52953345e-09 |
| delta_eq chapter 3 | 9.50273164e-09 |
| delta_eq full pendulum | 9.50263071e-09 |
| second harmonic cos 2psi | -5.72370752e-09 |
| Omega/omega_s chapter 3 | 0.99437500 |
| Omega/omega_s pendulum | 0.99437761 |
| Lambda | 0.00000000e+00 |

NF `d · J` = 9.52953345e-09, `d_grid · J` = 9.52953345e-09, torus mean = 9.52953345e-09. `J_s` = 2.80243870e-04 from `δ_m² / (2 σ_δ²)` with `σ_δ` = 4.22393223e-03.

| quantity | formula | NF torus | degree-3 map | exact track |
| --- | --- | --- | --- | --- |
| delta_eq | 9.50273164e-09 | 9.52953345e-09 | 9.54498645e-09 | 9.50125284e-09 |
| delta_m half-range | 1.00000000e-04 | 9.98593597e-05 | 9.99962098e-05 | 9.99962094e-05 |
| second harmonic cos 2psi | -5.72370752e-09 | -5.72370911e-09 | -5.72461929e-09 | -5.66769227e-09 |
| second harmonic amplitude | 5.72370752e-09 | 5.72370911e-09 | 5.72461929e-09 | 5.66769228e-09 |
| Q / Q_note | 0.99437500 | 1.00000005 | 0.99431377 | 0.99433522 |
| mean square of delta | 4.98593750e-09 | 4.97192446e-09 | 4.98532739e-09 | 4.98539996e-09 |
| phi_hat | 3.00000000e-01 | nan | 3.01128529e-01 | 3.01124980e-01 |
| gamma_eff | 1.14395294e+00 | 1.14395294e+00 | 1.14395294e+00 | 1.14395294e+00 |
| Delta E / E | 2.69041759e-09 | 2.69547584e-09 | 2.70032803e-09 | 2.69002039e-09 |

Map residual rms 8.12597653e-11, exact residual rms 1.44502925e-10. Map `δ_eq` split 1.36776944e-12, exact split 2.42291458e-12. Turns: map 16864, exact 16864. Degree-3 image of the origin, delta component 0.00000000e+00.

The run stopped on a physical disagreement. Numbers were not adjusted.
