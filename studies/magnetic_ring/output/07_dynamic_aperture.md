# Dynamic aperture

NoEdge (`FringeModel()`) on the step-5 own sextupoles. The one-turn object is `ring_lie_factor`. Guignard uses `factor_to_map` truncated at degree 2 and degree 3, then `jacobian_of_map` on the linear orbit of `F.R`. Survival is `exact_track_turn` at the ring start. The two are overlaid and not averaged.

Search window for the trace: 0 to 2.0 m. Aperture cut: 0.05 m. Turns: 300. Step: 0.1 m. q ≤ 100. Runtime 763.9 s.

## Lie readings

Tunes are `optics_along` (integer part included) and μ/2π from `periodic_optics` of the factor. `{Q_x}` is approximated by 6/13, the fraction with q ≤ 100 that minimises |q {Q_x} − p|.

| quantity | value |
| --- | --- |
| Qx | 2.461720e+00 |
| Qy | 2.461337e+00 |
| Qx fractional | 4.617201e-01 |
| Qy fractional | 4.613373e-01 |
| p/q | 6/13 |
| |q Qx_frac − p| | 2.361948e-03 |
| |q Qy_frac − nearest integer| | 2.614800e-03 |
| max |f3| | 1.163649e+02 |
| max |f4| | 7.160297e+02 |
| 2 cos(q μx) | 1.999780e+00 |
| Tr at x0 = 0, degree 3 | 1.999780e+00 |
| αx | -2.519118e-04 |
| βx [m] | 1.936752e+01 |
| coupling | 0.000000e+00 |
| δ_m | 1.000000e-04 |

Distance of the fractional tune to 1/3 and to 1/4, and to the nearest multiple of those families:

| plane | to 1/3 | to 1/4 | nearest third | nearest quarter |
| --- | --- | --- | --- | --- |
| x | 1.283868e-01 | 2.117201e-01 | 1.283868e-01 | 3.827985e-02 |
| y | 1.280040e-01 | 2.113373e-01 | 1.280040e-01 | 3.866268e-02 |

Distance of n_x Q_x + n_y Q_y to the nearest integer. The tunes are the full `optics_along` values.

| resonance | |n·Q − integer| |
| --- | --- |
| 3Qx | 3.851604e-01 |
| 4Qx | 1.531194e-01 |
| 3Qy | 3.840120e-01 |
| 4Qy | 1.546507e-01 |
| 2Qx+Qy | 3.847776e-01 |
| 2Qx-Qy | 4.621030e-01 |
| Qx+2Qy | 3.843948e-01 |
| Qx-2Qy | 4.609545e-01 |
| Qx+Qy | 7.694253e-02 |
| Qx-Qy | 3.828268e-04 |
| 2Qx+2Qy | 1.538851e-01 |
| 2Qx | 7.655970e-02 |
| 2Qy | 7.732535e-02 |
| 3Qx+Qy | 1.535022e-01 |
| 3Qx-Qy | 7.617687e-02 |
| Qx+3Qy | 1.542679e-01 |
| Qx-3Qy | 7.770818e-02 |

Largest |coefficients| of f3:

| monomial | coefficient |
| --- | --- |
| px^3 | 1.163649e+02 |
| px delta^2 | 1.015677e+02 |
| delta^3 | -1.597187e+01 |
| px py^2 | -1.536451e+01 |
| x delta^2 | 5.892759e+00 |
| x px delta | -3.798027e+00 |
| x px^2 | 2.150015e+00 |
| px^2 delta | 1.809104e+00 |

The magnetic factor was not accepted by `dragt_finn_nf6d`: one-turn matrix not stable. The longitudinal block is a shear, so this is the documented rejection (fewer than three elliptic pairs), not a second tune reading.

Optional normal form of `one_turn_factor_rf` at the linear-regime voltage 0.00678732526132278 MV, harmonic 1, BELL phase π. This h4 is not an input of the Guignard scan. `Q_s = |μ_s| / 2π = 1.067379e-03`. max |f3| before the homological step is 1.396936e+00; the cubic residual is 2.220446e-16; |f4 − h4| is 5.868084e-13.

| c_i1 | c_i2 | c_i3 |
| --- | --- | --- |
| -8.587587e+00 | 6.895292e+00 | -8.784327e-04 |
| 0.000000e+00 | -4.053072e-01 | 7.811269e-03 |
| 0.000000e+00 | 0.000000e+00 | 1.121570e-02 |

## Guignard boundary

The linear orbit is R applied to (x0, 0, 0, 0, 0, δ). Rprod is the product of the horizontal 2×2 blocks of the one-turn Jacobian. The boundary below is the stable side of the first crossing of |Tr| = 2. A sample that is stable again after that crossing is listed as re-entrant and does not move the boundary.

| degree | δ | sign | status | first x0 [m] | Tr stable side | largest stable sample [m] |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 0.000000e+00 | + | first crossing | 8.777416e-03 | 1.9999999999999998e+00 | 1.000000e+00 |
| 2 | 0.000000e+00 | − | first crossing | -8.860764e-03 | 1.9999999999999996e+00 | -9.500000e-01 |
| 2 | 1.000000e-04 | + | first crossing | 8.818307e-03 | 1.9999999999999998e+00 | 1.000000e+00 |
| 2 | 1.000000e-04 | − | first crossing | -8.899778e-03 | 1.9999999999999998e+00 | -9.500000e-01 |
| 2 | -1.000000e-04 | + | first crossing | 8.816194e-03 | 1.9999999999999996e+00 | 1.000000e+00 |
| 2 | -1.000000e-04 | − | first crossing | -8.901448e-03 | 1.9999999999999998e+00 | -9.500000e-01 |
| 3 | 0.000000e+00 | + | first crossing | 4.624863e-01 | -1.9999999999999925e+00 | same |
| 3 | 0.000000e+00 | − | first crossing | -4.621819e-01 | -1.9999999999999976e+00 | same |
| 3 | 1.000000e-04 | + | first crossing | 4.625543e-01 | -1.9999999999999987e+00 | same |
| 3 | 1.000000e-04 | − | first crossing | -4.622474e-01 | -1.9999999999999929e+00 | same |
| 3 | -1.000000e-04 | + | first crossing | 4.624149e-01 | -1.9999999999999918e+00 | same |
| 3 | -1.000000e-04 | − | first crossing | -4.621130e-01 | -1.9999999999999938e+00 | same |

The reported boundary is that first crossing, the edge of the stable interval that contains x0 = 0. Where the largest stable grid sample differs, the trace has become unstable and then stable again. That later sample is not used as the boundary.

## Tracking survival

Launches: 83. Survived all 300 turns: 83. Fraction: 1.000000e+00. The cut is |x| < 0.05 m and |y| < 0.05 m at the ring start, and a finite state vector. The initial condition is the Guignard launch (x0, 0, y0, 0, 0, δ), not a matched Courant–Snyder momentum.

Every launch survived. The outermost survivor is the edge of this grid. No dynamic-aperture loss was observed inside the aperture cut.

| y [m] | δ | sign | last survivor connected to 0 [m] | first loss [m] | re-entrant survivors |
| --- | --- | --- | --- | --- | --- |
| -2.000000e-02 | 0.000000e+00 | − | -2.000000e-02 | none | no |
| -2.000000e-02 | 0.000000e+00 | + | 2.000000e-02 | none | no |
| -1.000000e-02 | 0.000000e+00 | − | -2.000000e-02 | none | no |
| -1.000000e-02 | 0.000000e+00 | + | 2.000000e-02 | none | no |
| 0.000000e+00 | -1.000000e-04 | − | -4.500000e-02 | none | no |
| 0.000000e+00 | -1.000000e-04 | + | 4.500000e-02 | none | no |
| 0.000000e+00 | 0.000000e+00 | − | -4.500000e-02 | none | no |
| 0.000000e+00 | 0.000000e+00 | + | 4.500000e-02 | none | no |
| 0.000000e+00 | 1.000000e-04 | − | -4.500000e-02 | none | no |
| 0.000000e+00 | 1.000000e-04 | + | 4.500000e-02 | none | no |
| 1.000000e-02 | 0.000000e+00 | − | -2.000000e-02 | none | no |
| 1.000000e-02 | 0.000000e+00 | + | 2.000000e-02 | none | no |
| 2.000000e-02 | 0.000000e+00 | − | -2.000000e-02 | none | no |
| 2.000000e-02 | 0.000000e+00 | + | 2.000000e-02 | none | no |

Off-axis launches (y ≠ 0, δ = 0): 20, survived 20.

Phase portraits are the raw turn samples of the exact track at x0 = 5, 10, 20 and 40 mm, y = 0, δ = 0.

## Figures

- `07_da.svg` — degree-2 and degree-3 traces at δ = 0, and the exact survival samples at y = 0. Vertical lines are the |Tr| = 2 boundaries that fall inside the aperture window. Nothing is interpolated between a trace and a survival point.
- `07_phase.svg` — phase portraits.

Tables: `07_numbers.csv`, `07_guignard.csv`, `07_survival.csv`, `07_phase.csv`.
