# Linear optics of closed magnetic_2

Lattice file: `C:\GIT_REPOS\SCT_2\studies\magnetic_ring\lattices\magnetic_2_closed.bell`. Every sextupole has `k2 = 0`. Reference: deuteron, kinetic energy 270.0 MeV, charge number 1.0, harmonic 1.

gamma = 1.14395294, beta = 0.48563473, total energy = 2145.612943 MeV, brho = 3.475685 T m.

The one-turn matrix is `ring_matrix` (`linear_matrix` of each element, `R = exp(L S Hess H)`). The polynomial is `path_expansion` of `observables`, which builds `ring_lie_factor`. `periodic_optics` is taken on the matrix and on the Lie factor `R`. `optics_along` transports the matrix Twiss. Where two readings name the same quantity, both are printed and the run stops if they differ by more than 1.0e-8. Nothing was averaged and nothing was tuned.

Largest |matrix − Lie R| and the slip, compaction, and tune residuals are below 1.0e-08. Symplectic error of the matrix is 2.071e-14. `eta1` matches `eta1_formula` from the same expansion. `(1/C) ∫ D/ρ ds` differs from `alpha0` by 8.940e-08.

## One-turn numbers

`Qx` and `Qy` include the integer from `optics_along`. Fractional tunes are `mu/2π` from `periodic_optics` of the matrix. `alpha0`, `eta0`, `eta1`, and `gamma_t` are from `path_expansion`.

| quantity | value |
| --- | --- |
| C [m] | 1.17388296e+02 |
| Qx | 2.46172015e+00 |
| Qy | 2.46133732e+00 |
| Qx fractional | 4.61720150e-01 |
| Qy fractional | 4.61337323e-01 |
| alpha0 | 2.30470282e-01 |
| alpha1 | 3.39007701e-01 |
| eta0 | -5.33688628e-01 |
| eta1 | 1.01716073e+00 |
| gamma_t | 2.08301565e+00 |
| 1/gamma^2 | 7.64158910e-01 |
| (1/C) int D/rho ds | 2.30470372e-01 |
| beta_x(0) [m] | 1.93675224e+01 |
| alpha_x(0) | -2.51911777e-04 |
| beta_y(0) [m] | 4.05447047e+00 |
| alpha_y(0) | -1.88150489e-04 |
| D_x(0) [m] | 1.10941104e-02 |
| D_x'(0) | -5.06166051e-04 |
| D_y(0) [m] | 0.00000000e+00 |
| D_y'(0) | 0.00000000e+00 |
| max abs(D_y) [m] | 0.00000000e+00 |
| max abs(D_y') | 0.00000000e+00 |
| max beta_x [m] | 1.93677385e+01 |
| max beta_y [m] | 1.93646410e+01 |
| max D_x [m] | 1.17016103e+01 |
| min D_x [m] | 6.77386256e-03 |

Closure of `optics_along` (relative beta, absolute alpha, absolute dispersion, fractional tune):

| quantity | residual |
| --- | --- |
| betx | 4.58591657e-15 |
| alfx | 1.04063351e-14 |
| bety | 1.11721370e-14 |
| alfy | 2.07832341e-14 |
| D | 1.58884711e-12 |
| Dp | 4.81604773e-14 |
| Qx_frac | 7.82707232e-15 |
| Qy_frac | 1.12132525e-14 |

## Cell phase advance

An arc cell is the repeated block `QF BM QD (L = 0.12) QD BM QF` closed by one `L = 0.12` element, or by two `OS2` drifts where that element is split. Sextupoles are at `k2 = 0`. The in-ring phase is the advance of `optics_along` from the cell entrance to the cell exit. It is not set to 90 degrees. The periodic column is `mu` of `periodic_optics` on that cell's own matrix. Those two phases answer different questions: the periodic column assumes the cell is repeated by itself, and the in-ring column is the advance inside this ring, straights included in the periodic Twiss.

In-ring mean: mu_x = 1.57150893 rad (90.0408 deg), mu_y = 1.57122422 rad (90.0245 deg). Spread max−min: 1.347e-03 rad in x, 6.694e-03 rad in y.

| cell | s_start [m] | s_end [m] | mu_x [deg] | mu_y [deg] | periodic mu_x [deg] | periodic mu_y [deg] | elements |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.0600 | 11.7988 | 90.0650 | 90.2055 | 90.0443 | 90.0271 | QF BM QD OS1 QD BM QF SF2 |
| 2 | 11.7988 | 23.5377 | 90.0167 | 89.8435 | 90.0443 | 90.0271 | QF BM QD SD QD BM QF SF1 |
| 3 | 23.5377 | 35.2765 | 90.0794 | 90.2163 | 90.0443 | 90.0271 | QF BM QD SD QD BM QF SF2 |
| 4 | 35.2765 | 47.0153 | 90.0022 | 89.8328 | 90.0443 | 90.0271 | QF BM QD OS1 QD BM QF OS2 OS2 |
| 5 | 58.7541 | 70.4930 | 90.0650 | 90.2055 | 90.0443 | 90.0271 | QF BM QD OS1 QD BM QF SF2 |
| 6 | 70.4930 | 82.2318 | 90.0167 | 89.8435 | 90.0443 | 90.0271 | QF BM QD SD QD BM QF SF1 |
| 7 | 82.2318 | 93.9706 | 90.0794 | 90.2163 | 90.0443 | 90.0271 | QF BM QD SD QD BM QF SF2 |
| 8 | 93.9706 | 105.7095 | 90.0022 | 89.8328 | 90.0443 | 90.0271 | QF BM QD OS1 QD BM QF OS2 OS2 |

The element list is two copies of 42 elements. Phase over the first copy: mu_x = 7.73372194 rad (443.1096 deg), mu_y = 7.73251925 rad (443.0407 deg), equal to half the one-turn phase within 1e-8 rad.

Straights are the gaps longer than 1 m between arc cells. Their phase is included so that cells plus gaps sum to the one-turn phase.

| after cell | length [m] | mu_x [deg] | mu_y [deg] |
| --- | --- | --- | --- |
| 4 | 11.7388 | 82.9463 | 82.9427 |
| 8 | 11.7388 | 82.9463 | 82.9427 |

## Maxima

Local maxima of `D_x`, `beta_x`, and `beta_y` on the `optics_along` table, at least 0.8 of that curve's global maximum, and at least 0.5 m apart. A maximum on the periodic boundary `s = 0` is kept when it is above both neighbours. `dmu_to_next` is the phase to the next maximum of the same quantity, using `mux` for `D_x` and `beta_x` and `muy` for `beta_y`. The step that crosses `s = C` is included, so the advances of one quantity sum to the one-turn phase.

`find_optics_maxima` applies the same separation and fraction, and also requires `D_x > 0.5` m. The last column says whether that call returned the same peak. A step can exceed 360 degrees because the advance is betatron phase and the tune is above 2.

Curve maxima: 22. `find_optics_maxima` rows: 16. Curve maxima also returned by that call: 16.

| quantity | s [m] | value | beta_x [m] | beta_y [m] | D_x [m] | phase to next [deg] | find_optics_maxima |
| --- | --- | --- | --- | --- | --- | --- | --- |
| D_x | 23.4177 | 11.701610 | 19.3677 | 4.0554 | 11.7016 | 443.4646 | yes |
| D_x | 82.2318 | 11.701610 | 19.3677 | 4.0554 | 11.7016 | 442.7546 | yes |
| beta_x | 0.0600 | 19.367739 | 19.3677 | 4.0554 | 0.0111 | 89.6472 | no |
| beta_x | 11.6788 | 16.458656 | 16.4587 | 3.3248 | 5.8322 | 90.0794 | yes |
| beta_x | 23.4177 | 19.367715 | 19.3677 | 4.0554 | 11.7016 | 90.0167 | yes |
| beta_x | 35.1565 | 16.458626 | 16.4586 | 3.3248 | 5.8715 | 90.0650 | yes |
| beta_x | 46.8953 | 19.367739 | 19.3677 | 4.0554 | 0.0111 | 82.9463 | no |
| beta_x | 58.6341 | 19.367678 | 19.3677 | 4.0553 | 0.0111 | 90.0022 | no |
| beta_x | 70.3730 | 16.458656 | 16.4587 | 3.3248 | 5.8322 | 90.0794 | yes |
| beta_x | 82.1118 | 19.367715 | 19.3677 | 4.0554 | 11.7016 | 90.0167 | yes |
| beta_x | 93.8506 | 16.458626 | 16.4586 | 3.3248 | 5.8715 | 90.0650 | yes |
| beta_x | 105.5895 | 19.367739 | 19.3677 | 4.0554 | 0.0111 | 83.3013 | no |
| beta_y | 5.8094 | 18.783687 | 3.6607 | 18.7837 | 1.1278 | 96.0773 | yes |
| beta_y | 17.6682 | 18.780173 | 3.6597 | 18.7802 | 4.9218 | 83.9765 | yes |
| beta_y | 29.2871 | 18.780173 | 3.6597 | 18.7802 | 4.9218 | 96.0773 | yes |
| beta_y | 41.1459 | 18.783687 | 3.6607 | 18.7837 | 1.1278 | 83.2772 | yes |
| beta_y | 52.7647 | 19.364641 | 4.0544 | 19.3646 | 0.0068 | 83.6323 | no |
| beta_y | 64.5036 | 18.783687 | 3.6607 | 18.7837 | 1.1278 | 96.0773 | yes |
| beta_y | 76.3624 | 18.780173 | 3.6597 | 18.7802 | 4.9218 | 83.9765 | yes |
| beta_y | 87.9812 | 18.780173 | 3.6597 | 18.7802 | 4.9218 | 96.0773 | yes |
| beta_y | 99.8401 | 18.783687 | 3.6607 | 18.7837 | 1.1278 | 83.2772 | yes |
| beta_y | 111.4589 | 19.364641 | 4.0544 | 19.3646 | 0.0068 | 83.6323 | no |

Figure: `02_optics.svg`. Samples: `02_optics.csv` (ds = 0.02 m).

## RF regimes

Two stationary buckets share `h = 1` and `delta_m = 1.00000000e-04`. The linear voltage realises `phi_hat = 0.05`. The cubic-phase voltage realises `phi_hat = 0.30`. Full tables, the placeholder action, and the inequality values are in `02_rf_regimes.md` and `02_rf_regimes.json`.

| regime | h | V [MV] | delta_m | phi_hat | abs(eta1*dm/eta0) | abs(Lambda/(eta0*dm)) | omega_s T_s | Q_s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| linear | 1 | 6.78732526e-03 | 1.00000000e-04 | 5.00000000e-02 | 1.90590669e-04 | 9.54720568e-05 | 6.70652909e-03 | 1.06737726e-03 |
| cubic_phase | 1 | 1.88536813e-04 | 1.00000000e-04 | 3.00000000e-01 | 1.90590669e-04 | 9.54720568e-05 | 1.11775485e-03 | 1.77896209e-04 |

Natural chromaticity used only inside Lambda: xi_x = -1.461538, xi_y = -2.225488 (first derivative from `observables`, k2 = 0).
