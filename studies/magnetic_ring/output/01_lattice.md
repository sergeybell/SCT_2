# magnetic_2 lattice closed to 2π

Source file (read-only): `C:/GIT_REPOS/BELL/examples/lattices/magnetic_2.bell`.

The file token on the BM card is `0.3926990753258494`. Loaded as a Float64 it is 0.3926990753258494. Sixteen equal sector bends therefore sum to 6.283185205213589, which is -1.0196599742329227e-7 relative to 2π. `OBM` is a drift of the same length and is left unchanged. Every BM angle in the closed ring is the Float64 value of π/8, 0.39269908169872414. `total_bend_angle` of that ring is 6.283185307179588, which differs from 2π by 1.7763568394002505e-15. The product 16 * (π/8) differs from 2π by 0.0.

Closed lattice: `C:\GIT_REPOS\SCT_2\studies\magnetic_ring\lattices\magnetic_2_closed.bell`. Reloading that file reproduces the in-memory closed ring, including the reference particle. A separate `write_structure` round trip does not: element-wise match false, BM angle change -5.551115123125783e-17 rad on each of 16 bends, total_bend_angle change -4.440892098500626e-15 rad. That bend change is not 16 times the angle change: total_bend_angle reduces over the whole element list. Path-length change 0.0 m, other element mismatches 0, reference match true. The saved lattice is the edited copy, not the `write_structure` text.

Counts: 16 BM, 4 OBM, shared length 5.329414799999999 m, 84 elements in sequence. Reference: deuteron, kinetic energy 270.0 MeV, harmonic 1.

`Qx` and `Qy` are accumulated tunes from `optics_along` (integer part included). `Qx_frac` and `Qy_frac` are μ/2π from `periodic_optics` of `ring_matrix`. Open-ring fractional closure residuals: Qx 5.717648576819556e-15, Qy 1.1213252548714081e-14. Closed-ring residuals: Qx 7.827072323607354e-15, Qy 1.1213252548714081e-14. Symplectic error of the one-turn matrix: open 3.201432791984271e-14, closed 2.0714340087555213e-14.

Path length is the sum of element lengths. The arc length of BM is not changed, so the circumference does not change.

## Closure

| quantity | before | after | delta |
| --- | --- | --- | --- |
| bm_angle_rad | 0.3926990753258494 | 0.39269908169872414 | 6.372874727933464e-9 |
| total_bend_angle_rad | 6.283185205213589 | 6.283185307179588 | 1.019659991996491e-7 |
| total_bend_minus_2pi_rad | -1.0196599742329227e-7 | 1.7763568394002505e-15 | 1.019659991996491e-7 |
| path_length_m | 117.38829599999988 | 117.38829599999988 | 0.0 |
| Qx | 2.4617201386366054 | 2.4617201498667933 | 1.1230187890731713e-8 |
| Qy | 2.4613373230490914 | 2.4613373230490914 | 0.0 |
| Qx_frac | 0.4617201386365997 | 0.46172014986680116 | 1.1230201435452614e-8 |
| Qy_frac | 0.46133732304908015 | 0.46133732304908015 | 0.0 |

## Unique elements

| name | count | type | length_m | angle_open_rad | angle_closed_rad | k1 | k2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| OS2 | 8 | Drift | 0.06 |  |  |  |  |
| QF | 16 | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| BM | 16 | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| QD | 16 | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| OS1 | 6 | Drift | 0.12 |  |  |  |  |
| SF2 | 4 | Sextupole | 0.12 |  |  |  | 0.0 |
| SD | 4 | Sextupole | 0.12 |  |  |  | 0.0 |
| SF1 | 2 | Sextupole | 0.12 |  |  |  | 0.0 |
| QF1 | 4 | Quadrupole | 0.21 |  |  | 0.5565959259588208 |  |
| OBM | 4 | Drift | 5.329414799999999 |  |  |  |  |
| QD1 | 4 | Quadrupole | 0.21 |  |  | -0.5565959259588208 |  |

## Sequence

| index | name | type | length_m | angle_open_rad | angle_closed_rad | k1 | k2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | OS2 | Drift | 0.06 |  |  |  |  |
| 2 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 3 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 4 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 5 | OS1 | Drift | 0.12 |  |  |  |  |
| 6 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 7 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 8 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 9 | SF2 | Sextupole | 0.12 |  |  |  | 0.0 |
| 10 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 11 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 12 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 13 | SD | Sextupole | 0.12 |  |  |  | 0.0 |
| 14 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 15 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 16 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 17 | SF1 | Sextupole | 0.12 |  |  |  | 0.0 |
| 18 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 19 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 20 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 21 | SD | Sextupole | 0.12 |  |  |  | 0.0 |
| 22 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 23 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 24 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 25 | SF2 | Sextupole | 0.12 |  |  |  | 0.0 |
| 26 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 27 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 28 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 29 | OS1 | Drift | 0.12 |  |  |  |  |
| 30 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 31 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 32 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 33 | OS2 | Drift | 0.06 |  |  |  |  |
| 34 | OS2 | Drift | 0.06 |  |  |  |  |
| 35 | QF1 | Quadrupole | 0.21 |  |  | 0.5565959259588208 |  |
| 36 | OBM | Drift | 5.329414799999999 |  |  |  |  |
| 37 | QD1 | Quadrupole | 0.21 |  |  | -0.5565959259588208 |  |
| 38 | OS1 | Drift | 0.12 |  |  |  |  |
| 39 | QD1 | Quadrupole | 0.21 |  |  | -0.5565959259588208 |  |
| 40 | OBM | Drift | 5.329414799999999 |  |  |  |  |
| 41 | QF1 | Quadrupole | 0.21 |  |  | 0.5565959259588208 |  |
| 42 | OS2 | Drift | 0.06 |  |  |  |  |
| 43 | OS2 | Drift | 0.06 |  |  |  |  |
| 44 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 45 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 46 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 47 | OS1 | Drift | 0.12 |  |  |  |  |
| 48 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 49 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 50 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 51 | SF2 | Sextupole | 0.12 |  |  |  | 0.0 |
| 52 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 53 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 54 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 55 | SD | Sextupole | 0.12 |  |  |  | 0.0 |
| 56 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 57 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 58 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 59 | SF1 | Sextupole | 0.12 |  |  |  | 0.0 |
| 60 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 61 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 62 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 63 | SD | Sextupole | 0.12 |  |  |  | 0.0 |
| 64 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 65 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 66 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 67 | SF2 | Sextupole | 0.12 |  |  |  | 0.0 |
| 68 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 69 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 70 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 71 | OS1 | Drift | 0.12 |  |  |  |  |
| 72 | QD | Quadrupole | 0.21 |  |  | -0.5805856428570646 |  |
| 73 | BM | SBend | 5.329414799999999 | 0.3926990753258494 | 0.39269908169872414 | 0.0 |  |
| 74 | QF | Quadrupole | 0.21 |  |  | 0.5212441755015996 |  |
| 75 | OS2 | Drift | 0.06 |  |  |  |  |
| 76 | OS2 | Drift | 0.06 |  |  |  |  |
| 77 | QF1 | Quadrupole | 0.21 |  |  | 0.5565959259588208 |  |
| 78 | OBM | Drift | 5.329414799999999 |  |  |  |  |
| 79 | QD1 | Quadrupole | 0.21 |  |  | -0.5565959259588208 |  |
| 80 | OS1 | Drift | 0.12 |  |  |  |  |
| 81 | QD1 | Quadrupole | 0.21 |  |  | -0.5565959259588208 |  |
| 82 | OBM | Drift | 5.329414799999999 |  |  |  |  |
| 83 | QF1 | Quadrupole | 0.21 |  |  | 0.5565959259588208 |  |
| 84 | OS2 | Drift | 0.06 |  |  |  |  |

