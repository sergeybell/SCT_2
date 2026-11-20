# BELL finding: dispersion mask and shared sextupole slots

Observed while placing SXD, SXB, and SYB on the closed magnetic ring (`studies/magnetic_ring/lattices/magnetic_2_closed.bell`, every BM angle π/8, sextupoles at k2 = 0).

## What happened

`find_optics_maxima` keeps a local maximum of D_x, beta_x, or beta_y only on samples with D > 0.5 m. The largest beta_x sample on this ring is at s = 5.87541480e+01 m with D = 1.10637404e-02 m, and it is not in that list (`kept = false`).

`place_sextupole_families` assigns each kept maximum to the nearest sextupole or short drift. The first family to claim a slot keeps it. Later maxima that fall in that slot are skipped.

Notes returned on this lattice:

- beta_x max at s=23.418 m shares slot SF1 with family SXD
- beta_x max at s=82.112 m shares slot SF1 with family SXD

The two horizontal peaks that sit on the dispersion maxima therefore do not become SXB magnets. SXB is placed on the remaining horizontal peaks. SXD is at the dispersion maxima. SYB is at the vertical maxima.

## Reproducer

```julia
lat = load_structure("studies/magnetic_ring/lattices/magnetic_2_closed.bell")
C = path_length(lat)
po = periodic_optics(ring_matrix(lat), lat.ref, C)
opt = optics_along(lat, po; ds = 0.02)
maxima = find_optics_maxima(opt)
placed, fam, rows, notes = place_sextupole_families(lat, opt)
```

## Recommendation

Document that the D > dmin mask applies to the beta curves as well as to D_x, and that a maximum whose nearest slot is already taken is omitted rather than moved to the next free slot. The returned note is the right report for a shared slot. No change to the response or the solve is requested from this run: the study uses the families and the notes the functions returned.
