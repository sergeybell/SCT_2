# one_turn_factor_rf appends a kick the lattice may already contain

## Summary

`one_turn_factor_rf(lattice, cav, ref, path_length)` calls `ring_lie_factor(lattice)` and then concatenates `dragt_finn_factor(rf_kick_taylor(cav, ...))`. `ring_lie_factor` builds a factor for every element, and a thin `RFCavity` is the Lie map of `rf_potential`. On a lattice that already contains that cavity, the function applies the voltage twice.

## What this ring shows

Closed `magnetic_2`, sextupoles at `k2 = 0`, one thin cavity appended at the end. Voltage 6.78732526e-03 MV, harmonic 1, `phi_s = π`. The study factor is `ring_lie_factor` of that lattice. `one_turn_factor_rf` was not called.

Replacing the in-lattice cavity factor by one `rf_kick_taylor` factor changes the one-turn matrix by 2.842171e-14, f3 by 3.907985e-14, and f4 by 6.661338e-14. Concatenating that kick once more, which is the body of `one_turn_factor_rf` on this lattice, moves `|μ_s| / 2π` from 1.06737926e-03 to 1.50950505e-03. The ratio is 1.41421621e+00.

The normal form of the single in-lattice cavity returned `Q_s = 1.06737926e-03`.

## Recommendation

`one_turn_factor_rf` should take the magnetic lattice, via `magnetic_lattice`, before `ring_lie_factor`, or it should refuse a lattice that already contains an `RFCavity`. The docstring should state that a thin cavity already in the element list is a kick, so a second `rf_kick_taylor` doubles the voltage. Callers that have already appended the cavity should use `ring_lie_factor` and should not call `one_turn_factor_rf`.

## Reproducer

```julia
using Bell
using Bell.NormalForm: one_turn_factor_rf, symplectic_eigen_A
lat = load_structure(raw"C:/GIT_REPOS/SCT_2/studies/magnetic_ring/lattices/magnetic_2_closed.bell")
ref = lat.ref
C = path_length(lat)
cav = RFCavity("RF", 0.00678732526132278, 1, π)
elems = AbstractElement[e for e in lat.elements]
push!(elems, cav)
lat_rf = Lattice(lat.name * "+RF", elems, ref, lat.harmonic)
F = ring_lie_factor(lat_rf)                          # one kick
F2 = one_turn_factor_rf(lat_rf, cav, ref, C)         # two kicks
_, mu = symplectic_eigen_A(F.R)
_, mu2 = symplectic_eigen_A(F2.R)
abs(mu[3]) / 2π, abs(mu2[3]) / 2π
```

This note records the double count. It does not change BELL.
