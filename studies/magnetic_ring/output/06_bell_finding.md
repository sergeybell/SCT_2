# BELL finding: exact_track_element has no method for EdgeMap

Observed while tracking the closed magnetic ring with a dipole fringe model (`studies/magnetic_ring/lattices/magnetic_2_closed.bell`, every BM angle π/8, sextupoles at the step-5 own strengths).

## What happened

`fringe_sequence` inserts an `EdgeMap` before and after every dipole whose model is not `NoEdge`. `EdgeMap` is not an `AbstractElement`. `exact_track_element` is defined for `Marker`, `Drift`, `RFCavity`, and `AbstractElement`.

Calling it on the entrance edge of a BM with `HardEdge()` raises `MethodError`.

```
MethodError: no method matching exact_track_element(::Vector{Float64}, ::EdgeMap, ::RefParticle, ::Float64; ds::Float64)
The function `exact_track_element` exists, but no method is defined for this combination of argument types.

```

The edge item in that call is `BM_in`, parent `BM`, side `in`. The loaded BM card has L = 5.32941480e+00 m, e1 = 0.00000000e+00, e2 = 0.00000000e+00, gap = 0.00000000e+00 m, FINT = 0.00000000e+00.

The one-turn Lie factor already includes these edges: `ring_lie_factor(lattice, model)` builds a Lie factor of `edge_map` for each `EdgeMap`. There is no matching exact, untruncated tracker.

## Reproducer

```julia
using Bell
lat = load_structure("studies/magnetic_ring/lattices/magnetic_2_closed.bell")
model = FringeModel(dipole = HardEdge())
seq = fringe_sequence(lat, model)
edge = first(it for it in seq if it isa EdgeMap)
exact_track_element(zeros(6), edge, lat.ref, path_length(lat))
```

## What this study did

The study does not add a method. A turn of `fringe_sequence` applies `dragt_finn_factor(edge_map(edge, ref))` on each `EdgeMap` (`apply_factor`), `exact_track_element` on each magnet, and `rf_kick` on the thin cavity. The edge piece is the Lie factor at map order 3. The magnet body is the exact flow.

## Recommendation

Give `EdgeMap` an exact tracker that integrates the edge Hamiltonian without a polynomial cut, and accept that item in `exact_track_element` (and in `exact_track_turn`). The Lie factor of `edge_map` would stay the map used by `ring_lie_factor`. Until that method exists, a mixed turn has to branch on `EdgeMap` itself.
