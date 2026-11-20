# Magnetic ring

Storage-ring study on the current BELL API. Step 1 writes the `magnetic_2` element table and closes the sector bends onto 2π.

BELL is a path dependency at `C:/GIT_REPOS/BELL`. This study reads `examples/lattices/magnetic_2.bell` and does not modify that tree.

## Layout

- `src/MagneticRing.jl` loads the example, sets every BM angle to π/8, and reloads the closed ring.
- `lattices/magnetic_2_closed.bell` is the closed lattice later steps load with `load_closed_lattice()`.
- `scripts/01_lattice.jl` writes the element table, the sequence, and the before/after readings.
- `scripts/02_linear_optics.jl` reads the closed lattice: tunes, slip factor, cell phase, optics along `s`, and two stationary RF regimes.
- `scripts/02b_generator_chromaticity.jl` reads first and second chromaticity from the RF-off Lie generators.
- `scripts/06_fringe.jl` reads three dipole fringe models on the own sextupole solution.
- `output/` holds the markdown reports, the csv tables, and the figures.

## Setup

From the repository root:

```bat
julia --project=studies/magnetic_ring -e "using Pkg; Pkg.develop(path=\"C:/GIT_REPOS/BELL\")"
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/01_lattice.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/02_linear_optics.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/02b_generator_chromaticity.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/03_lie_invariants.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/04_equilibrium.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/05_sextupoles.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/06_fringe.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/07_dynamic_aperture.jl
```

## Step 1

The example has 16 sector bends named BM and 4 drifts named OBM. BM and OBM share one length. The BM card angle is not π/8, so the sixteen bends sum to slightly less than 2π. The closed copy sets all 16 BM angles to π/8. OBM stays a drift. Element lengths are unchanged, so the path length is unchanged.

`Qx` and `Qy` in the output are the accumulated tunes from `optics_along` (the integer part is included). `Qx_frac` and `Qy_frac` are μ/2π from `periodic_optics` of the one-turn matrix.

## Step 2

`02_linear_optics.jl` loads `lattices/magnetic_2_closed.bell` with sextupoles at `k2 = 0`. Tunes and Twiss come from `ring_matrix`, `periodic_optics`, and `optics_along`. `alpha0`, `eta0`, and `eta1` come from `path_expansion` of the one-turn Lie factor (`observables`). The arc-cell phase is the advance measured along the ring, not an assumed right angle.

`output/02_rf_regimes.json` stores two stationary buckets for later steps. They share the lattice harmonic and one momentum half-range. The linear regime targets a phase half-range 0.05. The cubic-phase regime targets 0.30.

## Step 2b

`02b_generator_chromaticity.jl` builds the normalising matrix of the RF-off one-turn factor: a symplectic dispersion shear, including the path-length row, then the Courant–Snyder matrix of each plane at the ring start. First chromaticity is the coefficient of `J δ` in the pulled cubic generator. Second chromaticity is the coefficient of `J δ²` in `f4 − ½[f3, F3]`, after the non-resonant cubic is removed. Both are compared with `chromaticity_orders` and `chromaticity_series`, and the first chromaticity is compared with the feed-down of one sextupole. The overall factor in front of the generator coefficient is the member of a fixed list that matches `chromaticity_orders` on the bare ring.

## Step 3

`03_lie_invariants.jl` reads `α0`, `α1`, `η0`, and `η1` from `path_expansion` of the RF-off factor and compares `η1` with `α1 − α0/γ² + (3γ² − 1)/(2γ⁴)`. Generator chromaticity, with the factor from step 2b, is placed next to `chromaticity_orders`. A thin cavity is then appended at the linear-regime voltage and harmonic, with BELL synchronous phase `π`. The one-turn factor is `ring_lie_factor` of that lattice. The 6D normal form supplies `μ`, `Q_s = |μ_s|/2π`, and the detuning matrix of `h4`. A short track uses `exact_track_element` on the magnets and `rf_kick` on the cavity, and writes normal-form actions beside the Courant–Snyder actions at the ring start. The default track is two synchrotron periods and at most 2000 turns. `--turns N` is the explicit longer run.

## Step 4

`04_equilibrium.jl` reads `δ_eq` for the two regimes in `output/02_rf_regimes.json`. The cavity phase in BELL is `π` (early-positive `l`); the notes' stationary phase below transition is `0`, and the reported phase is the notes' sign. A short track at `π` stays in the bucket. A control at `0` leaves it. Each physics case then has four columns: the note formula, the normal-form torus average, iteration of the degree-3 polynomial of `ring_lie_factor`, and `exact_track_element` plus `rf_kick`. The normal form is compared with chapters 1–2. The cubic-phase case compares the map and the exact track with `1 − φ̂²/32` and `Ω/ω_s = 1 − φ̂²/16`. The default track is three synchrotron periods. `--turns N` replaces that length.

## Step 5

`05_sextupoles.jl` places SXD, SXB, and SYB with `place_sextupole_families` on `optics_along` of the closed ring. The in-study response is a central difference of generator chromaticity and of `η1` from `path_expansion`. Its solution sets `ξx = ξy = η1 = 0`, and prints the `α1` that identity requires. `solve_A` is the neighbour: `ring_metrics` and `response` set `ξx = ξy = α1 = 0` and leave `η1`. Both settings, and `k2 = 0`, are read in the linear regime of `output/02_rf_regimes.json` with the cavity phase `π`. A short exact track gives `δ_eq` and `γ_eff` beside the formula. Strengths for later steps are `output/05_sext_strengths.json`. Shared slots and the `D > 0.5` m mask are recorded, not patched. `--turns N` replaces the track length.

## Step 6

`06_fringe.jl` rebuilds that placement and sets `k2_per_m3.own`. It does not solve the strengths again. Three dipole fringe models are read with `ring_lie_factor`: `FringeModel()` (NoEdge), `HardEdge()`, and `SoftEdge(profile = LinearRamp(0.05))`. Quadrupoles and sextupoles stay NoEdge. The ramp width 0.05 m is an assumption.

`Qx` and `Qy` are the body tunes from `optics_along` plus the fringe shift of μ/2π. Chromaticity is `generator_chromaticity` of that factor; the normalising matrix is built from the factor's own `periodic_optics`. `η1` is `path_expansion`. The linear-regime `δ_eq` uses that η1 and ξ with the step-4/5 values of δ_m, Jx, and Jy.

`exact_track_element` has no method for `EdgeMap`. The short track applies the Lie factor of `edge_map` on each edge, `exact_track_element` on each magnet, and `rf_kick` on the thin cavity. `--turns N` replaces the default length of one synchrotron period. The BELL finding is `output/06_bell_finding.md`.

## Step 7

`07_dynamic_aperture.jl` rebuilds the step-5 placement and sets `k2_per_m3.own`. The fringe model is `FringeModel()` (NoEdge). Tunes, resonance distances, and max `|f3|` are read from that magnetic factor. `dragt_finn_nf6d` on the magnetic factor is rejected because the longitudinal block is a shear. The optional `h4` matrix is `dragt_finn_nf6d` of `one_turn_factor_rf` at the linear-regime voltage, with BELL phase `π`. It is not an input of the aperture scan.

The Guignard boundary approximates `{Q_x}` by `p/q` with `q ≤ 100`, walks the linear orbit of `F.R` from `(x0, 0, 0, 0, 0, δ)`, and multiplies the horizontal 2×2 blocks of `jacobian_of_map` of the degree-2 map and of the degree-3 map. The boundary is the first `|x0|` at which `|Tr|` leaves the open interval `(−2, 2)`. The scan is repeated at `δ = 0` and at `±δ_m`.

Survival is `exact_track_turn`. A particle survives when `x` and `y` at the ring start stay finite and strictly inside 5 cm. The default is 300 turns. `LONG_DA=1` or `--long` uses 2000 turns and prints progress with flush. `07_da.svg` overlays the raw `|Tr| = 2` lines and the survival samples. A disagreement is left as it is.
