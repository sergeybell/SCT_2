# Magnetic ring

Storage-ring study on the current BELL API. Step 1 writes the `magnetic_2` element table and closes the sector bends onto 2π.

BELL is a path dependency at `C:/GIT_REPOS/BELL`. This study reads `examples/lattices/magnetic_2.bell` and does not modify that tree.

## Layout

- `src/MagneticRing.jl` loads the example, sets every BM angle to π/8, and reloads the closed ring.
- `lattices/magnetic_2_closed.bell` is the closed lattice later steps load with `load_closed_lattice()`.
- `scripts/01_lattice.jl` writes the element table, the sequence, and the before/after readings.
- `scripts/02_linear_optics.jl` reads the closed lattice: tunes, slip factor, cell phase, optics along `s`, and two stationary RF regimes.
- `scripts/02b_generator_chromaticity.jl` reads first and second chromaticity from the RF-off Lie generators.
- `output/` holds the markdown reports, the csv tables, and the figures.

## Setup

From the repository root:

```bat
julia --project=studies/magnetic_ring -e "using Pkg; Pkg.develop(path=\"C:/GIT_REPOS/BELL\")"
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/01_lattice.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/02_linear_optics.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/02b_generator_chromaticity.jl
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/03_lie_invariants.jl
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
