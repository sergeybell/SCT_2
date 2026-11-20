# Magnetic ring

Storage-ring study on the current BELL API. Step 1 writes the `magnetic_2` element table and closes the sector bends onto 2π.

BELL is a path dependency at `C:/GIT_REPOS/BELL`. This study reads `examples/lattices/magnetic_2.bell` and does not modify that tree.

## Layout

- `src/MagneticRing.jl` loads the example, sets every BM angle to π/8, and reloads the closed ring.
- `lattices/magnetic_2_closed.bell` is the closed lattice later steps load with `load_closed_lattice()`.
- `scripts/01_lattice.jl` writes the element table, the sequence, and the before/after readings.
- `output/` holds the markdown report, the csv tables, and the bend-deficit figure.

## Setup

From the repository root:

```bat
julia --project=studies/magnetic_ring -e "using Pkg; Pkg.develop(path=\"C:/GIT_REPOS/BELL\")"
julia --project=studies/magnetic_ring studies/magnetic_ring/scripts/01_lattice.jl
```

## What the step changes

The example has 16 sector bends named BM and 4 drifts named OBM. BM and OBM share one length. The BM card angle is not π/8, so the sixteen bends sum to slightly less than 2π. The closed copy sets all 16 BM angles to π/8. OBM stays a drift. Element lengths are unchanged, so the path length is unchanged.

`Qx` and `Qy` in the output are the accumulated tunes from `optics_along` (the integer part is included). `Qx_frac` and `Qy_frac` are μ/2π from `periodic_optics` of the one-turn matrix.
