# EquilibriumEnergy — study of δ_eq

Depends on Bell in `C:\GIT_REPOS\BELL\Julia`. Physics is in
`note/main.tex` (thin RF, Lie product of the ring, torus average of δ).
The OptiM file `C:\GIT_REPOS\SCT_2\OptiM\magnetic\magnetic_2.opt` is
converted on each run to `lattices/magnetic_2.bell` and loaded with
`load_structure`.

In this app the RF phase is reported with the opposite sign to Bell:
`φ_note > 0` means the particle is late.

```
julia --project=. setup.jl
julia --project=. scripts/run_v08.jl     # numbers and figures for the note
cd note; latexmk -xelatex main.tex
```

Notebook `notebooks/delta_eq.ipynb`: in Cursor pick the kernel
**Julia EquilibriumEnergy 1.13** (not Python). From a Julia REPL:

```julia
using IJulia
notebook(dir = "notebooks")
```

The notebook `notebooks/delta_eq.ipynb` walks the same calculation from one particle
to the four models. Markdown cells explain what the next code cell is about to show.
Short live tracking is in the notebook; the full $N=12000$ figures stay in `note/figs`.

Four models, all with a thin energy kick at the end of the turn:

- M1 analytic formula from time of flight
- M2 torus average of the Lie normal form (`one_turn_factor_rf`, `dragt_finn_nf6d`)
- M3 symplectic tracking of the unexpanded Hamiltonian
- M4 iteration of the degree-3 polynomial expansion of the same Lie product

Numbers in the note come from `note/gen/*.tex`, written by `run_v08.jl`.
