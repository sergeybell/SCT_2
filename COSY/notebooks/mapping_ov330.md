# Chromaticity from `OV 3 3 0` (magnetic_2)

Test mapping in the same layout as the SCT article scan
(`SCT/COSY/jobs/mapping.fox`, notebooks `SCT_mapping` / `Mapping_criterion`),
with one change: **do not mix two DA settings**.

## Two quantities that were both called chromaticity

| | \(\xi_B\) Berz / COSY `TP` | \(\xi_J\) particle / `MAP(5)` |
|---|---|---|
| COSY | `OV 3 2 1`, `SET_FOR_DEUTERONS_CHROM`, `TP` | `OV 3 3 0`, `SET_FOR_DEUTERONS`, ellipse of `MAP(5)` |
| What varies | energy of the **closed orbit** \(z(\delta)\), equivalently \(B\rho\) | a particle with **betatron amplitude** on a **fixed** lattice |
| Definition | \(\partial\nu/\partial\delta\) around \(z(\delta)\) | \(\Lambda/J = -2\pi\xi_J/C\) from time-of-flight |
| SCT article | this sat on the \(x\)-axis of \(\Delta\nu_s(\xi)\) | not used |
| This test | **not computed** | **this is `chrom_x`, `chrom_y`** |

They are not two approximations of one number. Berz asks how the working point
of the beam moves when the whole beam energy (or the magnet fields) changes.
The particle quantity asks how much extra path a transversely offset particle
accumulates, which is what the cavity sees and what shifts \(\delta_{\mathrm{eq}}\).

Amplitude-dependent orbital tunes (Berz §7.3, \(\partial\nu_x/\partial J\)) are
not used. `quad_Kx` in the file is still the **spin** DA coefficient
\(\partial^2\nu_s/\partial x^2\) from `TSS`.

## RF: not the lattice `RFFLAG`

`magnetic_2.fox` with `RFFLAG=1` is not this ring's cavity:

- `FREQ = h * REVFREQ(fACCLEN)` and `fACCLEN=141` in `header.fox`, while
  `SPOS = 117.39` m. Phase slip per turn \(\sim 60^\circ\). The expansion
  point is not a 6D closed orbit. The \(\delta\)-probe then gave a fake
  \(\Delta\nu_{s,z}\sim 3.8\times 10^{-6}\) (almost independent of sextupoles),
  and a linear fit of \(\Delta\nu_s(\xi_J,\eta_1)\) looked broken.
- Voltage \(5\,\mathrm{kV}\), \(h=1\) makes \(Q_s\) so small that COSY `GT`
  sees \(\cos\mu\approx 1\) with nonzero off-diagonals and quits
  (`unrepresentable map`). That is why `RFFLAG=1` cannot sit on the same
  map as `PATHXI`/`TSS`.

The job therefore follows `delta_eq.fox`: lattice with `RFFLAG=0`, then a
thin cavity \(h=5\), \(200\,\mathrm{kV}\), `FREQ = h*REVFREQ(SPOS)`.

`PATHXI` (ellipse of `MAP(5)`, \(\eta\), \(\alpha\)) is read **before** the
cavity: `ECM`/`MCM` on a 6D RF map also trip `GT`. Spin (`TSS`, P0/P1 probes)
is read **after** the cavity. Still one DA setting, `OV 3 3 0` only.

On four strong-sextupole nodes `TSS` hits a spin-orbit resonance (`NORMSU2`,
`MU=0`). Those points are dropped in the plots; they are not physics of
\(\xi_J\).

## Output

Under `COSY/dat/magnetic_2/`:

- `integrals_ov330.dat` — \(\xi_J\), \(\eta\), \(\alpha\), spin-DA `quad_K*`
- `particle_spin_tune_ov330.dat` — P0 \(\Delta\nu_s\)
- `particle_spin_tune_1_ov330.dat` — P1 \(\Delta\nu_s\)

Plots: `COSY/dat/magnetic_2/fig_ov330_*.png`. Notebook: `COSY/notebooks/mapping_ov330.ipynb`.

At zero sextupoles \(\xi_J=(-1.473,-2.230)\), same as Julia 6D Lie and
`delta_eq.fox`. Berz/passport \((-1.84,-1.85)\) is \(\xi_B\) and is not on
these axes.

On the test grid the ellipse and the centred values stay together.
\(\xi_J\) moves by tens of units at the SCT steps `SGx*=0.02`, while
\(\xi_B\) in the article stayed \(O(1)\). That is the two derivatives
reacting differently to the same sextupole, not a bug in `PATHXI`.
After dropping the four `TSS` failures, \(\Delta\nu_s\) of the probes is
again linear in \((\xi_J,\eta_1)\) on this \(3^3\) scan, because both
sides are linear in the three knobs. The previous \(R^2\sim 0.2\)–\(0.4\)
was the mistuned lattice RF plus the `MU=0` points left in the fit.
