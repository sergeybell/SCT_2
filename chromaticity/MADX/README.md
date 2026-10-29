# MAD-X 5.09.01 — хроматичность magnetic_2_md

Запуск: `run_madx.bat` (madx.exe из `C:\GIT_REPOS\NICA\MADX\2.Programms\Twiss`).
`magnetic_2_md.seq`: SBEND L=5.3294148 м, ANGLE=22.5°, E1=E2=0, FINT=0 (секторные, без краевой фокусировки),
квадруполи K1 = G[T/m]/Bρ (Bρ из BEAM), секступоли K2=0. BEAM: PARTICLE=DEUTERON, MASS=1.875612945 ГэВ,
ENERGY = MASS+0.270 ГэВ.

## (a) `chrom_plain.madx` → `output/madx_plain_results.txt`
TWISS (Q1, Q2, ALFA, DQ1, DQ2), TWISS с CHROM, и конечная разность TWISS при DELTAP=±1e-4.
В MAD-X DQ = dQ/dPT, PT = ΔE/(p0c). Так как dE = βc·dp, PT = β0·δ (1-й порядок), значит
**dQ/dδ = β0·dQ/dPT** (β0 = 0.485635). В файле: `DQ*_pt` (сырые) и `DQ*_dp` (по dp/p).
DELTAP в MAD-X — это dp/p, конечная разность даёт `DQ*_fd_dp` прямо по dp/p. eta = ALFA − 1/γ².

## (b) `chrom_ptc.madx` (по методу `Compaction.madx`)
PTC: `ptc_create_layout model=2 method=6 nst=20 exact time=false` (time=false ⇒ δ = dp/p),
`ptc_twiss icase=56 no=5` → `output/ptc_twiss_summary.tfs` (ALPHA_C, ALPHA_C_P…, ETA_C, Q1, Q2, DQ1, DQ2),
`ptc_normal icase=5 no=5` → `output/ptc_normal.tfs` (DQ порядков 1–3: производные dⁿQ/dδⁿ, ANHX/ANHY).
