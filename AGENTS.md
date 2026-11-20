# AGENTS.md

Rules for sessions that change SCT_2. The human workflow is [`CONTRIBUTING.md`](CONTRIBUTING.md): one issue, one branch, one pull request.

## What SCT_2 is

SCT_2 studies storage-ring physics that depends on BELL as a library: equilibrium momentum, chromaticity, sextupole layouts, dynamic aperture, and later spin coherence. The physics of magnets and maps lives in BELL. SCT_2 builds studies, tables, and figures on top of the current BELL API.

## Hard rules

- Do not change the BELL repository: no code, tests, docs, or commits there. A BELL gap is a finding: write a short English issue draft with a minimal reproducer and a recommendation. Open that issue in BELL; do not patch around it silently inside SCT_2 and do not implement the fix in BELL from this repo.
- Physics readings for optics, chromaticity, slip, and detuning come from Lie factors (Hamiltonian → Dragt–Finn). Tracking and DA are cross-checks only.
- Results come from real runs only. Never fit or tune numbers to match expectations. A physical disagreement (formula vs tracking, two BELL readings disagreeing) is a stop: report it, do not pick a winner.
- Do not force-push, and do not rewrite the history of `main`.
- Do not touch another person's uncommitted files.
- The working tree often has many unrelated untracked files. Stage only your paths. `git add -A` is forbidden.
- Use Conventional Commits, as in `CONTRIBUTING.md`.
- Write the repository in English: source, documentation, issues, pull requests, and commits.
- Do not name other accelerator codes in new files, issues, pull requests, or commits. Leave old folders as they are.
- Long tracking runs only behind an explicit flag; print progress unbuffered (`flush`).
- Do not commit raw particle dumps. Commit short-run tables and figures.

## Magnetic ring study

The active study lives under `studies/magnetic_ring/`. It uses the `magnetic_2` lattice through BELL. Notes `notes/equilibrium-energy-v3` sections 1–3 are the physics reference. Sections 4–6 (electrostatics) are out of scope for that study. Spin is out of scope for now.

Do not import or patch `Julia/equilibrium_energy`. Reuse ideas on the current BELL API.

## Orchestration

One study step is one issue, one agent task, one PR. The next step starts only after the previous numbers table is accepted.
