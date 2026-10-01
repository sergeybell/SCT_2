# SCT_2
The goal of this project is to understand the common features of storage ring classes for spin dynamics.

Continuation of the SCT study: the mapping method and spin coherence time (SCT) estimation.
Current scope: magnetic lattices.

## Layout

```text
OptiM/magnetic/          OptiM lattice projects (*.opt)
COSY/structures/<stem>/  COSY lattices: <stem>.fox (LATTICE) and <stem>_maps.fox (LATTICE + segment maps for Twiss)
COSY/jobs/               runnable .fox: mapping, Twiss, chromaticity, spin tunes
COSY/run/                launchers: pre_run_cosy.bat, run_cosy.bat, run_cosy.py
COSY/src/                COSY working directory: cosy.fox, utilities.fox, elements.fox, header.fox
COSY/dat/<stem>/         output of jobs (not tracked)
```

`INCLUDE 'name'` in COSY loads a module compiled by `SAVE 'name'` in the current directory,
so `cosy.exe` always runs with cwd = `COSY/src`; structures and jobs are passed by relative path.
Jobs write output to `../dat/<stem>/`.

`COSY/src/cosy.fox` is the COSY Infinity 10.1 core (MSU license, not distributed) and is not tracked:
copy it into `COSY/src/` locally. `cosy.exe` is taken from `COSY/src/` or `PATH`.

## Run

```bat
cd COSY\run
pre_run_cosy.bat                       :: once: cosy, utilities, elements, header
run_cosy.bat lattice magnetic_2        :: compile structure (lattice + maps)
run_cosy.bat mapping                   :: run COSY\jobs\mapping.fox
run_cosy.bat lattice magnetic_2 mapping
```

Python equivalent:

```bat
python COSY\run\run_cosy.py --pre
python COSY\run\run_cosy.py --lattice magnetic_2 mapping
python COSY\run\run_cosy.py --lattice magnetic_2 --twiss magnetic_2
```

The structure used by a job is selected by its `INCLUDE` line (and `structure :=` for the output folder);
Twiss initial values and `NUM_ELE` come from the `TWISS SETUP` block in the header of `<stem>.fox`.
