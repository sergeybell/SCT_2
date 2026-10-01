#!/usr/bin/env python3
"""Launch COSY Infinity with cwd=COSY/src (compiled modules live there)."""
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
COSY = REPO / "COSY"
COSY_SRC = COSY / "src"
JOBS = COSY / "jobs"
STRUCTURES = COSY / "structures"
DAT_ROOT = COSY / "dat"

PRE_RUN_FOX = ["cosy.fox", "utilities.fox", "elements.fox", "header.fox"]
TWISS_OUTPUT_FILES = ["BETAX", "BETAY", "DISPX", "DISPY"]


def cosy_exe() -> str:
    local = COSY_SRC / "cosy.exe"
    if local.is_file():
        return str(local)
    found = shutil.which("cosy.exe") or shutil.which("cosy")
    if found is None:
        raise FileNotFoundError(f"cosy.exe not found in {COSY_SRC} or PATH")
    return found


def resolve_job(name: str) -> Path:
    """Job name (mapping), file name (mapping.fox) or path to any .fox."""
    p = Path(name)
    for cand in (p, JOBS / p, JOBS / p.with_suffix(".fox")):
        if cand.is_file():
            return cand.resolve()
    raise FileNotFoundError(f"Job not found: {name}")


def structure_fox(stem: str, *, maps: bool = False) -> Path:
    p = STRUCTURES / stem / (f"{stem}_maps.fox" if maps else f"{stem}.fox")
    if not p.is_file():
        raise FileNotFoundError(f"Structure file not found: {p}")
    return p


def run_cosy(fox: Path, *, verbose: bool = True) -> None:
    fox = fox.resolve()
    if fox.is_relative_to(COSY_SRC):
        arg = fox.relative_to(COSY_SRC).as_posix()
    elif fox.is_relative_to(COSY):
        arg = "../" + fox.relative_to(COSY).as_posix()
    else:
        arg = fox.as_posix()
    if verbose:
        print("========================================")
        print(f"RUNNING: {arg}")
        print("========================================")
    (COSY_SRC / "tmp").mkdir(exist_ok=True)
    subprocess.run([cosy_exe(), arg], cwd=str(COSY_SRC), check=True)


def run_pre() -> None:
    for name in PRE_RUN_FOX:
        run_cosy(COSY_SRC / name)


def compile_structure(stem: str) -> None:
    run_cosy(structure_fox(stem))
    run_cosy(structure_fox(stem, maps=True))
    (DAT_ROOT / stem).mkdir(parents=True, exist_ok=True)


def collect_twiss(stem: str) -> Path:
    dest = DAT_ROOT / stem
    dest.mkdir(parents=True, exist_ok=True)
    for name in TWISS_OUTPUT_FILES:
        src = COSY_SRC / name
        if src.is_file():
            shutil.move(str(src), str(dest / name))
    return dest


def main() -> int:
    ap = argparse.ArgumentParser(description="Run COSY Infinity (cwd = COSY/src)")
    ap.add_argument("--pre", action="store_true", help="Compile base modules: cosy, utilities, elements, header")
    ap.add_argument("--lattice", metavar="STEM", help="Compile COSY/structures/STEM (lattice + maps)")
    ap.add_argument("--twiss", metavar="STEM", help="Run jobs/Twiss.fox and move BETAX/BETAY/DISPX to COSY/dat/STEM")
    ap.add_argument("jobs", nargs="*", help="Job names from COSY/jobs (e.g. mapping) or paths to .fox")
    args = ap.parse_args()

    if not (args.pre or args.lattice or args.twiss or args.jobs):
        ap.print_help()
        return 1

    if args.pre:
        run_pre()
    if args.lattice:
        compile_structure(args.lattice)
    for job in args.jobs:
        run_cosy(resolve_job(job))
    if args.twiss:
        run_cosy(JOBS / "Twiss.fox")
        print(f"OK: Twiss output in {collect_twiss(args.twiss)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
