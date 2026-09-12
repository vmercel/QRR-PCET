#!/usr/bin/env python3
"""QRR-PCET paper pipeline: experiments, figures, LaTeX."""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")
SRC_EXP = os.path.join(ROOT, "src", "experiments")
SRC_FIGS = os.path.join(ROOT, "src", "figures")
MSCR_DIR = os.path.join(ROOT, "manuscript")
SUPP_DIR = os.path.join(ROOT, "supplementary")

CANONICAL_FIGURES = {
    "fig1_schematic.pdf", "fig1_schematic.png",
    "fig2_harmonic.pdf", "fig2_harmonic.png",
    "fig4_eckart.pdf", "fig4_eckart.png",
    "fig5_landscape.pdf", "fig5_landscape.png",
}

EXPERIMENT_SCRIPTS = [
    ("run_eckart.py", ""),
    ("run_eckart_pimc.py", "--n_steps {n_pimc} --n_equil {n_equil_pimc}"),
    ("run_harmonic.py", "--n_steps {n_steps} --n_equil {n_equil}"),
    ("run_adw.py", "--n_steps {n_adw} --n_equil {n_equil_adw}"),
    ("run_2d.py", "--n_steps {n_2d} --n_equil {n_equil_2d}"),
]

FIGURE_SCRIPTS = [
    "fig1_schematic.py",
    "fig2_harmonic.py",
    "fig4_eckart.py",
    "fig5_landscape.py",
]


def run_cmd(cmd, cwd=None, label=""):
    label = label or " ".join(cmd[:2])
    print(f"  > {label} ...", flush=True)
    t0 = time.time()
    result = subprocess.run(cmd, cwd=cwd)
    dt = time.time() - t0
    if result.returncode != 0:
        print(f"  FAILED ({dt:.1f}s): {label}")
        return False
    print(f"  OK ({dt:.1f}s)")
    return True


def purge_stale_figures():
    if not os.path.isdir(FIGURES):
        return
    for name in list(os.listdir(FIGURES)):
        fpath = os.path.join(FIGURES, name)
        if os.path.isfile(fpath) and name not in CANONICAL_FIGURES:
            os.remove(fpath)
            print(f"  removed stale figure: {name}")


def stage_experiments(n_steps, n_equil):
    print("\n== Stage 1: Experiments ==")
    os.makedirs(RESULTS, exist_ok=True)
    ok = True
    n_adw = max(n_steps // 2, 20_000)
    n_2d = max(n_steps, 40_000)
    for script, tmpl in EXPERIMENT_SCRIPTS:
        args_str = tmpl.format(
            n_steps=n_steps, n_equil=n_equil,
            n_adw=n_adw, n_equil_adw=max(n_equil // 2, 4000),
            n_2d=n_2d, n_equil_2d=max(n_equil // 2, 8000),
            n_pimc=min(n_steps, 40_000), n_equil_pimc=min(n_equil, 8_000),
        )
        cmd = [sys.executable, os.path.join(SRC_EXP, script)]
        if args_str.strip():
            cmd += args_str.split()
        ok = run_cmd(cmd, cwd=ROOT, label=script) and ok
    return ok


def stage_figures():
    print("\n== Stage 2: Figures ==")
    os.makedirs(FIGURES, exist_ok=True)
    ok = True
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    for script in FIGURE_SCRIPTS:
        cmd = [sys.executable, os.path.join(SRC_FIGS, script)]
        # figures expect cwd = ROOT so results/ and figures/ resolve
        print(f"  > {script} ...", flush=True)
        t0 = time.time()
        r = subprocess.run(cmd, cwd=ROOT, env=env)
        dt = time.time() - t0
        if r.returncode != 0:
            print(f"  FAILED ({dt:.1f}s): {script}")
            ok = False
        else:
            print(f"  OK ({dt:.1f}s)")
    purge_stale_figures()
    missing = [f for f in CANONICAL_FIGURES if not os.path.isfile(os.path.join(FIGURES, f))]
    if missing:
        print(f"  missing figures: {missing}")
        ok = False
    else:
        print(f"  all {len(CANONICAL_FIGURES)} canonical figures present")
    return ok


def stage_compile():
    print("\n== Stage 3: LaTeX ==")
    ok = True
    for tex_dir, tex_name in [(MSCR_DIR, "manuscript"), (SUPP_DIR, "supplementary")]:
        tex_file = tex_name + ".tex"
        if not os.path.isfile(os.path.join(tex_dir, tex_file)):
            print(f"  skip {tex_file}")
            continue
        for i, cmd in enumerate([
            ["pdflatex", "-interaction=nonstopmode", tex_file],
            ["bibtex", tex_name],
            ["pdflatex", "-interaction=nonstopmode", tex_file],
            ["pdflatex", "-interaction=nonstopmode", tex_file],
        ]):
            if not run_cmd(cmd, cwd=tex_dir, label=f"{tex_name}: {' '.join(cmd[:2])}"):
                ok = False
                if i == 0:
                    break
        pdf = os.path.join(tex_dir, tex_name + ".pdf")
        if os.path.isfile(pdf):
            print(f"  PDF {pdf} ({os.path.getsize(pdf)//1024} KB)")
        else:
            ok = False
    return ok


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_steps", type=int, default=80_000)
    p.add_argument("--n_equil", type=int, default=15_000)
    p.add_argument("--skip_experiments", action="store_true")
    p.add_argument("--skip_compile", action="store_true")
    p.add_argument("--only", choices=["experiments", "figures", "compile"], default=None)
    args = p.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(FIGURES, exist_ok=True)

    if args.only == "experiments":
        stage_experiments(args.n_steps, args.n_equil)
    elif args.only == "figures":
        stage_figures()
    elif args.only == "compile":
        stage_compile()
    else:
        ok = True
        if not args.skip_experiments:
            ok &= stage_experiments(args.n_steps, args.n_equil)
        ok &= stage_figures()
        if not args.skip_compile:
            ok &= stage_compile()
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
