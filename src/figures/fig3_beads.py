#!/usr/bin/env python3
"""Figure 3: bead-number convergence of chi_p vs exact finite-P formula."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from _style import SINGLE, OKABE, apply
apply()


def make_figure(results_dir: str, out_dir: str):
    with open(os.path.join(results_dir, "harmonic_benchmark.json")) as f:
        d = json.load(f)
    bc = d["bead_convergence"]
    P = np.array(bc["P_values"])
    finH = np.array(bc["finite_P_H"])
    pimH = np.array(bc["pimc_H"])
    sH = np.array(bc["pimc_H_std"])
    cont = bc["continuum_H"]

    fig, ax = plt.subplots(figsize=(SINGLE, 2.55))
    ax.axhline(cont, color=OKABE["black"], ls="--", lw=1.1, label=r"$P\to\infty$")
    ax.plot(P, finH, color=OKABE["blue"], lw=1.6, marker="s", ms=5,
            label="exact finite $P$")
    ax.errorbar(P, pimH, yerr=sH, fmt="o", color=OKABE["vermillion"], ms=5,
                capsize=2, label="PIMC")
    ax.set_xscale("log", base=2)
    ax.set_xticks(P)
    ax.set_xticklabels([str(p) for p in P])
    ax.set_xlabel(r"number of beads $P$")
    ax.set_ylabel(r"$\langle\chi_p\rangle_H$ (Å$^2$) at 300 K")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig3_beads.pdf"))
    fig.savefig(os.path.join(out_dir, "fig3_beads.png"))
    plt.close(fig)
    print("Saved fig3_beads")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
