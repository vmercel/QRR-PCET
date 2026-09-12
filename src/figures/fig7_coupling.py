#!/usr/bin/env python3
"""Figure 7: 2D proton-solvent factorization diagnostic vs coupling."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from _style import DOUBLE, OKABE, apply
apply()


def make_figure(results_dir: str, out_dir: str):
    with open(os.path.join(results_dir, "coupling_2d.json")) as f:
        d = json.load(f)
    kc = np.array([e["kappa_c"] for e in d["scan"]])
    corr = np.array([e["corr_chi_qs"] for e in d["scan"]])
    fd = np.array([e["f_deloc_ds"] for e in d["scan"]])

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.55))
    ax = axes[0]
    ax.plot(kc, corr, "o-", color=OKABE["blue"], ms=6)
    ax.axhline(0.0, color="0.5", ls=":", lw=0.8)
    ax.set_xlabel(r"$\kappa_c$ (kJ mol$^{-1}$ Å$^{-2}$)")
    ax.set_ylabel(r"$\mathrm{corr}(\chi_p,q_s)\,|\,_{\ddagger}$")
    ax.text(0.03, 0.92, "(a)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    ax = axes[1]
    ax.plot(kc, fd, "s-", color=OKABE["vermillion"], ms=6)
    ax.set_xlabel(r"$\kappa_c$ (kJ mol$^{-1}$ Å$^{-2}$)")
    ax.set_ylabel(r"$f_{\mathrm{deloc}}$ at dividing surface")
    ax.text(0.03, 0.92, "(b)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig7_coupling.pdf"))
    fig.savefig(os.path.join(out_dir, "fig7_coupling.png"))
    plt.close(fig)
    print("Saved fig7_coupling")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
