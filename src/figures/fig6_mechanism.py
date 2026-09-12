#!/usr/bin/env python3
"""Figure 6: path-class fractions and Delta f_deloc vs barrier."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from _style import DOUBLE, OKABE, apply
apply()


def make_figure(results_dir: str, out_dir: str):
    with open(os.path.join(results_dir, "adw_landscapes.json")) as f:
        d = json.load(f)
    V, fH, fD, df, kie = [], [], [], [], []
    chiH, chiD = [], []
    for e in d["by_barrier"]:
        V.append(e["V0_over_kT"])
        fH.append(e["isotopes"]["H"]["f_deloc"])
        fD.append(e["isotopes"]["D"]["f_deloc"])
        df.append(e["Delta_f_deloc"])
        kie.append(e["KIE_cTST"])
        chiH.append(e["isotopes"]["H"]["chi_mean_ds"])
        chiD.append(e["isotopes"]["D"]["chi_mean_ds"])
    V = np.array(V)

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.55))
    ax = axes[0]
    ax.plot(V, fH, "o-", color=OKABE["vermillion"], ms=5, label=r"$f_{\mathrm{deloc}}$, H")
    ax.plot(V, fD, "s-", color=OKABE["blue"], ms=5, label=r"$f_{\mathrm{deloc}}$, D")
    ax.set_xlabel(r"$V_0/k_{\mathrm{B}}T$")
    ax.set_ylabel(r"delocalized fraction at $\bar\xi_p=0$")
    ax.legend(frameon=False)
    ax.set_ylim(-0.02, max(0.5, max(fH + fD) + 0.08))
    ax.text(0.03, 0.92, "(a)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    ax = axes[1]
    ax.plot(V, chiH, "o-", color=OKABE["vermillion"], ms=5, label=r"H")
    ax.plot(V, chiD, "s-", color=OKABE["blue"], ms=5, label=r"D")
    ax.set_xlabel(r"$V_0/k_{\mathrm{B}}T$")
    ax.set_ylabel(r"$\langle\chi_p\rangle_{\ddagger}$ (Å$^2$)")
    ax.legend(frameon=False)
    ax.text(0.03, 0.92, "(b)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig6_mechanism.pdf"))
    fig.savefig(os.path.join(out_dir, "fig6_mechanism.png"))
    plt.close(fig)
    print("Saved fig6_mechanism")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
