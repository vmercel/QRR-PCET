#!/usr/bin/env python3
"""Figure 4: Eckart exact kappa and KIE vs barrier height."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from _style import DOUBLE, OKABE, apply
apply()


def make_figure(results_dir: str, out_dir: str):
    with open(os.path.join(results_dir, "eckart_rates.json")) as f:
        d = json.load(f)
    V = np.array(d["V0_over_kT"])
    kH = np.array(d["kappa_exact_H"])
    kD = np.array(d["kappa_exact_D"])
    wH = np.array(d["kappa_wigner_H"])
    wD = np.array(d["kappa_wigner_D"])
    kie = np.array(d["KIE_exact"])
    kieT = np.array(d["KIE_tst"])
    kieW = np.array(d["KIE_wigner"])

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.55))
    ax = axes[0]
    ax.plot(V, kH, "o-", color=OKABE["vermillion"], ms=5, label="exact, H")
    ax.plot(V, kD, "s-", color=OKABE["blue"], ms=5, label="exact, D")
    ax.plot(V, wH, "o--", color=OKABE["orange"], ms=4, label="Wigner, H")
    ax.plot(V, wD, "s--", color=OKABE["sky"], ms=4, label="Wigner, D")
    ax.axhline(1.0, color="0.5", ls=":", lw=0.8)
    ax.set_xlabel(r"$V_0/k_{\mathrm{B}}T$")
    ax.set_ylabel(r"$\kappa = k/k_{\mathrm{TST}}$")
    ax.legend(frameon=False, loc="upper left")
    ax.text(0.03, 0.92, "(a)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    ax = axes[1]
    ax.plot(V, kie, "o-", color=OKABE["black"], ms=5, label="exact QM")
    ax.plot(V, kieW, "s--", color=OKABE["green"], ms=5, label="Wigner")
    ax.plot(V, kieT, "^:", color=OKABE["blue"], ms=5, label=r"classical TST ($\sqrt{2}$)")
    ax.set_xlabel(r"$V_0/k_{\mathrm{B}}T$")
    ax.set_ylabel(r"$k_{\mathrm{H}}/k_{\mathrm{D}}$")
    ax.legend(frameon=False, loc="upper left")
    ax.text(0.03, 0.92, "(b)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig4_eckart.pdf"))
    fig.savefig(os.path.join(out_dir, "fig4_eckart.png"))
    plt.close(fig)
    print("Saved fig4_eckart")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
