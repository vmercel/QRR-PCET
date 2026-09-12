#!/usr/bin/env python3
"""Figure 2: harmonic <chi_p>(T) and isotope ratio."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from _style import DOUBLE, OKABE, apply
apply()


def make_figure(results_dir: str, out_dir: str):
    with open(os.path.join(results_dir, "harmonic_benchmark.json")) as f:
        d = json.load(f)
    T = np.array(d["temperatures_K"])
    cH, cD = np.array(d["continuum"]["H"]), np.array(d["continuum"]["D"])
    pH, sH = np.array(d["pimc"]["H_mean"]), np.array(d["pimc"]["H_std"])
    pD, sD = np.array(d["pimc"]["D_mean"]), np.array(d["pimc"]["D_std"])
    rC = np.array(d["ratio_continuum"])
    rP = np.array(d["ratio_pimc"])

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.55))

    ax = axes[0]
    ax.plot(T, cH, color=OKABE["vermillion"], lw=1.6, label="H, continuum")
    ax.plot(T, cD, color=OKABE["blue"], lw=1.6, label="D, continuum")
    ax.errorbar(T, pH, yerr=sH, fmt="o", color=OKABE["vermillion"], ms=4,
                capsize=2, zorder=5, label="H, PIMC $P=64$")
    ax.errorbar(T, pD, yerr=sD, fmt="s", color=OKABE["blue"], ms=4,
                capsize=2, zorder=5, label="D, PIMC $P=32$")
    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$\langle\chi_p\rangle$ (Å$^2$)")
    ax.legend(frameon=False, loc="upper right")
    ax.set_xlim(80, 740)
    ax.text(0.03, 0.92, "(a)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    ax = axes[1]
    ax.plot(T, rC, color=OKABE["black"], lw=1.6, label="continuum")
    ax.plot(T, rP, "o", color=OKABE["green"], ms=5, label="PIMC")
    ax.axhline(np.sqrt(2.0), color=OKABE["vermillion"], ls="--", lw=1.0,
               label=r"$\sqrt{m_D/m_H}=\sqrt{2}$")
    ax.axhline(2.0, color=OKABE["blue"], ls=":", lw=1.0,
               label=r"$m_D/m_H=2$")
    ax.set_xlabel(r"$T$ (K)")
    ax.set_ylabel(r"$\langle\chi_p\rangle_H/\langle\chi_p\rangle_D$")
    ax.legend(frameon=False, loc="lower right")
    ax.set_xlim(80, 740)
    ax.set_ylim(1.35, 2.15)
    ax.text(0.03, 0.92, "(b)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig2_harmonic.pdf"))
    fig.savefig(os.path.join(out_dir, "fig2_harmonic.png"))
    plt.close(fig)
    print("Saved fig2_harmonic")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
