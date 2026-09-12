#!/usr/bin/env python3
"""Figure 1: schematic centroid PMF vs path-spread-resolved F. Not computed data."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
sys.path.insert(0, os.path.dirname(__file__))
from _style import DOUBLE, OKABE, apply
apply()


def pmf(xi, V0=8.0):
    return V0 * (xi**2 - 1.0) ** 2


def F2(xi, chi):
    # Pedagogical surface: double well in xi, harmonic confinement in chi
    # whose minimum rises at the barrier, so delocalized paths sit higher
    # in the wells but lower on the ridge than compact ones.
    well = 8.0 * (xi**2 - 1.0) ** 2
    chi0 = 0.010 + 0.012 * np.exp(-4.0 * xi**2)
    return well * np.exp(-chi / 0.030) + 400.0 * (chi - chi0) ** 2


def make_figure(results_dir: str, out_dir: str):
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.65))

    xi = np.linspace(-1.65, 1.65, 400)
    ax = axes[0]
    ax.plot(xi, pmf(xi), color=OKABE["blue"], lw=1.9)
    ax.axvline(0.0, color="0.55", ls=":", lw=0.8)
    ax.plot(0.0, pmf(0.0), "o", color=OKABE["vermillion"], ms=6, zorder=5)
    ax.annotate(r"$\Delta F^\ddagger$", xy=(0.0, pmf(0.0)),
                xytext=(0.45, 5.8), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.25"))
    ax.set_xlabel(r"centroid $\bar\xi_p$ (Å)")
    ax.set_ylabel(r"$F(\bar\xi_p)$ (kJ/mol)")
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-0.4, 10.2)
    ax.set_title("Centroid PMF (schematic)", fontsize=9, pad=4)
    ax.text(0.03, 0.92, "(a)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    ax = axes[1]
    xi2 = np.linspace(-1.5, 1.5, 200)
    chi2 = np.linspace(0.0, 0.040, 160)
    XI, CHI = np.meshgrid(xi2, chi2)
    Z = F2(XI, CHI)
    Z -= Z.min()
    cf = ax.contourf(xi2, chi2, Z, levels=16, cmap="cividis")
    ax.contour(xi2, chi2, Z, levels=8, colors="w", linewidths=0.25, alpha=0.45)
    ax.axvline(0.0, color="w", ls=":", lw=0.8, alpha=0.85)
    ax.plot(0.0, 0.010, "o", color="white", ms=5, zorder=6)
    ax.plot(0.0, 0.024, "*", color=OKABE["vermillion"], ms=10, zorder=6)
    ax.annotate("compact", xy=(0.08, 0.010), fontsize=7, color="white")
    ax.annotate("delocalized", xy=(0.08, 0.026), fontsize=7,
                color=OKABE["vermillion"])
    ax.set_xlabel(r"centroid $\bar\xi_p$ (Å)")
    ax.set_ylabel(r"path spread $\chi_p$ (Å$^2$)")
    ax.set_title(r"$F(\bar\xi_p,\chi_p)$ (schematic)", fontsize=9, pad=4)
    ax.text(0.03, 0.92, "(b)", transform=ax.transAxes, fontsize=9,
            fontweight="bold", color="white")
    cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("kJ/mol", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig1_schematic.pdf"))
    fig.savefig(os.path.join(out_dir, "fig1_schematic.png"))
    plt.close(fig)
    print("Saved fig1_schematic")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
