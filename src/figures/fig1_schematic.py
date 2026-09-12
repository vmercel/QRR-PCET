#!/usr/bin/env python3
"""Figure 1: centroid PMF vs path-spread-resolved free energy (schematic)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
sys.path.insert(0, os.path.dirname(__file__))
from _style import DOUBLE, OKABE, apply

apply()


def pmf(xi, V0=8.0):
    return V0 * (xi**2 - 1.0) ** 2


def F2(xi, chi, V0=8.0):
    # illustrative: delocalized paths see a lower ridge
    return V0 * (xi**2 - 1.0) ** 2 * np.exp(-chi / 0.025) + 18.0 * (chi - 0.012) ** 2 / 0.012


def make_figure(results_dir: str, out_dir: str):
    fig = plt.figure(figsize=(DOUBLE, 2.55))
    gs = GridSpec(1, 2, figure=fig, wspace=0.32)

    xi = np.linspace(-1.7, 1.7, 400)
    ax0 = fig.add_subplot(gs[0])
    ax0.plot(xi, pmf(xi), color=OKABE["blue"], lw=1.8)
    ax0.axvline(0.0, color="0.6", ls=":", lw=0.8)
    i_bar = np.argmin(np.abs(xi))
    ax0.plot(0.0, pmf(0.0), "o", color=OKABE["vermillion"], ms=6, zorder=5)
    ax0.annotate(r"$\Delta F^\ddagger$", xy=(0.0, pmf(0.0)), xytext=(0.55, 6.2),
                 fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8, color="0.3"),
                 color="0.2")
    ax0.set_xlabel(r"centroid $\bar\xi_p$ (Å)")
    ax0.set_ylabel(r"$F(\bar\xi_p)$ (kJ/mol)")
    ax0.set_xlim(-1.7, 1.7)
    ax0.set_ylim(-0.5, 10)
    ax0.set_title("Centroid PMF", fontsize=9, pad=4)
    ax0.text(0.03, 0.92, "(a)", transform=ax0.transAxes, fontsize=9, fontweight="bold")

    ax1 = fig.add_subplot(gs[1])
    xi2 = np.linspace(-1.6, 1.6, 160)
    chi2 = np.linspace(0.0, 0.045, 120)
    XI, CHI = np.meshgrid(xi2, chi2)
    Z = F2(XI, CHI)
    Z -= Z.min()
    cf = ax1.contourf(xi2, chi2, Z, levels=14, cmap="cividis")
    ax1.contour(xi2, chi2, Z, levels=7, colors="w", linewidths=0.3, alpha=0.4)
    # MFEP sketch
    xpath = np.linspace(-1.1, 1.1, 80)
    chipath = 0.012 + 0.010 * np.exp(-3.5 * xpath**2)
    ax1.plot(xpath, chipath, color="white", ls="--", lw=1.2)
    ax1.plot(0.0, 0.022, marker="*", color=OKABE["vermillion"], ms=9, zorder=5)
    ax1.axvline(0.0, color="w", ls=":", lw=0.7, alpha=0.7)
    ax1.set_xlabel(r"centroid $\bar\xi_p$ (Å)")
    ax1.set_ylabel(r"path spread $\chi_p$ (Å$^2$)")
    ax1.set_title("Path-spread-resolved $F$", fontsize=9, pad=4)
    ax1.text(0.03, 0.92, "(b)", transform=ax1.transAxes, fontsize=9,
             fontweight="bold", color="white")
    cbar = fig.colorbar(cf, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label("kJ/mol", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)

    fig.savefig(os.path.join(out_dir, "fig1_schematic.pdf"))
    fig.savefig(os.path.join(out_dir, "fig1_schematic.png"))
    plt.close(fig)
    print("Saved fig1_schematic")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
