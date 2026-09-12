#!/usr/bin/env python3
"""Figure 5: ADW centroid PMF, dividing-surface p(chi), and isotope contrast."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from _style import DOUBLE, OKABE, apply
apply()

KB_KJMOL = 0.008314462618
T_K = 300.0
CHI_C = 0.010


def p_chi_at_ds(iso):
    xi = np.array(iso["xi_cent"])
    chi = np.array(iso["chi_cent"])
    F2 = np.array(iso["F2d"], dtype=float)
    i_div = int(np.argmin(np.abs(xi)))
    F = F2[i_div, :]
    # empty bins were set to a huge floor; drop them
    cap = np.percentile(F[np.isfinite(F)], 80)
    mask = np.isfinite(F) & (F < cap + 1e-9)
    if mask.sum() < 4:
        mask = np.isfinite(F)
    kT = KB_KJMOL * T_K
    w = np.zeros_like(F)
    w[mask] = np.exp(-(F[mask] - np.min(F[mask])) / kT)
    # chi bin width
    dchi = np.median(np.diff(chi)) if len(chi) > 1 else 1.0
    Z = np.sum(w) * dchi
    pdf = w / Z if Z > 0 else w
    return chi, pdf


def make_figure(results_dir: str, out_dir: str):
    with open(os.path.join(results_dir, "adw_landscapes.json")) as f:
        d = json.load(f)
    by = {e["V0_over_kT"]: e for e in d["by_barrier"]}
    hi = by[8.0]

    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 2.50))

    ax = axes[0]
    for iso, col, ls in (("H", OKABE["vermillion"], "-"), ("D", OKABE["blue"], "--")):
        z = hi["isotopes"][iso]
        xi = np.array(z["xi_cent"])
        F = np.array(z["F_xi"])
        ax.plot(xi, F, color=col, ls=ls, lw=1.6, label=iso)
    ax.axvline(0.0, color="0.6", ls=":", lw=0.8)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-0.5, 22)
    ax.set_xlabel(r"$\bar\xi_p$ (Å)")
    ax.set_ylabel(r"$F(\bar\xi_p)$ (kJ/mol)")
    ax.legend(frameon=False, loc="upper right")
    ax.set_title(r"$V_0=8\,k_{\mathrm{B}}T$", fontsize=8, pad=3)
    ax.text(0.04, 0.92, "(a)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    ax = axes[1]
    for iso, col, ls in (("H", OKABE["vermillion"], "-"), ("D", OKABE["blue"], "--")):
        chi, pdf = p_chi_at_ds(hi["isotopes"][iso])
        ax.plot(chi, pdf, color=col, ls=ls, lw=1.6, label=iso)
    ax.axvline(2.0 * CHI_C, color="0.35", ls=":", lw=1.0)
    ax.set_xlim(0.0, 0.040)
    ax.set_xlabel(r"$\chi_p$ (Å$^2$)")
    ax.set_ylabel(r"$p^\ddagger(\chi_p)$ (Å$^{-2}$)")
    ax.legend(frameon=False)
    ax.set_title(r"$p^\ddagger$ at $\bar\xi_p=0$", fontsize=8, pad=3)
    ax.text(2.0 * CHI_C + 0.0015, 0.92, r"$2\chi_c$", transform=ax.get_xaxis_transform(),
            fontsize=7, color="0.3")
    ax.text(0.04, 0.92, "(b)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    ax = axes[2]
    V, fH, fD = [], [], []
    for e in d["by_barrier"]:
        V.append(e["V0_over_kT"])
        fH.append(e["isotopes"]["H"]["f_deloc"])
        fD.append(e["isotopes"]["D"]["f_deloc"])
    ax.plot(V, fH, "o-", color=OKABE["vermillion"], ms=5, label="H")
    ax.plot(V, fD, "s--", color=OKABE["blue"], ms=5, label="D")
    ax.set_xlabel(r"$V_0/k_{\mathrm{B}}T$")
    ax.set_ylabel(r"$f_{\mathrm{deloc}}$")
    ax.set_ylim(-0.02, 0.32)
    ax.legend(frameon=False)
    ax.text(0.04, 0.92, "(c)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig5_landscape.pdf"))
    fig.savefig(os.path.join(out_dir, "fig5_landscape.png"))
    plt.close(fig)
    print("Saved fig5_landscape")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
