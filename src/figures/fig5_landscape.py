#!/usr/bin/env python3
"""Figure 5: ADW F(xi, chi) landscapes for H at two barrier heights."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from _style import DOUBLE, apply
apply()


def make_figure(results_dir: str, out_dir: str):
    with open(os.path.join(results_dir, "adw_landscapes.json")) as f:
        d = json.load(f)
    # pick V0 = 4 and 8 kT if present
    by = {e["V0_over_kT"]: e for e in d["by_barrier"]}
    keys = sorted(by.keys())
    pick = [keys[0], keys[-1]] if len(keys) >= 2 else keys

    fig, axes = plt.subplots(1, len(pick), figsize=(DOUBLE, 2.65), sharey=True)
    if len(pick) == 1:
        axes = [axes]
    for ax, vkt in zip(axes, pick):
        iso = by[vkt]["isotopes"]["H"]
        xi = np.array(iso["xi_cent"])
        chi = np.array(iso["chi_cent"])
        # histogram2d(xi, chi) -> (n_xi, n_chi); contourf needs (n_chi, n_xi)
        Z = np.array(iso["F2d"]).T.astype(float)
        finite = Z[np.isfinite(Z)]
        cap = np.percentile(finite, 90) if finite.size else 40.0
        Z = np.clip(Z, 0.0, cap)
        # hide nearly unsampled (very high F) as white
        Z = np.ma.masked_where(Z >= 0.98 * cap, Z)
        cf = ax.contourf(xi, chi, Z, levels=10, cmap="cividis")
        ax.set_ylim(0.0, 0.045)
        ax.axvline(0.0, color="w", ls=":", lw=0.7, alpha=0.7)
        ax.set_xlabel(r"$\bar\xi_p$ (Å)")
        ax.set_title(rf"$V_0={vkt:.0f}\,k_{{\mathrm{{B}}}}T$ (H)", fontsize=9)
        fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.03).set_label("kJ/mol", fontsize=7)
    axes[0].set_ylabel(r"$\chi_p$ (Å$^2$)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig5_landscape.pdf"))
    fig.savefig(os.path.join(out_dir, "fig5_landscape.png"))
    plt.close(fig)
    print("Saved fig5_landscape")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
