#!/usr/bin/env python3
"""Figure 4: Eckart exact vs Wigner vs centroid QTST (PIMC) vs classical.

Unreliable PIMC points (reliable=false in eckart_pimc.json) are omitted
from the rate and KIE panels. Path-class fractions are shown for every
barrier; D at 12 kT is consistent with the other D points even though
its centroid barrier is not.
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from _style import DOUBLE, OKABE, apply
apply()


def _reliable(iso):
    return bool(iso.get("reliable", True))


def make_figure(results_dir: str, out_dir: str):
    pimc_path = os.path.join(results_dir, "eckart_pimc.json")
    ana_path = os.path.join(results_dir, "eckart_rates.json")
    with open(ana_path) as f:
        ana = json.load(f)
    V_ana = np.array(ana["V0_over_kT"])
    kexH = np.array(ana["kappa_exact_H"])
    kexD = np.array(ana["kappa_exact_D"])
    kwH = np.array(ana["kappa_wigner_H"])
    kwD = np.array(ana["kappa_wigner_D"])
    kie_ex = np.array(ana["KIE_exact"])
    kie_w = np.array(ana["KIE_wigner"])
    kie_cl = np.array(ana["KIE_tst"])

    pc = None
    if os.path.isfile(pimc_path):
        with open(pimc_path) as f:
            pc = json.load(f)

    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 2.60))
    legkw = dict(frameon=False, fontsize=6.5, handlelength=1.4,
                 borderpad=0.2, labelspacing=0.28, handletextpad=0.4)

    ax = axes[0]
    ax.axhline(1.0, color="0.45", ls=":", lw=0.8, label="TST")
    ax.plot(V_ana, kexH, "o-", color=OKABE["vermillion"], ms=5, label="exact H")
    ax.plot(V_ana, kexD, "s-", color=OKABE["blue"], ms=5, label="exact D")
    ax.plot(V_ana, kwH, "^--", color=OKABE["orange"], ms=4, label="Wigner H")
    ax.plot(V_ana, kwD, "v--", color=OKABE["sky"], ms=4, label="Wigner D")
    if pc is not None:
        VpH, cH, VpD, cD = [], [], [], []
        for e in pc["by_barrier"]:
            v = e["V0_over_kT"]
            if _reliable(e["isotopes"]["H"]):
                VpH.append(v)
                cH.append(e["isotopes"]["H"]["kappa_ctst"])
            if _reliable(e["isotopes"]["D"]):
                VpD.append(v)
                cD.append(e["isotopes"]["D"]["kappa_ctst"])
        ax.plot(VpH, cH, "D", color=OKABE["vermillion"], ms=7, zorder=6,
                label="cTST H")
        ax.plot(VpD, cD, "D", color=OKABE["blue"], ms=7, zorder=6, mfc="white",
                label="cTST D")
    ax.set_xlabel(r"$V_0/k_{\mathrm{B}}T$")
    ax.set_ylabel(r"$\kappa = k/k_{\mathrm{TST}}$")
    ax.set_ylim(0.95, 2.85)
    ax.legend(**legkw, loc="upper left", ncol=1)
    ax.text(0.72, 0.92, "(a)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    ax = axes[1]
    ax.plot(V_ana, kie_ex, "o-", color=OKABE["black"], ms=5, label="exact")
    ax.plot(V_ana, kie_w, "s--", color=OKABE["green"], ms=5, label="Wigner")
    ax.plot(V_ana, kie_cl, "^:", color=OKABE["blue"], ms=5, label=r"classical $\sqrt{2}$")
    if pc is not None:
        Vp, kie_c = [], []
        for e in pc["by_barrier"]:
            if _reliable(e["isotopes"]["H"]) and _reliable(e["isotopes"]["D"]):
                Vp.append(e["V0_over_kT"])
                kie_c.append(e["KIE_ctst"])
        ax.plot(Vp, kie_c, "D", color=OKABE["vermillion"], ms=7, label="cTST")
    ax.set_xlabel(r"$V_0/k_{\mathrm{B}}T$")
    ax.set_ylabel(r"$k_{\mathrm{H}}/k_{\mathrm{D}}$")
    ax.set_ylim(1.32, 2.40)
    ax.legend(**legkw, loc="upper left")
    ax.text(0.86, 0.92, "(b)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    ax = axes[2]
    if pc is not None:
        Vp, fH, fD = [], [], []
        for e in pc["by_barrier"]:
            Vp.append(e["V0_over_kT"])
            fH.append(e["isotopes"]["H"]["f_deloc"])
            fD.append(e["isotopes"]["D"]["f_deloc"])
        ax.plot(Vp, fH, "o-", color=OKABE["vermillion"], ms=5, label="H")
        ax.plot(Vp, fD, "s--", color=OKABE["blue"], ms=5, label="D")
    ax.set_xlabel(r"$V_0/k_{\mathrm{B}}T$")
    ax.set_ylabel(r"$f_{\mathrm{deloc}}$ at barrier")
    ax.legend(frameon=False, fontsize=7.0)
    ax.set_ylim(-0.02, 0.45)
    ax.text(0.04, 0.92, "(c)", transform=ax.transAxes, fontsize=9, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig4_eckart.pdf"))
    fig.savefig(os.path.join(out_dir, "fig4_eckart.png"))
    plt.close(fig)
    print("Saved fig4_eckart")


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    make_figure("results", "figures")
