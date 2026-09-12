#!/usr/bin/env python3
"""
Asymmetric double-well: WHAM centroid PMF, 2D F(xi, chi), path-class fractions.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from models.qrr_core import (  # noqa: E402
    KB_KJMOL,
    adw_V,
    make_rng,
    pimc_1d,
    wham_1d,
)

T_K = 300.0
Q0 = 0.70
ALPHA = 0.025
V0_KT = [2.0, 4.0, 6.0, 8.0]
P_H, P_D = 32, 24
SEED = 42
CHI_C = 0.010  # Å^2, ~ harmonic ZPE spread at 1000 cm^-1 (not 0.10)


def run(n_steps: int = 80_000, n_equil: int = 12_000) -> dict:
    os.makedirs("results", exist_ok=True)
    rng = make_rng(SEED)
    kT = KB_KJMOL * T_K
    windows = np.linspace(-1.00, 1.00, 11)
    k_umb = 80.0  # kJ/mol/Å^2

    results = {
        "description": "ADW WHAM PMF and path-spread landscapes",
        "T_K": T_K,
        "q0_angstrom": Q0,
        "alpha_kJ_per_A": ALPHA,
        "V0_over_kT": V0_KT,
        "P_H": P_H,
        "P_D": P_D,
        "chi_c_ang2": CHI_C,
        "n_steps": n_steps,
        "n_equil": n_equil,
        "windows": windows.tolist(),
        "k_umb": k_umb,
        "seed": SEED,
        "by_barrier": [],
    }

    xi_bins = np.linspace(-1.25, 1.25, 41)
    chi_bins = np.linspace(0.0, 0.08, 25)
    xi_cent = 0.5 * (xi_bins[:-1] + xi_bins[1:])
    chi_cent = 0.5 * (chi_bins[:-1] + chi_bins[1:])
    dxi = xi_cent[1] - xi_cent[0]

    for vkt in V0_KT:
        V0 = vkt * kT
        print(f"ADW V0 = {vkt} kT")
        entry = {"V0_over_kT": vkt, "isotopes": {}}
        for iso, mass, P in (("H", 1.0, P_H), ("D", 2.0, P_D)):
            def V(q, _V0=V0):
                return adw_V(q, _V0, Q0, ALPHA)

            xi_all = []
            chi_all = []
            win_id = []
            counts = np.zeros((len(windows), len(xi_cent)))
            N_win = np.zeros(len(windows))
            U_bias = np.zeros((len(windows), len(xi_cent)))
            for i, x0 in enumerate(windows):
                U_bias[i] = 0.5 * k_umb * (xi_cent - x0) ** 2
                samp = pimc_1d(
                    V, T_K, mass, P, n_steps, n_equil, x0, rng,
                    k_umb=k_umb, q_umb=x0, centroid_step=0.03,
                )
                xi_s, chi_s = samp["xi"], samp["chi"]
                hist, _ = np.histogram(xi_s, bins=xi_bins)
                counts[i] = hist
                N_win[i] = len(xi_s)
                xi_all.append(xi_s)
                chi_all.append(chi_s)
                win_id.append(np.full(len(xi_s), i))
                print(f"    {iso} window {x0:+.2f}  <xi>={xi_s.mean():+.3f}  <chi>={chi_s.mean():.4f}")

            F, Pxi, f_i = wham_1d(counts, N_win, U_bias, T_K, dxi)
            # barrier from WHAM PMF
            i_div = int(np.argmin(np.abs(xi_cent)))
            left = xi_cent < -0.15
            i_min = int(np.argmin(F[left])) if np.any(left) else 0
            # map i_min in left-subset back
            left_idx = np.where(left)[0]
            i_well = int(left_idx[i_min]) if left_idx.size else 0
            dF = float(F[i_div] - F[i_well])

            xi_cat = np.concatenate(xi_all)
            chi_cat = np.concatenate(chi_all)
            wid = np.concatenate(win_id)
            beta = 1.0 / kT
            Vbias_s = 0.5 * k_umb * (xi_cat - windows[wid]) ** 2
            w = np.exp(beta * (Vbias_s + f_i[wid]))
            w /= w.sum()

            # 2D weighted histogram
            H2, _, _ = np.histogram2d(xi_cat, chi_cat, bins=[xi_bins, chi_bins], weights=w)
            H2 = np.maximum(H2, 1e-16)
            F2 = -kT * np.log(H2)
            F2 -= np.nanmin(F2)

            # dividing-surface path-class distribution: |xi| < 0.08 Å
            mask_ds = np.abs(xi_cat) < 0.08
            if mask_ds.sum() < 50:
                mask_ds = np.abs(xi_cat) < 0.15
            w_ds = w[mask_ds]
            chi_ds = chi_cat[mask_ds]
            w_ds = w_ds / w_ds.sum()
            f_comp = float(np.sum(w_ds[chi_ds < 0.5 * CHI_C]))
            f_del = float(np.sum(w_ds[chi_ds >= 2.0 * CHI_C]))
            f_int = 1.0 - f_comp - f_del
            chi_mean_ds = float(np.sum(w_ds * chi_ds))

            entry["isotopes"][iso] = {
                "DeltaF_kJmol": round(dF, 4),
                "DeltaF_kT": round(dF / kT, 4),
                "f_compact": round(max(f_comp, 0.0), 4),
                "f_intermediate": round(max(f_int, 0.0), 4),
                "f_deloc": round(max(f_del, 0.0), 4),
                "chi_mean_ds": round(chi_mean_ds, 6),
                "xi_cent": np.round(xi_cent, 5).tolist(),
                "F_xi": np.round(F, 5).tolist(),
                "chi_cent": np.round(chi_cent, 5).tolist(),
                "F2d": np.round(F2, 4).tolist(),
            }
            print(f"    {iso}  ΔF={dF:.3f} kJ/mol  f_deloc={f_del:.3f}  <chi>_ds={chi_mean_ds:.4f}")

        # KIE proxy from centroid TST using same well frequency scale
        dF_H = entry["isotopes"]["H"]["DeltaF_kJmol"]
        dF_D = entry["isotopes"]["D"]["DeltaF_kJmol"]
        # k ~ omega exp(-β ΔF); omega ~ 1/sqrt(m) at the well
        kie_ctst = np.sqrt(2.0) * np.exp(-(dF_H - dF_D) / kT)
        entry["KIE_cTST"] = round(float(kie_ctst), 4)
        entry["Delta_f_deloc"] = round(
            entry["isotopes"]["H"]["f_deloc"] - entry["isotopes"]["D"]["f_deloc"], 4
        )
        results["by_barrier"].append(entry)

    path = os.path.join("results", "adw_landscapes.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print("Saved", path)
    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--n_steps", type=int, default=80_000)
    p.add_argument("--n_equil", type=int, default=12_000)
    a = p.parse_args()
    run(a.n_steps, a.n_equil)
