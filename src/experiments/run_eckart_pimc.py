#!/usr/bin/env python3
"""
Centroid QTST on the symmetric Eckart barrier via umbrella PIMC.

k_cTST = sqrt(kT/2πm) exp(-β F(0)), with F(∞)=0 from WHAM.
Compare to analytic k_exact, k_cl, k_Wigner.
Also histogram χ_p at the dividing-surface window.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from models.qrr_core import (  # noqa: E402
    AMU_SI,
    KB_KJMOL,
    KB_SI,
    eckart_classical_tst,
    eckart_thermal_rate,
    make_rng,
    pimc_1d,
    wham_1d,
    wigner_kappa,
    eckart_barrier_frequency,
)

T_K = 300.0
A_ANG = 0.50
CHI_C = 0.010
V0_KT = [4.0, 8.0, 12.0]
P_H, P_D = 32, 24
SEED = 42


def eckart_V_kJ(q_A, V0_kJ, a=A_ANG):
    x = np.asarray(q_A, dtype=float) / a
    # sech^2 x = 1/cosh^2 x
    return V0_kJ / np.cosh(x) ** 2


def run(n_steps: int = 40_000, n_equil: int = 8_000) -> dict:
    os.makedirs("results", exist_ok=True)
    rng = make_rng(SEED)
    kT = KB_KJMOL * T_K
    windows = np.linspace(-2.0, 2.0, 11)
    k_umb = 60.0
    xi_bins = np.linspace(-2.4, 2.4, 49)
    xi_cent = 0.5 * (xi_bins[:-1] + xi_bins[1:])
    dxi = xi_cent[1] - xi_cent[0]

    results = {
        "description": "Eckart centroid QTST from umbrella PIMC vs exact analytic rates",
        "T_K": T_K,
        "a_angstrom": A_ANG,
        "V0_over_kT": V0_KT,
        "P_H": P_H,
        "P_D": P_D,
        "n_steps": n_steps,
        "windows": windows.tolist(),
        "k_umb": k_umb,
        "chi_c_ang2": CHI_C,
        "seed": SEED,
        "by_barrier": [],
    }
    print("Eckart PIMC centroid QTST")

    for vkt in V0_KT:
        V0_kJ = vkt * kT
        V0_J = vkt * KB_SI * T_K
        entry = {"V0_over_kT": vkt, "isotopes": {}}
        print(f"  V0 = {vkt} kT")
        for iso, mass, P in (("H", 1.0, P_H), ("D", 2.0, P_D)):

            def V(q, _V0=V0_kJ):
                return eckart_V_kJ(q, _V0)

            counts = np.zeros((len(windows), len(xi_cent)))
            N_win = np.zeros(len(windows))
            U_bias = np.zeros((len(windows), len(xi_cent)))
            chi_ds = []
            xi_win0_mean = 0.0
            for i, x0 in enumerate(windows):
                U_bias[i] = 0.5 * k_umb * (xi_cent - x0) ** 2
                samp = pimc_1d(
                    V, T_K, mass, P, n_steps, n_equil, x0, rng,
                    k_umb=k_umb, q_umb=x0, centroid_step=0.04,
                )
                xi_s, chi_s = samp["xi"], samp["chi"]
                hist, _ = np.histogram(xi_s, bins=xi_bins)
                counts[i] = hist
                N_win[i] = len(xi_s)
                if abs(x0) < 1e-12:
                    chi_ds = chi_s[np.abs(xi_s) < 0.12].copy()
                    xi_win0_mean = float(xi_s.mean())
                print(f"    {iso} x0={x0:+.2f}  <xi>={xi_s.mean():+.3f}  <chi>={chi_s.mean():.4f}")

            F, Pxi, f_i = wham_1d(counts, N_win, U_bias, T_K, dxi)
            # Zero at the edges (reactant/product asymptote)
            edge = (np.abs(xi_cent) > 1.6)
            F_inf = float(np.mean(F[edge])) if np.any(edge) else float(F[0])
            F0 = float(F[int(np.argmin(np.abs(xi_cent)))])
            dF = F0 - F_inf  # kJ/mol quantum barrier

            k_cl = eckart_classical_tst(T_K, V0_J, mass * AMU_SI)
            k_ex = eckart_thermal_rate(T_K, V0_J, A_ANG * 1e-10, mass * AMU_SI)
            w_b = eckart_barrier_frequency(V0_J, A_ANG * 1e-10, mass * AMU_SI)
            k_w = k_cl * wigner_kappa(T_K, w_b)
            # centroid QTST: same prefactor, quantum barrier
            k_ctst = k_cl * np.exp(-(dF - V0_kJ) / kT)

            if len(chi_ds) < 20:
                chi_ds = np.array([0.01])
            f_del = float(np.mean(chi_ds >= 2.0 * CHI_C))
            chi_mean = float(np.mean(chi_ds))

            # QRR with κ=1 equals k_cTST identically; store both.
            # Flag the point if the DS umbrella drifted off ξ=0.
            reliable = abs(xi_win0_mean) < 0.10
            entry["isotopes"][iso] = {
                "DeltaF_kJmol": round(dF, 4),
                "V0_kJmol": round(V0_kJ, 4),
                "k_cl_mps": k_cl,
                "k_wigner_mps": k_w,
                "k_ctst_mps": float(k_ctst),
                "k_qrr_mps": float(k_ctst),  # κ = 1
                "k_exact_mps": k_ex,
                "kappa_cl": 1.0,
                "kappa_wigner": float(k_w / k_cl),
                "kappa_ctst": float(k_ctst / k_cl),
                "kappa_exact": float(k_ex / k_cl),
                "f_deloc": round(f_del, 4),
                "chi_mean_ds": round(chi_mean, 6),
                "n_ds": int(len(chi_ds)),
                "xi_win0_mean": round(xi_win0_mean, 4),
                "xi_cent": np.round(xi_cent, 4).tolist(),
                "F_xi": np.round(F - F_inf, 4).tolist(),
                "chi_ds_mean": round(chi_mean, 6),
                "reliable": reliable,
            }
            print(f"    {iso}  ΔF={dF:.3f} (V0={V0_kJ:.3f})  "
                  f"κ_cTST={k_ctst/k_cl:.3f}  κ_ex={k_ex/k_cl:.3f}  "
                  f"κ_W={k_w/k_cl:.3f}  f_deloc={f_del:.3f}")

        h = entry["isotopes"]["H"]
        dd = entry["isotopes"]["D"]
        # Unphysical isotope order of the centroid barrier.
        if dd["DeltaF_kJmol"] < h["DeltaF_kJmol"]:
            dd["reliable"] = False
        entry["KIE_exact"] = h["k_exact_mps"] / dd["k_exact_mps"]
        entry["KIE_ctst"] = h["k_ctst_mps"] / dd["k_ctst_mps"]
        entry["KIE_cl"] = h["k_cl_mps"] / dd["k_cl_mps"]
        entry["KIE_wigner"] = h["k_wigner_mps"] / dd["k_wigner_mps"]
        results["by_barrier"].append(entry)

    path = os.path.join("results", "eckart_pimc.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print("Saved", path)
    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--n_steps", type=int, default=40_000)
    p.add_argument("--n_equil", type=int, default=8_000)
    a = p.parse_args()
    run(a.n_steps, a.n_equil)
