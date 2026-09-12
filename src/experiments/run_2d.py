#!/usr/bin/env python3
"""
2D proton–solvent model: factorization diagnostic Cov(chi, q_s | dividing surface).
V = V_ADW(q_p) + (1/2) k_s q_s^2 + kappa_c * q_p * q_s
Proton is a ring polymer; solvent is classical.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from models.qrr_core import (  # noqa: E402
    HBAR2_OVER_AMU,
    KB_KJMOL,
    adw_V,
    make_rng,
)

T_K = 300.0
V0_KT = 4.0
K_S = 2.0          # kJ/mol/Å^2  (soft solvent)
KAPPA = [0.0, 0.5, 1.0, 2.0, 4.0]  # kJ/mol/Å^2
P = 32
MASS = 1.0
SEED = 42
CHI_C = 0.010


def _levy_bridge(q0, qend, L, lam, rng):
    interior = np.empty(L - 1)
    prev = q0
    for i in range(1, L):
        n_rem = L - i
        mu = (n_rem * prev + qend) / (n_rem + 1)
        sig = np.sqrt(lam * n_rem / (n_rem + 1))
        qn = rng.normal(mu, sig)
        interior[i - 1] = qn
        prev = qn
    return interior


def run(n_steps: int = 120_000, n_equil: int = 20_000) -> dict:
    os.makedirs("results", exist_ok=True)
    rng = make_rng(SEED)
    kT = KB_KJMOL * T_K
    beta = 1.0 / kT
    V0 = V0_KT * kT
    lam = (HBAR2_OVER_AMU / MASS) * beta / P
    L = max(2, P // 4)

    results = {
        "description": "2D proton-solvent factorization diagnostic",
        "T_K": T_K,
        "V0_over_kT": V0_KT,
        "k_s": K_S,
        "kappa_c": KAPPA,
        "P": P,
        "n_steps": n_steps,
        "chi_c_ang2": CHI_C,
        "scan": [],
    }
    print("2D coupling scan")

    for kc in KAPPA:
        beads = np.full(P, -0.70) + rng.normal(0.0, 0.03, P)
        qs = 0.0
        chi_ds = []
        qs_ds = []
        chi_all = []

        def energy(b, s):
            return float(np.sum(adw_V(b, V0))) / P + 0.5 * K_S * s**2 + kc * b.mean() * s

        E = energy(beads, qs)
        for step in range(n_equil + n_steps):
            i0 = int(rng.integers(0, P))
            i_end = (i0 + L) % P
            trial = beads.copy()
            interior = _levy_bridge(beads[i0], beads[i_end], L, lam, rng)
            for j, qn in enumerate(interior):
                trial[(i0 + 1 + j) % P] = qn
            Et = energy(trial, qs)
            if Et <= E or rng.random() < np.exp(-beta * (Et - E)):
                beads = trial
                E = Et
            # solvent Metropolis
            qs_t = qs + rng.normal(0.0, 0.08)
            Et = energy(beads, qs_t)
            if Et <= E or rng.random() < np.exp(-beta * (Et - E)):
                qs = qs_t
                E = Et
            # centroid kick
            b2 = beads + rng.normal(0.0, 0.03)
            Et = energy(b2, qs)
            if Et <= E or rng.random() < np.exp(-beta * (Et - E)):
                beads = b2
                E = Et

            if step >= n_equil:
                c = beads.mean()
                ch = float(np.mean((beads - c) ** 2))
                chi_all.append(ch)
                if abs(c) < 0.10:
                    chi_ds.append(ch)
                    qs_ds.append(qs)

        chi_ds = np.asarray(chi_ds)
        qs_ds = np.asarray(qs_ds)
        chi_all = np.asarray(chi_all)
        if len(chi_ds) < 30:
            cov = 0.0
            corr = 0.0
            f_del = float(np.mean(chi_all >= 2 * CHI_C))
        else:
            cov = float(np.mean((chi_ds - chi_ds.mean()) * (qs_ds - qs_ds.mean())))
            sx = chi_ds.std()
            sy = qs_ds.std()
            corr = float(cov / (sx * sy)) if sx > 0 and sy > 0 else 0.0
            f_del = float(np.mean(chi_ds >= 2 * CHI_C))
        entry = {
            "kappa_c": kc,
            "n_ds": int(len(chi_ds)),
            "cov_chi_qs": round(cov, 8),
            "corr_chi_qs": round(corr, 5),
            "f_deloc_ds": round(f_del, 4),
            "chi_mean_ds": round(float(chi_ds.mean()) if len(chi_ds) else float(chi_all.mean()), 6),
        }
        results["scan"].append(entry)
        print(f"  kappa={kc:4.1f}  n_ds={len(chi_ds):5d}  corr={corr:+.3f}  f_deloc={f_del:.3f}")

    path = os.path.join("results", "coupling_2d.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print("Saved", path)
    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--n_steps", type=int, default=120_000)
    p.add_argument("--n_equil", type=int, default=20_000)
    a = p.parse_args()
    run(a.n_steps, a.n_equil)
