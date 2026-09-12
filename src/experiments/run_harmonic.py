#!/usr/bin/env python3
"""Harmonic-oscillator benchmark for chi_p: continuum, exact finite-P, PIMC."""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from models.qrr_core import (  # noqa: E402
    block_mean_err,
    chi_harmonic_continuum,
    chi_harmonic_finite_P,
    make_rng,
    omega_cm1_to_rad_s,
    pimc_1d,
)

OMEGA_CM1 = 1000.0
OMEGA_H = omega_cm1_to_rad_s(OMEGA_CM1)
OMEGA_D = OMEGA_H / np.sqrt(2.0)
T_RANGE = np.array([100, 150, 200, 250, 300, 350, 400, 500, 600, 700], dtype=float)
P_VALUES = [4, 8, 16, 32, 64, 128]
SEED = 42


def V_harm_H(q):
    # 0.5 m ω^2 q^2 in kJ/mol, q in Å
    from models.qrr_core import AMU_SI, KJMOL_SI
    k = 1.0 * AMU_SI * OMEGA_H**2              # N/m = J/m^2
    k_kJ_A2 = k / KJMOL_SI / 1e20              # kJ/mol / Å^2
    return 0.5 * k_kJ_A2 * np.asarray(q) ** 2


def V_harm_D(q):
    from models.qrr_core import AMU_SI, KJMOL_SI
    k = 2.0 * AMU_SI * OMEGA_D**2
    k_kJ_A2 = k / KJMOL_SI / (1e20)
    return 0.5 * k_kJ_A2 * np.asarray(q) ** 2


def run(n_steps: int = 200_000, n_equil: int = 40_000) -> dict:
    os.makedirs("results", exist_ok=True)
    rng = make_rng(SEED)

    results = {
        "description": "Harmonic oscillator chi_p: continuum, exact finite-P, Levy-staging PIMC",
        "omega_p_cm1": OMEGA_CM1,
        "m_H_amu": 1.0,
        "m_D_amu": 2.0,
        "seed": SEED,
        "n_steps": n_steps,
        "n_equil": n_equil,
        "temperatures_K": T_RANGE.tolist(),
        "P_H": 64,
        "P_D": 32,
        "continuum": {"H": [], "D": []},
        "finite_P": {"H": [], "D": []},
        "pimc": {"H_mean": [], "H_std": [], "D_mean": [], "D_std": []},
        "ratio_continuum": [],
        "ratio_finite_P": [],
        "ratio_pimc": [],
        "bead_convergence": {
            "T_K": 300.0,
            "P_values": P_VALUES,
            "continuum_H": None,
            "continuum_D": None,
            "finite_P_H": [],
            "finite_P_D": [],
            "pimc_H": [],
            "pimc_H_std": [],
            "pimc_D": [],
            "pimc_D_std": [],
        },
    }

    print("Harmonic benchmark")
    print(f"  omega = {OMEGA_CM1} cm^-1   n_steps={n_steps}")

    for T in T_RANGE:
        cH = chi_harmonic_continuum(T, 1.0, OMEGA_H)
        cD = chi_harmonic_continuum(T, 2.0, OMEGA_D)
        fH = chi_harmonic_finite_P(T, 1.0, OMEGA_H, 64)
        fD = chi_harmonic_finite_P(T, 2.0, OMEGA_D, 32)
        results["continuum"]["H"].append(round(cH, 7))
        results["continuum"]["D"].append(round(cD, 7))
        results["finite_P"]["H"].append(round(fH, 7))
        results["finite_P"]["D"].append(round(fD, 7))
        results["ratio_continuum"].append(round(cH / cD, 5))
        results["ratio_finite_P"].append(round(fH / fD, 5))

        samp_H = pimc_1d(V_harm_H, T, 1.0, 64, n_steps, n_equil, 0.0, rng)
        samp_D = pimc_1d(V_harm_D, T, 2.0, 32, n_steps, n_equil, 0.0, rng)
        mH, sH = block_mean_err(samp_H["chi"])
        mD, sD = block_mean_err(samp_D["chi"])
        results["pimc"]["H_mean"].append(round(mH, 6))
        results["pimc"]["H_std"].append(round(sH, 6))
        results["pimc"]["D_mean"].append(round(mD, 6))
        results["pimc"]["D_std"].append(round(sD, 6))
        results["ratio_pimc"].append(round(mH / mD, 5))
        print(f"  T={T:4.0f}  cont_H={cH:.5f}  finP_H={fH:.5f}  PIMC_H={mH:.5f}±{sH:.5f}")

    print("Bead convergence at 300 K")
    results["bead_convergence"]["continuum_H"] = round(chi_harmonic_continuum(300, 1.0, OMEGA_H), 7)
    results["bead_convergence"]["continuum_D"] = round(chi_harmonic_continuum(300, 2.0, OMEGA_D), 7)
    for P in P_VALUES:
        fH = chi_harmonic_finite_P(300, 1.0, OMEGA_H, P)
        fD = chi_harmonic_finite_P(300, 2.0, OMEGA_D, P)
        samp_H = pimc_1d(V_harm_H, 300.0, 1.0, P, n_steps, n_equil, 0.0, rng)
        samp_D = pimc_1d(V_harm_D, 300.0, 2.0, P, n_steps, n_equil, 0.0, rng)
        mH, sH = block_mean_err(samp_H["chi"])
        mD, sD = block_mean_err(samp_D["chi"])
        results["bead_convergence"]["finite_P_H"].append(round(fH, 7))
        results["bead_convergence"]["finite_P_D"].append(round(fD, 7))
        results["bead_convergence"]["pimc_H"].append(round(mH, 6))
        results["bead_convergence"]["pimc_H_std"].append(round(sH, 6))
        results["bead_convergence"]["pimc_D"].append(round(mD, 6))
        results["bead_convergence"]["pimc_D_std"].append(round(sD, 6))
        print(f"  P={P:3d}  finP_H={fH:.5f}  PIMC_H={mH:.5f}±{sH:.5f}")

    # RMSE of PIMC vs exact finite-P (the quantity PIMC actually samples)
    fin = np.array(results["finite_P"]["H"])
    pmc = np.array(results["pimc"]["H_mean"])
    results["rmse_pimc_vs_finiteP_H_percent"] = round(float(np.sqrt(np.mean(((pmc - fin) / fin) ** 2)) * 100), 3)
    finD = np.array(results["finite_P"]["D"])
    pmcD = np.array(results["pimc"]["D_mean"])
    results["rmse_pimc_vs_finiteP_D_percent"] = round(float(np.sqrt(np.mean(((pmcD - finD) / finD) ** 2)) * 100), 3)
    print(f"RMSE PIMC vs finite-P: H {results['rmse_pimc_vs_finiteP_H_percent']}%  D {results['rmse_pimc_vs_finiteP_D_percent']}%")

    path = os.path.join("results", "harmonic_benchmark.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print("Saved", path)
    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--n_steps", type=int, default=200_000)
    p.add_argument("--n_equil", type=int, default=40_000)
    a = p.parse_args()
    run(a.n_steps, a.n_equil)
