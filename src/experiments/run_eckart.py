#!/usr/bin/env python3
"""Analytic Eckart-barrier thermal rates, KIEs, and Wigner corrections."""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from models.qrr_core import (  # noqa: E402
    AMU_SI,
    KB_SI,
    eckart_barrier_frequency,
    eckart_classical_tst,
    eckart_thermal_rate,
    wigner_kappa,
)

T_K = 300.0
A_ANG = 0.50
A_M = A_ANG * 1e-10
V0_KT = [4.0, 6.0, 8.0, 10.0, 12.0, 16.0]


def run() -> dict:
    os.makedirs("results", exist_ok=True)
    kT = KB_SI * T_K
    results = {
        "description": "Symmetric Eckart barrier: exact thermal rates and KIEs",
        "T_K": T_K,
        "a_angstrom": A_ANG,
        "V0_over_kT": V0_KT,
        "m_H_amu": 1.0,
        "m_D_amu": 2.0,
        "k_exact_H_mps": [],
        "k_exact_D_mps": [],
        "k_tst_H_mps": [],
        "k_tst_D_mps": [],
        "k_wigner_H_mps": [],
        "k_wigner_D_mps": [],
        "kappa_exact_H": [],
        "kappa_exact_D": [],
        "kappa_wigner_H": [],
        "kappa_wigner_D": [],
        "KIE_exact": [],
        "KIE_tst": [],
        "KIE_wigner": [],
        "hbar_omega_b_over_kT_H": [],
        "hbar_omega_b_over_kT_D": [],
    }
    print(f"Eckart benchmark  T={T_K:.0f} K  a={A_ANG:.2f} Å")
    for vkt in V0_KT:
        V0 = vkt * kT
        kH = eckart_thermal_rate(T_K, V0, A_M, 1.0 * AMU_SI)
        kD = eckart_thermal_rate(T_K, V0, A_M, 2.0 * AMU_SI)
        tH = eckart_classical_tst(T_K, V0, 1.0 * AMU_SI)
        tD = eckart_classical_tst(T_K, V0, 2.0 * AMU_SI)
        wH = eckart_barrier_frequency(V0, A_M, 1.0 * AMU_SI)
        wD = eckart_barrier_frequency(V0, A_M, 2.0 * AMU_SI)
        kapW_H = wigner_kappa(T_K, wH)
        kapW_D = wigner_kappa(T_K, wD)
        from models.qrr_core import HBAR_SI
        results["k_exact_H_mps"].append(kH)
        results["k_exact_D_mps"].append(kD)
        results["k_tst_H_mps"].append(tH)
        results["k_tst_D_mps"].append(tD)
        results["k_wigner_H_mps"].append(tH * kapW_H)
        results["k_wigner_D_mps"].append(tD * kapW_D)
        results["kappa_exact_H"].append(kH / tH)
        results["kappa_exact_D"].append(kD / tD)
        results["kappa_wigner_H"].append(kapW_H)
        results["kappa_wigner_D"].append(kapW_D)
        results["KIE_exact"].append(kH / kD)
        results["KIE_tst"].append(tH / tD)
        results["KIE_wigner"].append((tH * kapW_H) / (tD * kapW_D))
        results["hbar_omega_b_over_kT_H"].append(HBAR_SI * wH / kT)
        results["hbar_omega_b_over_kT_D"].append(HBAR_SI * wD / kT)
        print(f"  V0={vkt:4.1f} kT  kappa_H={kH/tH:.3f}  KIE_exact={kH/kD:.3f}  KIE_TST={tH/tD:.3f}")

    # round for JSON readability while keeping full precision in lists via python floats
    def rnd(xs, n=8):
        return [round(float(x), n) for x in xs]

    for key, val in list(results.items()):
        if isinstance(val, list) and val and isinstance(val[0], float):
            results[key] = rnd(val, 8)

    path = os.path.join("results", "eckart_rates.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print("Saved", path)
    return results


if __name__ == "__main__":
    run()
