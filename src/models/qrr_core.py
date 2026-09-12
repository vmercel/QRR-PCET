"""
qrr_core.py
===========
Shared constants, exact formulae, and path-integral samplers for QRR-PCET.

Internal PIMC units
-------------------
length : angstrom
energy : kJ/mol
mass   : amu
temperature : kelvin
"""

from __future__ import annotations

import numpy as np
from numpy.random import Generator, default_rng
from scipy.integrate import quad

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
HBAR_SI = 1.054_571_817e-34          # J s
KB_SI = 1.380_649e-23                # J/K
AMU_SI = 1.660_539_066_60e-27        # kg
ANG_SI = 1.0e-10                     # m
NA = 6.022_140_76e23
KJMOL_SI = 1.0e3 / NA                # J per molecule per (kJ/mol)

# ħ in kJ mol^-1 ps ; used only for documentation
HBAR_KJPS = HBAR_SI / (KJMOL_SI * 1e-12)

# k_B in kJ mol^-1 K^-1
KB_KJMOL = KB_SI / KJMOL_SI          # 0.008314462618

# ħ^2 / amu in (kJ/mol) * Å^2
# ħ^2 / m  = (HBAR_SI^2 / AMU_SI) / KJMOL_SI / ANG_SI^2
HBAR2_OVER_AMU = (HBAR_SI**2 / AMU_SI) / KJMOL_SI / (ANG_SI**2)

CM1_TO_RAD_S = 2.0 * np.pi * 2.997_924_58e10


def omega_cm1_to_rad_s(nu_cm1: float) -> float:
    return nu_cm1 * CM1_TO_RAD_S


# ---------------------------------------------------------------------------
# Harmonic oscillator: continuum and exact finite-P
# ---------------------------------------------------------------------------
def chi_harmonic_continuum(T_K: float, mass_amu: float, omega_rad_s: float) -> float:
    """<chi_p> in Å^2 for a harmonic oscillator, P -> infinity."""
    m = mass_amu * AMU_SI
    beta = 1.0 / (KB_SI * T_K)
    x = 0.5 * beta * HBAR_SI * omega_rad_s
    if x > 40.0:
        coth = 1.0
    elif x < 1e-12:
        coth = 1.0 / x
    else:
        coth = 1.0 / np.tanh(x)
    zpm = HBAR_SI / (2.0 * m * omega_rad_s) * coth
    centroid = KB_SI * T_K / (m * omega_rad_s**2)
    return max(zpm - centroid, 0.0) / ANG_SI**2


def chi_harmonic_finite_P(T_K: float, mass_amu: float, omega_rad_s: float, P: int) -> float:
    """Exact primitive <chi_p> in Å^2 for a P-bead harmonic ring polymer."""
    m = mass_amu * AMU_SI
    beta = 1.0 / (KB_SI * T_K)
    k = np.arange(1, P)
    sin_k = np.sin(np.pi * k / P)
    lam_k = (4.0 * m * P) / (HBAR_SI**2 * beta**2) * sin_k**2 + m * omega_rad_s**2 / P
    chi_m2 = np.sum(1.0 / (beta * lam_k)) / P
    return chi_m2 / ANG_SI**2


# ---------------------------------------------------------------------------
# Symmetric Eckart barrier (analytic transmission)
# V(x) = V0 sech^2(x/a)
# ---------------------------------------------------------------------------
def eckart_T(E_J: np.ndarray, V0_J: float, a_m: float, mass_kg: float) -> np.ndarray:
    """
    Symmetric Eckart transmission (Connor).
    V(x) = V0 sech^2(x/a),  alpha = a sqrt(2 m E)/ħ,
    delta = (1/2) sqrt(8 m V0 a^2/ħ^2 - 1).
    T = sinh^2(pi alpha) / [sinh^2(pi alpha) + cosh^2(pi delta)].
    Evaluated in log space to avoid overflow.
    """
    E = np.asarray(E_J, dtype=float)
    T = np.zeros_like(E, dtype=float)
    pos = E > 0.0
    if not np.any(pos):
        return T
    alpha = (a_m / HBAR_SI) * np.sqrt(2.0 * mass_kg * E[pos])
    disc = 8.0 * mass_kg * V0_J * a_m**2 / HBAR_SI**2 - 1.0

    def _log_cosh(x):
        ax = np.abs(x)
        return ax + np.log1p(np.exp(-2.0 * ax)) - np.log(2.0)

    def _log_sinh(x):
        x = np.maximum(x, 1e-16)
        return x + np.log1p(-np.exp(-2.0 * x)) - np.log(2.0)

    if disc >= 0.0:
        delta = 0.5 * np.sqrt(disc)
        log_ratio = 2.0 * _log_cosh(np.pi * delta) - 2.0 * _log_sinh(np.pi * alpha)
        T[pos] = 1.0 / (1.0 + np.exp(log_ratio))
    else:
        delta = 0.5 * np.sqrt(-disc)
        # cosh^2 -> cos^2, which can vanish; fall back to direct eval
        sinh2 = np.sinh(np.pi * alpha) ** 2
        den = sinh2 + np.cos(np.pi * delta) ** 2
        T[pos] = sinh2 / np.maximum(den, 1e-300)
    return np.clip(T, 0.0, 1.0)


def eckart_thermal_rate(T_K: float, V0_J: float, a_m: float, mass_kg: float) -> float:
    """
    One-dimensional thermal rate (velocity units, m/s):
        k = [1 / (2 pi ħ Q_r)] ∫ dE e^{-βE} T(E)
    with Q_r = sqrt(m / (2 pi ħ^2 β)) the reactant translational
    partition function per unit length.  Equivalently
        k = (1/h) (2 pi ħ^2 / m)^{1/2} β^{-1/2}  *  (β ∫ dE e^{-βE} T(E)).
    Returned in m/s.
    """
    beta = 1.0 / (KB_SI * T_K)
    # integrand of ∫ dE e^{-βE} T(E)
    def integrand(E):
        if E <= 0.0:
            return 0.0
        t = float(eckart_T(np.array([E]), V0_J, a_m, mass_kg)[0])
        return np.exp(-beta * E) * t

    # integrate to a safe cutoff
    E_max = V0_J + 40.0 / beta
    integ, _ = quad(integrand, 0.0, E_max, epsabs=1e-20, limit=400)
    k = integ / (2.0 * np.pi * HBAR_SI * np.sqrt(mass_kg / (2.0 * np.pi * HBAR_SI**2 * beta)))
    return k


def eckart_classical_tst(T_K: float, V0_J: float, mass_kg: float) -> float:
    """k_cl = sqrt(kT / 2 pi m) exp(-β V0), m/s."""
    return np.sqrt(KB_SI * T_K / (2.0 * np.pi * mass_kg)) * np.exp(-V0_J / (KB_SI * T_K))


def eckart_barrier_frequency(V0_J: float, a_m: float, mass_kg: float) -> float:
    """ω_b = sqrt(-V''(0)/m) = sqrt(2 V0 / (m a^2)), rad/s."""
    return np.sqrt(2.0 * V0_J / (mass_kg * a_m**2))


def wigner_kappa(T_K: float, omega_b: float) -> float:
    """Leading Wigner tunneling factor 1 + (β ħ ω_b)^2 / 24."""
    x = HBAR_SI * omega_b / (KB_SI * T_K)
    return 1.0 + x**2 / 24.0


# ---------------------------------------------------------------------------
# Asymmetric double well (1D)
# V(q) = V0 [(q/q0)^2 - 1]^2 + alpha q
# ---------------------------------------------------------------------------
def adw_V(q_A: np.ndarray, V0_kJ: float, q0_A: float = 0.70, alpha_kJ_A: float = 0.025) -> np.ndarray:
    q = np.asarray(q_A, dtype=float)
    return V0_kJ * ((q / q0_A) ** 2 - 1.0) ** 2 + alpha_kJ_A * q


# ---------------------------------------------------------------------------
# Levy staging PIMC
# ---------------------------------------------------------------------------
def _levy_bridge(q0: float, qend: float, L: int, lam: float, rng: Generator) -> np.ndarray:
    """
    Interior beads of an L-slice free-particle bridge (L-1 values).
    Sequential Levy construction; lam = ħ^2 β / (m P) in Å^2.
    """
    interior = np.empty(L - 1)
    prev = q0
    for i in range(1, L):
        n_rem = L - i
        mu = (n_rem * prev + qend) / (n_rem + 1)
        sig = np.sqrt(lam * n_rem / (n_rem + 1))
        q_new = rng.normal(mu, sig)
        interior[i - 1] = q_new
        prev = q_new
    return interior


def pimc_1d(
    V_func,
    T_K: float,
    mass_amu: float,
    P: int,
    n_steps: int,
    n_equil: int,
    q_init: float,
    rng: Generator,
    stage_length: int | None = None,
    k_umb: float = 0.0,
    q_umb: float = 0.0,
    centroid_step: float = 0.04,
) -> dict:
    """
    1D ring-polymer Monte Carlo with Levy staging + centroid moves.

    V_func(q_A) returns energy in kJ/mol.
    Umbrella (optional): (k_umb/2)(q_centroid - q_umb)^2, k_umb in kJ/mol/Å^2.

    Returns dict with arrays of centroid and chi (Å, Å^2) of length n_steps.
    """
    beta = 1.0 / (KB_KJMOL * T_K)
    lam = (HBAR2_OVER_AMU / mass_amu) * beta / P          # Å^2
    if stage_length is None:
        stage_length = max(2, P // 4)
    L = min(stage_length, P - 1)

    beads = np.full(P, q_init, dtype=float) + rng.normal(0.0, 0.02, P)
    xi = np.empty(n_steps)
    chi = np.empty(n_steps)

    def on_site(q):
        e = float(np.sum(V_func(q))) / P
        if k_umb != 0.0:
            e += 0.5 * k_umb * (q.mean() - q_umb) ** 2
        return e

    V_curr = on_site(beads)

    for step in range(n_equil + n_steps):
        # Levy staging of L consecutive beads (endpoints of the segment fixed)
        i0 = int(rng.integers(0, P))
        i_end = (i0 + L) % P
        trial = beads.copy()
        interior = _levy_bridge(beads[i0], beads[i_end], L, lam, rng)
        for j, qn in enumerate(interior):
            trial[(i0 + 1 + j) % P] = qn
        V_try = on_site(trial)
        dV = V_try - V_curr
        if dV <= 0.0 or rng.random() < np.exp(-beta * dV):
            beads = trial
            V_curr = V_try

        # rigid centroid move (helps the k = 0 mode)
        dq = rng.normal(0.0, centroid_step)
        trial_c = beads + dq
        V_try = on_site(trial_c)
        dV = V_try - V_curr
        if dV <= 0.0 or rng.random() < np.exp(-beta * dV):
            beads = trial_c
            V_curr = V_try

        if step >= n_equil:
            c = beads.mean()
            xi[step - n_equil] = c
            chi[step - n_equil] = np.mean((beads - c) ** 2)

    return {"xi": xi, "chi": chi}


def block_mean_err(x: np.ndarray, n_blocks: int = 10) -> tuple[float, float]:
    n = len(x)
    bs = n // n_blocks
    if bs < 1:
        return float(np.mean(x)), float(np.std(x))
    means = np.array([x[b * bs:(b + 1) * bs].mean() for b in range(n_blocks)])
    return float(means.mean()), float(means.std(ddof=1) / np.sqrt(n_blocks))


# ---------------------------------------------------------------------------
# 1D WHAM
# ---------------------------------------------------------------------------
def wham_1d(
    counts: np.ndarray,
    N_win: np.ndarray,
    U_bias: np.ndarray,
    T_K: float,
    dxi: float,
    n_iter: int = 400,
    tol: float = 1e-10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    counts : (n_win, n_bin)
    U_bias : (n_win, n_bin) bias energy at bin centres, kJ/mol
    Returns F (kJ/mol, min 0), P(xi), f_i (kJ/mol).
    """
    beta = 1.0 / (KB_KJMOL * T_K)
    n_win, n_bin = counts.shape
    f = np.zeros(n_win)
    num = counts.sum(axis=0)
    log_N = np.log(np.maximum(N_win, 1.0))
    for _ in range(n_iter):
        # log den_j = log sum_i N_i exp(-β(U_ij - f_i))
        log_terms = log_N[:, None] - beta * (U_bias - f[:, None])
        m = np.max(log_terms, axis=0)
        den = np.exp(m) * np.exp(log_terms - m).sum(axis=0)
        P = num / np.maximum(den, 1e-300)
        P = np.maximum(P, 0.0)
        s = P.sum() * dxi
        if s > 0:
            P /= s
        f_new = np.empty(n_win)
        for i in range(n_win):
            z = np.sum(P * dxi * np.exp(-beta * U_bias[i]))
            f_new[i] = -KB_KJMOL * T_K * np.log(max(z, 1e-300))
        f_new -= f_new[0]
        if np.max(np.abs(f_new - f)) < tol:
            f = f_new
            break
        f = f_new
    F = -KB_KJMOL * T_K * np.log(np.maximum(P, 1e-40))
    F -= np.min(F)
    return F, P, f


def make_rng(seed: int = 42) -> Generator:
    return default_rng(seed)
