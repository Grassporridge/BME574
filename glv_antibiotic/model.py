"""Random diagonally-stable gLV community with direct antibiotic killing.

No beta-lactamase, inhibitor, or explicit resource state is included. Ecological
interactions are represented by the gLV matrix; antibiotic is a decaying scalar.
"""
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


@dataclass
class Config:
    seed: int = 42
    n_species: int = 25
    t_end: float = 300.0
    n_times: int = 501
    initial_biomass: float = 0.01
    equilibrium_abundance: float = 0.20
    interaction_strength: float = 0.25
    antibiotic_initial: float = 1.0
    antibiotic_decay: float = 0.02
    antibiotic_kmax: float = 0.50
    antibiotic_hill: float = 2.0
    ic50_min: float = 0.05
    ic50_max: float = 2.0
    resistance_cost_max: float = 0.20
    extinction_fraction: float = 1e-6


def sample_community(c: Config):
    """Sample A=W^-1(K-P), with positive-definite P and a stability certificate.

    Rows are scaled so A_ii=-1. The certificate is W=diag(1/scale), giving
    W A + A.T W = -2P (negative definite). Retry only to ensure r=-A x*>0.
    """
    if c.n_species < 2 or c.interaction_strength < 0:
        raise ValueError("Need at least two species and nonnegative interaction strength")
    rng = np.random.default_rng(c.seed)
    n = c.n_species
    for _ in range(1000):
        # P is symmetric positive definite; its diagonal is normalized away below.
        L = rng.normal(0.0, c.interaction_strength / np.sqrt(n), size=(n, n))
        P = np.eye(n) + L @ L.T
        G = rng.normal(0.0, c.interaction_strength / np.sqrt(n), size=(n, n))
        K = 0.5 * (G - G.T)
        scale = 1.0 / np.diag(P)
        A = scale[:, None] * (K - P)
        W = np.diag(1.0 / scale)
        x_star = rng.lognormal(np.log(c.equilibrium_abundance), 0.12, n)
        r = -A @ x_star
        if np.all(r > 0):
            break
    else:
        raise RuntimeError("Could not sample a feasible positive-growth community")

    resistance_rng = np.random.default_rng(c.seed + 200000)
    resistance = resistance_rng.uniform(0.0, 1.0, n)
    ic50 = c.ic50_min * (c.ic50_max / c.ic50_min) ** resistance
    species = pd.DataFrame({
        "species": [f"S{i+1:02d}" for i in range(n)],
        "intrinsic_growth_r": r,
        "resistance_rho": resistance,
        "IC50": ic50,
        "fitness_cost_at_max": c.resistance_cost_max * resistance,
        "baseline_equilibrium_abundance": x_star,
        "initial_abundance": c.initial_biomass / n,
    })
    return {"A": A, "P": P, "K": K, "W": W, "r": r,
            "x_star": x_star, "resistance": resistance, "IC50": ic50,
            "ic50_min": c.ic50_min, "ic50_max": c.ic50_max,
            "resistance_cost_max": c.resistance_cost_max,
            "species": species}


def antibiotic_kill(A, ic50, kmax, hill):
    A = max(float(A), 0.0)
    return kmax * A**hill / (ic50**hill + A**hill)


def simulate(c: Config, community=None, *, initial_biomass=None,
             antibiotic_initial=None, antibiotic_decay=None,
             antibiotic_kmax=None, antibiotic_hill=None,
             resistance_cost_max=None):
    p = community or sample_community(c)
    n = c.n_species
    biomass0 = c.initial_biomass if initial_biomass is None else initial_biomass
    A0 = c.antibiotic_initial if antibiotic_initial is None else antibiotic_initial
    kdeg = c.antibiotic_decay if antibiotic_decay is None else antibiotic_decay
    kmax = c.antibiotic_kmax if antibiotic_kmax is None else antibiotic_kmax
    hill = c.antibiotic_hill if antibiotic_hill is None else antibiotic_hill
    cost = c.resistance_cost_max if resistance_cost_max is None else resistance_cost_max
    X0 = np.full(n, biomass0 / n)
    r_eff = p["r"] * (1.0 - cost * p["resistance"])

    def rhs(_t, y):
        X = np.maximum(y[:n], 0.0)
        drug = max(y[n], 0.0)
        kill = antibiotic_kill(drug, p["IC50"], kmax, hill)
        dX = X * (r_eff + p["A"] @ X - kill)
        return np.r_[dX, -kdeg * drug]

    t = np.linspace(0.0, c.t_end, c.n_times)
    sol = solve_ivp(rhs, (0.0, c.t_end), np.r_[X0, A0], t_eval=t,
                    method="LSODA", rtol=2e-6, atol=1e-10)
    if not sol.success:
        raise RuntimeError(sol.message)
    X = np.maximum(sol.y[:n].T, 0.0)
    drug = np.maximum(sol.y[n], 0.0)
    cutoff = max(biomass0 / n * c.extinction_fraction, 1e-14)
    # Numerical values below this threshold are treated as biologically absent.
    X[X < cutoff] = 0.0
    return {"t": sol.t, "X": X, "drug": drug, "cutoff": cutoff,
            "r_eff": r_eff, "kill": antibiotic_kill(drug[-1], p["IC50"], kmax, hill)}


def shannon(x):
    x = np.asarray(x, dtype=float)
    total = x.sum()
    if total <= 0:
        return 0.0
    q = x[x > 0] / total
    return float(-np.sum(q * np.log(q)))


def metrics(c, community, result, *, initial_biomass=None, dose=None, decay=None,
            kmax=None, hill=None, cost=None):
    X = result["X"]
    final = X[-1]
    base = {"dose": c.antibiotic_initial if dose is None else dose,
            "decay": c.antibiotic_decay if decay is None else decay,
            "kmax": c.antibiotic_kmax if kmax is None else kmax,
            "hill": c.antibiotic_hill if hill is None else hill,
            "resistance_cost": c.resistance_cost_max if cost is None else cost,
            "initial_biomass": c.initial_biomass if initial_biomass is None else initial_biomass,
            "final_biomass": float(final.sum()),
            "survivors": int(np.count_nonzero(final > result["cutoff"])),
            "final_shannon": shannon(final),
            "final_weighted_IC50": float(np.average(community["IC50"], weights=final)) if final.sum() else np.nan,
            "final_weighted_resistance": float(np.average(community["resistance"], weights=final)) if final.sum() else np.nan,
            "final_drug": float(result["drug"][-1]),
            "min_total_biomass": float(X.sum(axis=1).min())}
    return base


def species_outcomes(community, result, label):
    final = result["X"][-1]
    initial = result["X"][0]
    return community["species"].assign(
        scenario=label, initial_abundance=initial, final_abundance=final,
        log10_final_abundance=np.log10(np.maximum(final, 1e-300)),
        survived=final > result["cutoff"],
        relative_change=final / np.maximum(initial, 1e-300))


def stability_diagnostics(community):
    WA = community["W"] @ community["A"]
    sym = WA + WA.T
    return {"max_eigenvalue_symmetric_part": float(np.linalg.eigvalsh(sym).max()),
            "min_eigenvalue_P": float(np.linalg.eigvalsh(community["P"]).min()),
            "max_abs_diagonal_error": float(np.max(np.abs(np.diag(community["A"]) + 1))),
            "max_equilibrium_residual": float(np.max(np.abs(community["r"] + community["A"] @ community["x_star"])))}
