"""
Multispecies consumer-resource model with cross-feeding (baseline, no antibiotics).

Notation
--------
S  number of species        (index i)
M  number of resources      (index j or k)

State vector  y = [N_1..N_S, R_1..R_M]

Species:
    dN_i/dt = N_i * ( growth_i(R) - mortality_i - self_limitation_i * N_i )
    growth_i(R) = sum_j mu_max_i * U_ij * R_j / (K_ij + R_j)

Resources:
    dR_j/dt = supply_j
              - sum_i N_i * uptake_ij(R)                (consumption)
              + sum_i N_i * leakage_i * P_ij * sum_k uptake_ik(R)   (cross-feeding)
              - decay_j * R_j

with per-biomass uptake
    uptake_ij(R) = (mu_max_i * U_ij / e_ij) * R_j / (K_ij + R_j)

so that  e_ij * uptake_ij = mu_max_i * U_ij * Monod_ij,  i.e. e_ij is the biomass
yield (biomass per unit resource) and mu_max_i is the genuine maximum growth rate
when every resource is saturating (because each row of U sums to 1).

Mass bookkeeping per unit resource taken up by species i from resource j:
    e_ij              -> new biomass
    leakage_i * P_ik  -> metabolite k released to the environment
    remainder         -> lost (respiration / waste), must be >= 0 (checked)

Future (not used here)
----------------------
`FutureTraits.susceptibility` is generated, stored and validated, but NOTHING in
`rhs` reads it. Antibiotic dynamics should be added as a separate parameter
dataclass and extra state variable(s), leaving EcologicalParams untouched.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp


# --------------------------------------------------------------------------- #
# Parameter containers
# --------------------------------------------------------------------------- #
@dataclass
class EcologicalParams:
    """All parameters of the no-antibiotic consumer-resource model."""

    species_names: list
    resource_names: list

    # --- species-level vectors, shape (S,) ---
    mu_max: np.ndarray            # maximum growth rate (1/time)
    mortality: np.ndarray         # density-independent loss rate (1/time)
    self_limitation: np.ndarray   # density-dependent loss coefficient (1/(biomass*time))
    leakage: np.ndarray           # fraction of uptake released as byproducts (0-1)

    # --- species x resource matrices, shape (S, M) ---
    preference: np.ndarray        # U_ij: resource-use preference, rows sum to 1
    half_saturation: np.ndarray   # K_ij: Monod half-saturation constant (resource units)
    efficiency: np.ndarray        # e_ij: biomass yield per unit resource consumed
    production: np.ndarray        # P_ij: fraction of species i's byproducts that is metabolite j, rows sum to 1

    # --- resource-level vectors, shape (M,) ---
    supply: np.ndarray            # external input rate (resource units/time)
    decay: np.ndarray             # first-order removal rate (1/time)

    # --- initial conditions ---
    initial_abundance: np.ndarray  # (S,)
    initial_resources: np.ndarray  # (M,)

    @property
    def n_species(self) -> int:
        return len(self.species_names)

    @property
    def n_resources(self) -> int:
        return len(self.resource_names)


@dataclass
class FutureTraits:
    """Traits generated now but NOT used by the baseline dynamics."""

    # 0 = strongly resistant, 1 = highly susceptible. Shape (S,).
    susceptibility: np.ndarray


@dataclass
class Community:
    ecology: EcologicalParams
    future: FutureTraits
    seed: Optional[int] = None
    metadata: dict = field(default_factory=dict)


@dataclass
class SimulationConfig:
    t_end: float = 2000.0
    n_points: int = 2001
    method: str = "LSODA"
    rtol: float = 1e-8
    atol: float = 1e-12
    extinction_threshold: float = 1e-6   # used by diagnostics, not by the solver


@dataclass
class SimulationResult:
    t: np.ndarray            # (T,)
    N: np.ndarray            # (T, S)
    R: np.ndarray            # (T, M)
    success: bool
    message: str

    @property
    def y(self) -> np.ndarray:
        return np.hstack([self.N, self.R])


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #
def validate_params(community: Community, tol: float = 1e-9) -> None:
    """Raise ValueError if shapes/ranges are inconsistent."""
    p = community.ecology
    S, M = p.n_species, p.n_resources

    shapes = {
        "mu_max": (S,), "mortality": (S,), "self_limitation": (S,), "leakage": (S,),
        "preference": (S, M), "half_saturation": (S, M),
        "efficiency": (S, M), "production": (S, M),
        "supply": (M,), "decay": (M,),
        "initial_abundance": (S,), "initial_resources": (M,),
    }
    for name, shape in shapes.items():
        arr = np.asarray(getattr(p, name))
        if arr.shape != shape:
            raise ValueError(f"{name} has shape {arr.shape}, expected {shape}")
    if np.asarray(community.future.susceptibility).shape != (S,):
        raise ValueError("susceptibility must have shape (S,)")

    for name in ("mu_max", "mortality", "self_limitation", "leakage", "preference",
                 "half_saturation", "efficiency", "production", "supply", "decay",
                 "initial_abundance", "initial_resources"):
        if np.any(np.asarray(getattr(p, name)) < 0):
            raise ValueError(f"{name} contains negative values")
    if np.any(p.half_saturation <= 0) or np.any(p.efficiency <= 0):
        raise ValueError("half_saturation and efficiency must be strictly positive")
    if not np.allclose(p.preference.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError("each row of preference must sum to 1")
    if not np.allclose(p.production.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError("each row of production must sum to 1")
    # Mass balance: biomass fraction + leaked fraction cannot exceed 1.
    if np.any(p.efficiency.max(axis=1) + p.leakage > 1.0 + tol):
        raise ValueError("efficiency + leakage exceeds 1 for some species (mass not conserved)")
    sus = np.asarray(community.future.susceptibility)
    if np.any(sus < 0) or np.any(sus > 1):
        raise ValueError("susceptibility must lie in [0, 1]")


# --------------------------------------------------------------------------- #
# Model terms (all matrix-based, no per-species code)
# --------------------------------------------------------------------------- #
def monod_saturation(R: np.ndarray, half_saturation: np.ndarray) -> np.ndarray:
    """Monod term R_j / (K_ij + R_j) for every species/resource pair -> (S, M).

    Negative R (solver round-off) is clipped to 0 so the term never goes negative
    or singular.
    """
    Rp = np.maximum(R, 0.0)[None, :]
    return Rp / (half_saturation + Rp)


def uptake_rates(R: np.ndarray, p: EcologicalParams) -> np.ndarray:
    """Per-biomass resource uptake, shape (S, M)."""
    vmax = p.mu_max[:, None] * p.preference / p.efficiency
    return vmax * monod_saturation(R, p.half_saturation)


def growth_rates(R: np.ndarray, p: EcologicalParams) -> np.ndarray:
    """Per-capita growth from resources, shape (S,):  sum_j mu_max_i U_ij Monod_ij."""
    return (p.mu_max[:, None] * p.preference * monod_saturation(R, p.half_saturation)).sum(axis=1)


def resource_consumption(N: np.ndarray, R: np.ndarray, p: EcologicalParams) -> np.ndarray:
    """Total consumption of each resource by the community, shape (M,)."""
    return N @ uptake_rates(R, p)


def crossfeeding_production(N: np.ndarray, R: np.ndarray, p: EcologicalParams) -> np.ndarray:
    """Total metabolite release into each resource pool, shape (M,)."""
    total_uptake = uptake_rates(R, p).sum(axis=1)       # (S,) per-biomass
    release = N * p.leakage * total_uptake              # (S,) metabolite flux per species
    return release @ p.production                       # (M,)


def crossfeeding_potential(p: EcologicalParams) -> np.ndarray:
    """Derived (S, S) matrix A[i, k]: how strongly species i's byproducts match
    species k's preferences (leakage_i * sum_m P_im U_km). For inspection only."""
    return p.leakage[:, None] * (p.production @ p.preference.T)


def species_derivatives(N: np.ndarray, R: np.ndarray, p: EcologicalParams) -> np.ndarray:
    net = growth_rates(R, p) - p.mortality - p.self_limitation * N
    return N * net


def resource_derivatives(N: np.ndarray, R: np.ndarray, p: EcologicalParams) -> np.ndarray:
    return (p.supply
            - resource_consumption(N, R, p)
            + crossfeeding_production(N, R, p)
            - p.decay * R)


def rhs(t: float, y: np.ndarray, p: EcologicalParams) -> np.ndarray:
    """ODE right-hand side for state y = [N, R]."""
    S = p.n_species
    N, R = y[:S], y[S:]
    return np.concatenate([species_derivatives(N, R, p), resource_derivatives(N, R, p)])


# --------------------------------------------------------------------------- #
# Simulation
# --------------------------------------------------------------------------- #
def simulate(community: Community, config: Optional[SimulationConfig] = None) -> SimulationResult:
    """Integrate the baseline (no-antibiotic) model."""
    config = config or SimulationConfig()
    validate_params(community)
    p = community.ecology

    y0 = np.concatenate([p.initial_abundance, p.initial_resources]).astype(float)
    t_eval = np.linspace(0.0, config.t_end, config.n_points)
    sol = solve_ivp(rhs, (0.0, config.t_end), y0, args=(p,), t_eval=t_eval,
                    method=config.method, rtol=config.rtol, atol=config.atol)

    S = p.n_species
    y = sol.y.T
    return SimulationResult(t=sol.t, N=y[:, :S], R=y[:, S:],
                            success=bool(sol.success), message=str(sol.message))


# --------------------------------------------------------------------------- #
# Diagnostics helpers
# --------------------------------------------------------------------------- #
def richness(N: np.ndarray, threshold: float = 1e-6) -> np.ndarray:
    """Number of species above threshold at each time point. N: (T, S)."""
    return (N > threshold).sum(axis=1)


def shannon_diversity(N: np.ndarray, threshold: float = 1e-6) -> np.ndarray:
    """Shannon index H = -sum p ln p over species above threshold, per time point."""
    Nc = np.where(N > threshold, N, 0.0)
    total = Nc.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        pr = np.where(total > 0, Nc / total, 0.0)
        logp = np.where(pr > 0, np.log(pr), 0.0)
    return -(pr * logp).sum(axis=1)


# --------------------------------------------------------------------------- #
# Serialization (JSON)
# --------------------------------------------------------------------------- #
def community_to_dict(c: Community) -> dict:
    def conv(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        return o
    return {
        "seed": c.seed,
        "metadata": c.metadata,
        "ecology": {k: conv(v) for k, v in asdict(c.ecology).items()},
        "future": {k: conv(v) for k, v in asdict(c.future).items()},
    }


def community_from_dict(d: dict) -> Community:
    eco = dict(d["ecology"])
    for k, v in eco.items():
        if k not in ("species_names", "resource_names"):
            eco[k] = np.asarray(v, dtype=float)
    fut = {k: np.asarray(v, dtype=float) for k, v in d["future"].items()}
    return Community(ecology=EcologicalParams(**eco), future=FutureTraits(**fut),
                     seed=d.get("seed"), metadata=d.get("metadata", {}))


def save_community_json(c: Community, path: str) -> None:
    with open(path, "w") as f:
        json.dump(community_to_dict(c), f, indent=2)


def load_community_json(path: str) -> Community:
    with open(path) as f:
        return community_from_dict(json.load(f))
