"""
Reproducibly generate the baseline 25-species, 3-resource community.

All random draws come from a single numpy Generator seeded with `seed`, so the
same seed always produces the same community. Every distribution below is a
PLACEHOLDER chosen to give a stable, interpretable toy ecosystem; replace with
measured values when available (see README).

Usage:
    python generate_community.py                  # seed 42 -> community.json (+ CSVs)
    python generate_community.py --seed 7 --out my_community.json
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

from consumer_resource_model import (Community, EcologicalParams, FutureTraits,
                                     save_community_json, validate_params)

N_SPECIES = 25
N_RESOURCES = 3


def generate_community(seed: int = 42,
                       n_species: int = N_SPECIES,
                       n_resources: int = N_RESOURCES,
                       *,
                       mu_range=(0.4, 1.0),
                       mortality_range=(0.05, 0.15),
                       self_limitation: float = 0.02,
                       preference_alpha: float = 0.6,
                       k_median: float = 0.5,
                       k_sigma: float = 0.5,
                       efficiency_range=(0.2, 0.5),
                       leakage_range=(0.1, 0.3),
                       production_alpha: float = 0.5,
                       supply=(2.0, 0.2, 0.2),
                       decay=(0.05, 0.05, 0.05),
                       init_log10_range=(-2.0, -1.0),
                       init_resources=(1.0, 0.1, 0.1),
                       susceptibility_beta=(2.0, 2.0)) -> Community:
    """Build a Community from a seed. Keyword arguments expose every placeholder."""
    rng = np.random.default_rng(seed)
    S, M = n_species, n_resources

    # Draw order is fixed so that adding new traits later (appended at the end)
    # does not change the values of existing ones for a given seed.
    mu_max = rng.uniform(*mu_range, size=S)
    mortality = rng.uniform(*mortality_range, size=S)
    preference = rng.dirichlet(np.full(M, preference_alpha), size=S)          # rows sum to 1
    half_saturation = k_median * rng.lognormal(0.0, k_sigma, size=(S, M))
    efficiency = rng.uniform(*efficiency_range, size=(S, M))
    leakage = rng.uniform(*leakage_range, size=S)
    production = rng.dirichlet(np.full(M, production_alpha), size=S)          # rows sum to 1
    initial_abundance = 10.0 ** rng.uniform(*init_log10_range, size=S)
    susceptibility = rng.beta(*susceptibility_beta, size=S)                   # bounded in [0, 1]

    eco = EcologicalParams(
        species_names=[f"sp{i + 1:02d}" for i in range(S)],
        resource_names=[f"R{j + 1}" for j in range(M)],
        mu_max=mu_max,
        mortality=mortality,
        self_limitation=np.full(S, self_limitation),
        leakage=leakage,
        preference=preference,
        half_saturation=half_saturation,
        efficiency=efficiency,
        production=production,
        supply=np.asarray(supply, dtype=float),
        decay=np.asarray(decay, dtype=float),
        initial_abundance=initial_abundance,
        initial_resources=np.asarray(init_resources, dtype=float),
    )
    community = Community(
        ecology=eco,
        future=FutureTraits(susceptibility=susceptibility),
        seed=seed,
        metadata={
            "description": "Baseline no-antibiotic consumer-resource community (placeholder parameters)",
            "susceptibility_distribution": f"Beta{tuple(susceptibility_beta)}",
            "susceptibility_convention": "0 = strongly resistant, 1 = highly susceptible; unused in baseline",
        },
    )
    validate_params(community)
    return community


def save_community_csvs(community: Community, prefix: str) -> None:
    """Write human-readable CSV tables next to the JSON file."""
    p = community.ecology
    species = pd.DataFrame({
        "species": p.species_names,
        "mu_max": p.mu_max,
        "mortality": p.mortality,
        "self_limitation": p.self_limitation,
        "leakage": p.leakage,
        "initial_abundance": p.initial_abundance,
        "susceptibility": community.future.susceptibility,
    })
    species.to_csv(f"{prefix}_species.csv", index=False)
    for name in ("preference", "half_saturation", "efficiency", "production"):
        df = pd.DataFrame(getattr(p, name), index=p.species_names, columns=p.resource_names)
        df.index.name = "species"
        df.to_csv(f"{prefix}_{name}.csv")
    pd.DataFrame({"resource": p.resource_names, "supply": p.supply, "decay": p.decay,
                  "initial_concentration": p.initial_resources}).to_csv(
        f"{prefix}_resources.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="community.json")
    args = ap.parse_args()

    community = generate_community(seed=args.seed)
    save_community_json(community, args.out)
    prefix = os.path.splitext(args.out)[0]
    save_community_csvs(community, prefix)
    print(f"Saved {args.out} and CSV tables ({prefix}_*.csv); seed={args.seed}")


if __name__ == "__main__":
    main()
