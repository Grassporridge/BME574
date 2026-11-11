"""
Run the no-antibiotic baseline simulation.

Usage:
    python run_baseline.py                          # generate seed-42 community if community.json is absent
    python run_baseline.py --community community.json --t-end 1000
    python run_baseline.py --seed 7 --regenerate

Outputs:
    baseline_timecourse.csv   columns: time, <species...>, <resources...>
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

from consumer_resource_model import (SimulationConfig, load_community_json, resource_derivatives,
                                     richness, save_community_json, shannon_diversity,
                                     simulate, species_derivatives)
from generate_community import generate_community, save_community_csvs


def load_or_generate(path: str, seed: int, regenerate: bool):
    if os.path.exists(path) and not regenerate:
        print(f"Loading community from {path}")
        return load_community_json(path)
    print(f"Generating community (seed={seed}) -> {path}")
    community = generate_community(seed=seed)
    save_community_json(community, path)
    save_community_csvs(community, os.path.splitext(path)[0])
    return community


def print_diagnostics(community, result, config: SimulationConfig) -> None:
    p = community.ecology
    thr = config.extinction_threshold
    N, R = result.N, result.R
    final = N[-1]

    print("\n=== Baseline diagnostics ===")
    print(f"Solver success : {result.success} ({result.message})")
    print(f"Time span      : 0 - {result.t[-1]:g}  ({len(result.t)} points)")
    print(f"Min value seen : N={N.min():.3e}, R={R.min():.3e}  (should be >= ~0)")
    print(f"Final biomass  : {final.sum():.4f}")
    print(f"Final richness : {int(richness(N[-1:], thr)[0])} / {p.n_species} (threshold {thr:g})")
    print(f"Final Shannon  : {shannon_diversity(N[-1:], thr)[0]:.3f}")
    print("Final resources: " + ", ".join(f"{n}={v:.4f}" for n, v in zip(p.resource_names, R[-1])))

    survivors = np.where(final > thr)[0]
    order = survivors[np.argsort(final[survivors])[::-1]]
    print("\nSurvivors (descending abundance):")
    for i in order:
        print(f"  {p.species_names[i]}  N={final[i]:.4f}  ({final[i] / final.sum():.1%})")
    extinct = [p.species_names[i] for i in range(p.n_species) if final[i] <= thr]
    print(f"\nExtinct/below threshold ({len(extinct)}): {', '.join(extinct) if extinct else 'none'}")

    # Steady-state check: largest |dy/dt| at the final state, and change over last 10% of run.
    dN = species_derivatives(final, R[-1], p)
    dR = resource_derivatives(final, R[-1], p)
    tail = slice(int(0.9 * len(result.t)), None)
    drift = np.abs(N[tail][-1] - N[tail][0]).max()
    print(f"\nMax |dN/dt| at end: {np.abs(dN).max():.3e}   Max |dR/dt| at end: {np.abs(dR).max():.3e}")
    print(f"Max abundance change over final 10% of run: {drift:.3e}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--community", default="community.json")
    ap.add_argument("--seed", type=int, default=42, help="seed used if the community is generated")
    ap.add_argument("--regenerate", action="store_true")
    ap.add_argument("--t-end", type=float, default=2000.0)
    ap.add_argument("--n-points", type=int, default=2001)
    ap.add_argument("--out", default="baseline_timecourse.csv")
    args = ap.parse_args()

    community = load_or_generate(args.community, args.seed, args.regenerate)
    config = SimulationConfig(t_end=args.t_end, n_points=args.n_points)
    result = simulate(community, config)

    p = community.ecology
    df = pd.DataFrame(np.hstack([result.t[:, None], result.N, result.R]),
                      columns=["time"] + list(p.species_names) + list(p.resource_names))
    df.to_csv(args.out, index=False)
    print(f"Saved time course -> {args.out}")
    print_diagnostics(community, result, config)


if __name__ == "__main__":
    main()
