"""
grid_runner.py

Parallel grid execution over (alpha, beta_min, c) for the two scenarios
described in sweep_engine.py:

  "treatment"  -- the paper's Fig. 2C dose (a0=2.15, i=10), where the six
                  requested endpoints (diversity, time-to-SS, f_R,
                  extinction, survival, recovery time, peak loss) are
                  computed via sweep_engine.compute_endpoints().
  "criterion"  -- highest doses in the paper's range (a0=100, i=10), where
                  sweep_engine's criterion_check() is run to validate the
                  paper's simplified selection criterion (Eq. 2) across a
                  dense structured grid (their Fig. 3C used 10,000 random
                  strains).

Default ranges are the Supplementary Information's randomized ranges:
alpha 0.75-1, beta_min 0-1, c 0-1.

Uses multiprocessing.Pool since this machine has 2 CPUs; each worker
process re-imports sweep_engine, which is import-cheap (no heavy state).
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count

from .sweep_engine import compute_endpoints, criterion_check, no_drug_reference


def _treatment_worker(args):
    alpha, beta_min, c, a0, i_dose, ns0_ref, t_span, n_points = args
    ep = compute_endpoints(alpha, beta_min, c, a0=a0, i_dose=i_dose,
                            ns0_no_drug=ns0_ref, t_span=t_span, n_points=n_points)
    return ep._asdict()


def run_treatment_grid(
    n_per_axis: int = 21,
    alpha_range=(0.75, 1.0),
    beta_min_range=(0.0, 1.0),
    c_range=(0.0, 1.0),
    a0: float = 2.15,
    i_dose: float = 10.0,
    t_span=(0.0, 150.0),
    n_points: int = 300,
    n_workers: int | None = None,
) -> pd.DataFrame:
    """Runs compute_endpoints() over a full n_per_axis^3 structured grid,
    in parallel. Returns a tidy DataFrame, one row per (alpha, beta_min, c)."""
    ns0_ref = no_drug_reference()

    alphas = np.linspace(*alpha_range, n_per_axis)
    beta_mins = np.linspace(*beta_min_range, n_per_axis)
    cs = np.linspace(*c_range, n_per_axis)

    tasks = [
        (alpha, beta_min, c, a0, i_dose, ns0_ref, t_span, n_points)
        for alpha in alphas for beta_min in beta_mins for c in cs
    ]

    n_workers = n_workers or max(1, cpu_count())
    with Pool(n_workers) as pool:
        results = pool.map(_treatment_worker, tasks, chunksize=32)

    df = pd.DataFrame(results)
    df.attrs["ns0_no_drug_reference"] = ns0_ref
    df.attrs["a0"] = a0
    df.attrs["i_dose"] = i_dose
    df.attrs["t_span"] = t_span
    return df


def _criterion_worker(args):
    alpha, beta_min, c, a0, i_dose, t_span, n_points = args
    return criterion_check(alpha, beta_min, c, a0=a0, i_dose=i_dose,
                            t_span=t_span, n_points=n_points)


def run_criterion_grid(
    n_per_axis: int = 41,
    alpha_range=(0.75, 1.0),
    beta_min_range=(0.0, 1.0),
    c_range=(0.0, 1.0),
    a0: float = 100.0,
    i_dose: float = 10.0,
    t_span=(0.0, 600.0),
    n_points: int = 50,
    n_workers: int | None = None,
) -> pd.DataFrame:
    """Runs criterion_check() over a full n_per_axis^3 structured grid, in
    parallel. Returns a tidy DataFrame, one row per (alpha, beta_min, c)."""
    alphas = np.linspace(*alpha_range, n_per_axis)
    beta_mins = np.linspace(*beta_min_range, n_per_axis)
    cs = np.linspace(*c_range, n_per_axis)

    tasks = [
        (alpha, beta_min, c, a0, i_dose, t_span, n_points)
        for alpha in alphas for beta_min in beta_mins for c in cs
    ]

    n_workers = n_workers or max(1, cpu_count())
    with Pool(n_workers) as pool:
        results = pool.map(_criterion_worker, tasks, chunksize=64)

    return pd.DataFrame(results)
