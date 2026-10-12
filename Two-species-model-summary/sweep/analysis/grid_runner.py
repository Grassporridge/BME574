"""
grid_runner.py -- one generic, parallel grid runner used by every sweep.

Axes are named by Params fields (alpha, beta_min, c, gamma, h_a, h_i, kappa_b,
d_a, d_b, xi, phi_max) or by the two dose variables "a0" and "i". Anything not
on an axis stays at the SI base value and the Fig. 2C dose (a0 = 2.15, i = 10).

`axis(name)` builds the standard axis for a parameter so that every sweep uses
the same rule: 21 evenly spaced values over the SI randomized range (log-spaced
for the two doses), plus the base value where it does not already fall on the
grid (gamma = 1.38 and a0 = 2.15). Every grid therefore contains the base point.
"""

from __future__ import annotations
import itertools
import os
import multiprocessing as mp

import numpy as np
import pandas as pd

from ..model.ma2024_model import (
    DEFAULT_PARAMS, PARAM_NAMES, SI_RANDOMIZED_RANGES, GROUP_CHOSEN_RANGES,
    SI_DOSE_RANGES, BASE_A0, BASE_I,
)
from .sweep_engine import compute_endpoints, no_drug_reference

N_PER_AXIS = 21
DOSE_KEYS = ("a0", "i")
LOG_AXES = ("a0", "i")
BASE_VALUES = {**{k: getattr(DEFAULT_PARAMS, k) for k in PARAM_NAMES}, "a0": BASE_A0, "i": BASE_I}
RANGES = {**SI_RANDOMIZED_RANGES, **GROUP_CHOSEN_RANGES, **SI_DOSE_RANGES}


def axis(name: str, n: int = N_PER_AXIS, include_base: bool = True) -> np.ndarray:
    """Standard sweep axis for one parameter (see module docstring)."""
    lo, hi = RANGES[name]
    values = np.geomspace(lo, hi, n) if name in LOG_AXES else np.linspace(lo, hi, n)
    values = np.round(values, 10)
    base = BASE_VALUES[name]
    if include_base and not np.any(np.isclose(values, base, rtol=1e-9, atol=1e-12)):
        values = np.sort(np.append(values, base))
    return values


def _worker(args):
    point, a0, i_dose, ns_ref, with_floor = args
    overrides = {k: v for k, v in point.items() if k not in DOSE_KEYS}
    row = compute_endpoints(overrides, a0=point.get("a0", a0), i_dose=point.get("i", i_dose),
                            ns_ref=ns_ref, with_floor=with_floor)
    return {**point, **row}


def run_points(points: list[dict], a0: float = BASE_A0, i_dose: float = BASE_I,
               n_workers: int | None = None, with_floor: bool = True) -> pd.DataFrame:
    """Run a list of parameter points (dicts) in parallel; one row per point."""
    ns_ref = no_drug_reference()
    tasks = [(pt, a0, i_dose, ns_ref, with_floor) for pt in points]
    n_workers = n_workers or os.cpu_count() or 1
    if n_workers == 1 or len(tasks) < 64:
        rows = [_worker(t) for t in tasks]
    else:
        with mp.Pool(n_workers) as pool:
            rows = pool.map(_worker, tasks, chunksize=max(8, len(tasks) // (n_workers * 40)))
    df = pd.DataFrame(rows)
    if "a0" not in df.columns:
        df["a0"] = a0
    if "i" not in df.columns:
        df["i"] = i_dose
    return df


def run_grid(axes: dict[str, np.ndarray], **kwargs) -> pd.DataFrame:
    """Full factorial grid over `axes` (name -> values)."""
    names = list(axes)
    points = [dict(zip(names, map(float, combo))) for combo in itertools.product(*axes.values())]
    return run_points(points, **kwargs)


def run_one_at_a_time(names=None, n: int = 41, **kwargs) -> pd.DataFrame:
    """Vary each parameter alone over its range (n values, base included) with
    everything else at base. Returns long format: parameter, value, endpoints."""
    names = list(RANGES) if names is None else list(names)
    points, labels = [], []
    for name in names:
        for v in axis(name, n):
            points.append({name: float(v)})
            labels.append((name, float(v)))
    df = run_points(points, **kwargs)
    df.insert(0, "parameter", [l[0] for l in labels])
    df.insert(1, "value", [l[1] for l in labels])
    return df[["parameter", "value"] + [c for c in df.columns
                                        if c not in names and c not in ("parameter", "value")]]


def run_global_random(n: int = 10000, seed: int = 574, **kwargs) -> pd.DataFrame:
    """Paper-style randomized strains: every strain parameter drawn uniformly
    from its range at once (Latin hypercube), dose fixed at the base dose."""
    from scipy.stats import qmc
    names = [k for k in RANGES if k not in DOSE_KEYS]
    unit = qmc.LatinHypercube(d=len(names), seed=seed).random(n)
    lo = [RANGES[k][0] for k in names]
    hi = [RANGES[k][1] for k in names]
    sample = qmc.scale(unit, lo, hi)
    points = [dict(zip(names, map(float, row))) for row in sample]
    return run_points(points, **kwargs)


TRAITS = ("alpha", "beta_min", "c")


def coarse_axis(name: str, step: int = 4) -> np.ndarray:
    """Every `step`-th value of the standard 21-point axis, plus the base value."""
    values = axis(name, include_base=False)[::step]
    base = BASE_VALUES[name]
    if not np.any(np.isclose(values, base)):
        values = np.sort(np.append(values, base))
    return values


def run_dose_by_trait(design: str = "single", **kwargs) -> pd.DataFrame:
    """The full dose plane (a0 x i, standard axes) repeated for many resistant strains.

    design = "single": each trait (alpha, beta_min, c) varied alone over its
                       standard 21-value axis, the other two at base.
    design = "joint":  all three traits varied together on a coarse grid
                       (every fourth standard value plus the base value).
    Returns one row per strain x dose; "varied" names the trait in the single design.
    """
    doses = list(itertools.product(map(float, axis("a0")), map(float, axis("i"))))
    base = {k: float(BASE_VALUES[k]) for k in TRAITS}
    strains, varied = [], []
    if design == "single":
        for name in TRAITS:
            for v in axis(name):
                strains.append({**base, name: float(v)})
                varied.append(name)
    elif design == "joint":
        for combo in itertools.product(*(coarse_axis(k) for k in TRAITS)):
            strains.append(dict(zip(TRAITS, map(float, combo))))
            varied.append("joint")
    else:
        raise ValueError("design must be 'single' or 'joint'")
    points, labels = [], []
    for strain, lab in zip(strains, varied):
        for a0, i_dose in doses:
            points.append({**strain, "a0": a0, "i": i_dose})
            labels.append(lab)
    df = run_points(points, **kwargs)
    df.insert(0, "varied", labels)
    return df


FINE_CORNER = {"phi_max": (0.0, 1.0), "kappa_b": (0.0, 0.35), "d_a": (0.0, 0.1)}


def fine_axis(name: str, n: int = N_PER_AXIS) -> np.ndarray:
    """21 evenly spaced values over the low-degradation corner (FINE_CORNER).
    The corner runs from the bottom of the standard range up to the base value
    for phi_max and kappa_b, so the base point is the top corner of the grid."""
    lo, hi = FINE_CORNER[name]
    return np.round(np.linspace(lo, hi, n), 10)


def run_fine_corner(**kwargs) -> pd.DataFrame:
    """phi_max x kappa_b x d_a on the fine low-degradation grid (xi at base)."""
    return run_grid({k: fine_axis(k) for k in FINE_CORNER}, **kwargs)
