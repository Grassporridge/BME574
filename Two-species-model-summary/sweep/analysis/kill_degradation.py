"""
kill_degradation.py -- killing rate (gamma) against antibiotic degradation by
living resistant cells (phi_max).

gamma    maximum lysis coefficient: how hard the antibiotic kills (y-axis).
phi_max  maximum rate at which living R cells degrade the antibiotic (x-axis).
         With the inhibitor present the effective rate is phi = phi_max * (1 - c * iota).

Unlike private protection, phi_max does not appear in the equation for the
R/S ratio; it acts only through how long the antibiotic stays. So there is no
closed-form neutral line here. The test is empirical: for each gamma, find the
phi_max that holds an endpoint at a fixed level (an iso-line), and ask whether
phi_max on that line is a straight-line function of gamma.
"""

from __future__ import annotations
import itertools
import math
import multiprocessing as mp
import os
from dataclasses import replace

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from ..model.ma2024_model import DEFAULT_PARAMS, BASE_I, SI_INITIAL_CONDITIONS
from .sweep_engine import integrate, endpoints, no_drug_reference
from .kill_protection import linear_fit  # noqa: F401  (re-exported for the notebook)

GAMMA_AXIS = np.sort(np.append(np.round(np.linspace(1.1, 3.1, 21), 10), DEFAULT_PARAMS.gamma))
PHI_AXIS = np.round(np.linspace(0.0, 5.0, 21), 10)          # SI range, base value 1 is on the grid
DOSES = (2.15, 10.0)
PHI_SEARCH_MAX = 40.0                                         # iso-lines are followed beyond the SI range
LOSS_LEVELS = (0.5, 0.9, 0.99, 0.9999)                        # peak S loss levels for the iso-lines


def _row(gamma, phi_max, a0):
    p = replace(DEFAULT_PARAMS, gamma=float(gamma), phi_max=float(phi_max))
    row = endpoints(integrate(p, a0, BASE_I), no_drug_reference())
    row["log_kill_S"] = math.log(SI_INITIAL_CONDITIONS["ns0"] / max(row["ns_min"], 1e-300))
    row["log_kill_R"] = math.log(SI_INITIAL_CONDITIONS["nr0"] / max(row["nr_min"], 1e-300))
    return row


def _run(args):
    gamma, phi_max, a0 = args
    return {"gamma": gamma, "phi_max": phi_max, "a0": a0, **_row(gamma, phi_max, a0)}


def run_plane(n_workers: int | None = None) -> pd.DataFrame:
    """gamma x phi_max for each dose (22 x 21 x 2 = 924 runs)."""
    tasks = [(float(g), float(f), d) for d in DOSES for g, f in itertools.product(GAMMA_AXIS, PHI_AXIS)]
    with mp.Pool(n_workers or os.cpu_count() or 1) as pool:
        return pd.DataFrame(pool.map(_run, tasks, chunksize=8))


def _iso(args):
    gamma, a0, column, level = args
    f = lambda phi: _row(gamma, phi, a0)[column] - level
    lo, hi = f(0.0), f(PHI_SEARCH_MAX)
    phi = brentq(f, 0.0, PHI_SEARCH_MAX, xtol=1e-6) if lo * hi < 0 else float("nan")
    return {"gamma": gamma, "a0": a0, "endpoint": column, "level": level, "phi_iso": phi}


def run_iso_lines(base: dict, n_workers: int | None = None) -> pd.DataFrame:
    """For every gamma: the phi_max that holds an endpoint at a fixed level.
    Levels: the listed peak-S-loss levels, plus the base-point value of peak S
    loss and of the final resistant fraction (`base` maps dose -> those values)."""
    tasks = []
    for d in DOSES:
        targets = [("peak_loss_S", lv) for lv in LOSS_LEVELS]
        targets += [("peak_loss_S", base[d]["peak_loss_S"]), ("f_R_final", base[d]["f_R_final"])]
        tasks += [(float(g), d, col, float(lv)) for col, lv in targets for g in GAMMA_AXIS]
    with mp.Pool(n_workers or os.cpu_count() or 1) as pool:
        return pd.DataFrame(pool.map(_iso, tasks, chunksize=2))
