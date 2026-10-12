"""
kill_protection.py -- killing rate (gamma) against private protection of R.

Definitions used in notebook 10
-------------------------------
killing rate        gamma, the maximum lysis coefficient (y-axis). The per-capita
                    lysis rate of S is  gamma * H(a) * g(s).
protection          p = 1 - beta, the fraction of that lysis the resistant strain
                    avoids with the inhibitor present (x-axis). The per-capita
                    lysis rate of R is  (1 - p) * gamma * H(a) * g(s).
                    p is set through beta_min:  p = (1 - beta_min) * (1 - c * iota).

Exact identity from the model equations (S1, S2):

    d ln(n_r / n_s) / dt = g(s) * [ gamma * H(a) * p - (1 - alpha) ]

so while the antibiotic saturates (H -> 1) R gains exactly when gamma * p > 1 - alpha.
The neutral line is therefore p* = (1 - alpha) / gamma: a straight line against
1 / gamma, a hyperbola against gamma.
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

from ..model.ma2024_model import DEFAULT_PARAMS, BASE_I, hill, growth_rate, inhibitor_effect, SI_INITIAL_CONDITIONS
from .sweep_engine import integrate, endpoints, no_drug_reference

IOTA = inhibitor_effect(BASE_I, DEFAULT_PARAMS)
P_MAX = 1.0 - DEFAULT_PARAMS.c * IOTA                    # largest protection reachable at base c (beta_min = 0)
BASE_PROTECTION = (1.0 - DEFAULT_PARAMS.beta_min) * P_MAX

GAMMA_AXIS = np.sort(np.append(np.round(np.linspace(1.1, 3.1, 21), 10), DEFAULT_PARAMS.gamma))
PROTECTION_AXIS = np.sort(np.append(np.round(np.linspace(0.0, 0.30, 21), 10), round(BASE_PROTECTION, 10)))
ALPHAS = (0.85, 0.90, 0.95)
DOSES = (2.15, 100.0)


def beta_min_for(protection: float) -> float:
    """beta_min that gives the requested protection p = 1 - beta at base c and i."""
    return 1.0 - protection / P_MAX


def _params(gamma, protection, alpha):
    return replace(DEFAULT_PARAMS, gamma=float(gamma), alpha=float(alpha), beta_min=float(beta_min_for(protection)))


def _run(args):
    gamma, protection, alpha, a0 = args
    p = _params(gamma, protection, alpha)
    tr = integrate(p, a0, BASE_I)
    row = endpoints(tr, no_drug_reference())
    g = growth_rate(tr.s)
    lysis_s = gamma * hill(tr.a, p.h_a) * g
    g0 = growth_rate(SI_INITIAL_CONDITIONS["s0"])
    k0 = gamma * float(hill(a0, p.h_a)) * g0
    return {"gamma": gamma, "protection": protection, "alpha": alpha, "a0": a0, "beta_min": p.beta_min,
            "lysis0_S": k0, "lysis0_R": (1 - protection) * k0, "protection_rate0": protection * k0,
            "peak_death_S": float(np.max(lysis_s - g)),
            "peak_death_R": float(np.max((1 - protection) * lysis_s - alpha * g)),
            "log_kill_S": math.log(SI_INITIAL_CONDITIONS["ns0"] / max(row["ns_min"], 1e-300)),
            "log_kill_R": math.log(SI_INITIAL_CONDITIONS["nr0"] / max(row["nr_min"], 1e-300)),
            **row}


def run_plane(n_workers: int | None = None) -> pd.DataFrame:
    """gamma x protection for each alpha and dose (22 x 22 x 3 x 2 = 2,904 runs)."""
    tasks = [(float(g), float(q), a, d) for d in DOSES for a in ALPHAS
             for g, q in itertools.product(GAMMA_AXIS, PROTECTION_AXIS)]
    with mp.Pool(n_workers or os.cpu_count() or 1) as pool:
        rows = pool.map(_run, tasks, chunksize=16)
    return pd.DataFrame(rows)


def _delta(protection, gamma, alpha, a0, readout):
    p = _params(gamma, protection, alpha)
    if readout == "final":
        tr = integrate(p, a0, BASE_I)
    else:
        tr = integrate(p, a0, BASE_I, stages=((float(readout), 5.0),), adaptive=False)
    return float((tr.ln_nr[-1] - tr.ln_ns[-1]) - (tr.ln_nr[0] - tr.ln_ns[0]))


def _neutral(args):
    gamma, alpha, a0, readout = args
    lo, hi = 0.0, P_MAX
    f_lo, f_hi = _delta(lo, gamma, alpha, a0, readout), _delta(hi, gamma, alpha, a0, readout)
    p_star = brentq(_delta, lo, hi, args=(gamma, alpha, a0, readout), xtol=1e-7) if f_lo * f_hi < 0 else float("nan")
    return {"gamma": gamma, "alpha": alpha, "a0": a0, "readout": str(readout), "p_neutral": p_star,
            "p_criterion": (1 - alpha) / gamma}


def run_neutral_lines(n_workers: int | None = None) -> pd.DataFrame:
    """For every gamma: the protection at which R neither gains nor loses share,
    found by root-finding on the simulation. Read-outs: steady state for both
    doses, and t = 100 (still under treatment) for the high dose."""
    tasks = [(float(g), a, d, r) for (d, r) in ((2.15, "final"), (100.0, "final"), (100.0, 100))
             for a in ALPHAS for g in GAMMA_AXIS]
    with mp.Pool(n_workers or os.cpu_count() or 1) as pool:
        rows = pool.map(_neutral, tasks, chunksize=4)
    return pd.DataFrame(rows)


def linear_fit(x, y) -> dict:
    """Least-squares line y = slope * x + intercept, with R^2 and the largest
    deviation from the line as a fraction of the range of y."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3 or np.ptp(y) == 0:
        return {"n": int(len(x)), "slope": np.nan, "intercept": np.nan, "r2": np.nan, "max_dev": np.nan}
    slope, intercept = np.polyfit(x, y, 1)
    res = y - (slope * x + intercept)
    return {"n": int(len(x)), "slope": float(slope), "intercept": float(intercept),
            "r2": float(1 - res.var() / y.var()), "max_dev": float(np.abs(res).max() / np.ptp(y))}


def regrowth_prediction(gamma, protection, alpha, a0) -> dict:
    """Exact integrated form of the identity over a whole run.

    Integrating S1 gives  integral(g * gamma * H) = G - dlnS,  with
    G = integral of g(s) dt (total growth) and dlnS = ln(n_s final / n_s start).
    R neither gains nor loses over the run when  p * (G - dlnS) = (1 - alpha) * G, so

        p_neutral = (1 - alpha) / (1 - dlnS / G)

    When killing and regrowth are large, G >> dlnS and p_neutral -> 1 - alpha,
    whatever gamma is."""
    tr = integrate(_params(gamma, protection, alpha), a0, BASE_I)
    G = float(np.trapezoid(growth_rate(tr.s), tr.t))
    dlns = float(tr.ln_ns[-1] - tr.ln_ns[0])
    return {"G": G, "dlnS": dlns, "p_predicted": (1 - alpha) / (1 - dlns / G)}
