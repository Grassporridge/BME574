"""
dose_window.py -- how wide is the antibiotic-dose window that removes R but keeps S?

For one strain and one inhibitor dose, raising the antibiotic dose a0 pushes the
lowest density reached by each strain down. Each strain is lost once its lowest
density falls below one cell (the extinction floor):

    a0_R = lowest antibiotic dose at which R drops below one cell
    a0_S = lowest antibiotic dose at which S drops below one cell
    window ratio = a0_S / a0_R

A ratio above 1 means there are doses that eradicate R and spare S; the ratio is
how wide that window is (1.2 = the dose can be up to 20% above a0_R). A ratio at
or below 1 means S is lost first, so no dose can meet the goal.

Both critical doses are interpolated (log-log) from the lowest densities in the
plain model on the standard a0 axis, so the ratio is a smooth measure that does
not depend on whether a grid dose happens to fall inside a narrow window. It is
checked against the explicit extinction-floor runs in notebook 08.

Open-ended windows: if R is eradicated on the dose axis but S never is (up to
the highest dose, a0 = 100), the true ratio is unknown but at least
(highest dose) / a0_R. That lower bound is reported and the case is flagged
`censored`.
"""

from __future__ import annotations
import numpy as np
import pandas as pd

from .sweep_engine import EXTINCTION_FLOOR

TRAITS = ["alpha", "beta_min", "c"]


def critical_dose(a0: np.ndarray, lowest_density: np.ndarray, floor: float = EXTINCTION_FLOOR) -> float:
    """Lowest a0 at which `lowest_density` crosses below the floor (log-log interpolation).
    NaN if it never crosses on the axis."""
    la = np.log10(a0)
    lm = np.log10(np.clip(lowest_density, 1e-300, None)) - np.log10(floor)
    if lm[0] <= 0:                      # already below the floor at the lowest dose
        return float(a0[0])
    idx = np.where((lm[:-1] > 0) & (lm[1:] <= 0))[0]
    if len(idx) == 0:
        return np.nan
    k = idx[0]
    return float(10 ** (la[k] + (la[k + 1] - la[k]) * lm[k] / (lm[k] - lm[k + 1])))


def windows(df: pd.DataFrame, by=TRAITS) -> pd.DataFrame:
    """One row per strain and inhibitor dose: critical doses, window ratio, whether
    the ratio is only a lower bound (`censored`), and the number of grid doses
    that met the goal in the extinction-floor runs."""
    rows = []
    for key, s in df.groupby(list(by) + ["i"]):
        s = s.sort_values("a0")
        a_r = critical_dose(s.a0.values, s.nr_min.values)
        a_s = critical_dose(s.a0.values, s.ns_min.values)
        censored = bool(np.isfinite(a_r) and not np.isfinite(a_s))
        ratio = float(s.a0.max()) / a_r if censored else a_s / a_r
        rows.append({**dict(zip(list(by) + ["i"], key)), "a0_R": a_r, "a0_S": a_s,
                     "window_ratio": ratio, "censored": censored,
                     "goal_doses": int(s.goal_met_floor.sum())})
    return pd.DataFrame(rows)


def strain_summary(df: pd.DataFrame, by=TRAITS) -> pd.DataFrame:
    """One row per strain: how often and how robustly the goal can be met.
    best_window_ratio is the widest window over all inhibitor doses and is a lower
    bound when best_is_lower_bound is True (open-ended window);
    widest_measured_ratio is the widest window for which both critical doses lie
    on the dose axis."""
    w = windows(df, by)
    out = []
    for key, s in w.groupby(list(by)):
        key = key if isinstance(key, tuple) else (key,)
        best = s.loc[s.window_ratio.idxmax()] if s.window_ratio.notna().any() else None
        out.append({**dict(zip(by, key)),
                    "goal_doses": int(s.goal_doses.sum()),
                    "inhibitor_doses_with_window": int((s.window_ratio > 1).sum()),
                    "widest_measured_ratio": float(s[~s.censored].window_ratio.max()) if (~s.censored).any() else np.nan,
                    "best_window_ratio": np.nan if best is None else float(best.window_ratio),
                    "best_is_lower_bound": False if best is None else bool(best.censored),
                    "i_at_best_window": np.nan if best is None else float(best.i),
                    "a0_R_at_best_window": np.nan if best is None else float(best.a0_R)})
    return pd.DataFrame(out)
