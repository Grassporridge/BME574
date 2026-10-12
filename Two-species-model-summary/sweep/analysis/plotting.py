"""
plotting.py -- shared figure helpers so every notebook draws the same kind of
figure the same way.

Colour rules: one perceptually uniform ramp (viridis) for every heat map
(signed quantities get a white dashed line at zero; outcome classes are steps
of the same ramp), and the colour-blind-safe Okabe-Ito palette for lines.
The red cross marks the base point and is drawn only in the panel that
contains it.
"""

from __future__ import annotations
import itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

from .grid_runner import BASE_VALUES, LOG_AXES
from .sweep_engine import OUTCOME_FLOOR_ORDER

OKABE_ITO = {"blue": "#0072B2", "orange": "#E69F00", "vermillion": "#D55E00",
             "green": "#009E73", "sky": "#56B4E9", "purple": "#CC79A7",
             "yellow": "#F0E442", "black": "#000000", "grey": "#999999"}

def _viridis(x: float) -> str:
    from matplotlib import colormaps
    from matplotlib.colors import to_hex
    return to_hex(colormaps["viridis"](x))


# One colour scale for every heat map: outcomes are steps of viridis, ordered
# from worst (dark purple) to the goal "R eradicated, S survives" (yellow).
OUTCOME_COLORS = {
    "community eradicated": _viridis(0.0),
    "community collapsed": _viridis(0.0),
    "S eradicated, R survives": _viridis(0.2),
    "S not viable, R persists": _viridis(0.38),
    "R eradicated, S not viable": _viridis(0.48),
    "S viable, R persists": _viridis(0.62),
    "R suppressed": _viridis(0.8),
    "R cleared": _viridis(0.9),
    "R eradicated, S survives": _viridis(1.0),
    "not converged": OKABE_ITO["grey"],
}

LABELS = {
    "alpha": "α (growth of R relative to S)", "beta_min": "β_min (1 = no private protection)",
    "c": "c (inhibitor efficacy inside R)", "gamma": "γ (maximum lysis)",
    "h_a": "h_a (antibiotic Hill coefficient)", "h_i": "h_i (inhibitor Hill coefficient)",
    "kappa_b": "κ_b (degradation by free enzyme)", "d_a": "d_a (basal antibiotic decay)",
    "d_b": "d_b (enzyme inactivation by inhibitor)", "xi": "ξ (nutrient recycling)",
    "phi_max": "φ_max (degradation by living R)", "a0": "a0 (antibiotic dose)",
    "i": "i (inhibitor dose)",
    "f_R_final": "Final resistant fraction f_R", "diversity_final": "Final Shannon diversity H",
    "diversity_min": "Lowest Shannon diversity H", "peak_loss_S": "Peak S loss",
    "peak_loss_total": "Peak total loss", "t_steady_state": "Time to steady state",
    "t_regain_S": "S recovery time (back to start)", "t_90_final_S": "Time for S to reach 90% of final",
    "total_final": "Final total density", "ns_final": "Final S density", "nr_final": "Final R density",
    "delta_log_ratio": "Change in ln(R/S)", "log10_f_R": "log10 final f_R",
    "log10_ns_min": "log10 lowest S density", "log10_nr_min": "log10 lowest R density",
}
SHORT = {"alpha": "α", "beta_min": "β_min", "c": "c", "gamma": "γ", "h_a": "h_a", "h_i": "h_i",
         "kappa_b": "κ_b", "d_a": "d_a", "d_b": "d_b", "xi": "ξ", "phi_max": "φ_max",
         "a0": "a0", "i": "i"}


def label(name: str) -> str:
    return LABELS.get(name, name)


def _coord(name, values):
    values = np.asarray(values, float)
    return np.log10(values) if name in LOG_AXES else values


def _format_axis(ax, name, which):
    if name in LOG_AXES:
        lo, hi = (0.1, 10) if name == "i" else (1, 100)
        ticks = [t for t in (0.1, 1, 10, 100) if lo <= t <= hi]
        getattr(ax, f"set_{which}ticks")(np.log10(ticks))
        getattr(ax, f"set_{which}ticklabels")([f"{t:g}" for t in ticks])
    getattr(ax, f"set_{which}label")(label(name))


def nearest(df: pd.DataFrame, name: str, target: float) -> float:
    values = np.sort(df[name].unique())
    return float(values[np.argmin(np.abs(values - target))])


def select(df: pd.DataFrame, fixed: dict) -> pd.DataFrame:
    """Rows where each parameter in `fixed` equals the nearest grid value."""
    out = df
    for name, target in fixed.items():
        out = out[np.isclose(out[name], nearest(df, name, target))]
    return out


def _draw(ax, sub, x, y, value, cmap, norm=None, vmin=None, vmax=None):
    piv = sub.pivot(index=y, columns=x, values=value)
    return ax.pcolormesh(_coord(x, piv.columns), _coord(y, piv.index), piv.values,
                         cmap=cmap, norm=norm, vmin=None if norm else vmin,
                         vmax=None if norm else vmax, shading="nearest")


def _zero_line(ax, sub, x, y, value):
    """White dashed line where a signed quantity changes sign."""
    piv = sub.pivot(index=y, columns=x, values=value)
    if np.nanmin(piv.values) < 0 < np.nanmax(piv.values):
        ax.contour(_coord(x, piv.columns), _coord(y, piv.index), piv.values, levels=[0.0],
                   colors="white", linewidths=1.5, linestyles="--")


def heatmap_slices(df, x, y, value, slice_by, slice_values=None, fixed=None, cmap="viridis",
                   vmin=None, vmax=None, diverging=False, title=None, mark_base=True):
    """`value` on an (x, y) plane, one panel per value of `slice_by`."""
    fixed = fixed or {}
    if slice_values is None:
        vals = np.sort(df[slice_by].unique())
        slice_values = [vals[0], BASE_VALUES[slice_by], vals[-1]]
    slice_values = sorted({nearest(df, slice_by, v) for v in slice_values})
    data = select(df, fixed)
    shown = data[np.isin(data[slice_by], slice_values)][value]
    vmin = np.nanmin(shown) if vmin is None else vmin
    vmax = np.nanmax(shown) if vmax is None else vmax
    norm = None
    fig, axes = plt.subplots(1, len(slice_values), figsize=(3.6 * len(slice_values) + 1.2, 3.6),
                             sharey=True, constrained_layout=True, squeeze=False)
    for ax, sv in zip(axes[0], slice_values):
        sub = data[np.isclose(data[slice_by], sv)]
        mesh = _draw(ax, sub, x, y, value, cmap, norm, vmin, vmax)
        if diverging:
            _zero_line(ax, sub, x, y, value)
        tag = " (base)" if np.isclose(sv, BASE_VALUES[slice_by]) else ""
        ax.set_title(f"{SHORT[slice_by]} = {sv:g}{tag}", fontsize=10)
        _format_axis(ax, x, "x")
        on_base_slice = bool(np.isclose(sv, BASE_VALUES[slice_by])) and all(
            np.isclose(nearest(df, k, v), BASE_VALUES[k]) for k, v in fixed.items())
        if mark_base and on_base_slice:          # the cross is the base point, so only in the base panel
            ax.plot(*_coord(x, [BASE_VALUES[x]]), *_coord(y, [BASE_VALUES[y]]), marker="x",
                    color="red", ms=7, mew=2)
    _format_axis(axes[0][0], y, "y")
    fig.colorbar(mesh, ax=axes, shrink=0.9, pad=0.02).set_label(label(value))
    extra = ", ".join(f"{SHORT[k]} = {nearest(df, k, v):g}" for k, v in fixed.items())
    fig.suptitle((title or label(value)) + (f"   [{extra}]" if extra else ""), fontsize=11)
    return fig


def pair_heatmaps(df, params, value, cmap="viridis", vmin=None, vmax=None, diverging=False):
    """`value` on every pair of `params`, the remaining ones held at base."""
    pairs = list(itertools.combinations(params, 2))
    ncol = min(3, len(pairs))
    nrow = int(np.ceil(len(pairs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 3.7 * nrow), constrained_layout=True,
                             squeeze=False)
    for ax in axes.flat[len(pairs):]:
        ax.axis("off")
    for ax, (x, y) in zip(axes.flat, pairs):
        fixed = {k: BASE_VALUES[k] for k in params if k not in (x, y)}
        sub = select(df, fixed)
        lo = np.nanmin(sub[value]) if vmin is None else vmin
        hi = np.nanmax(sub[value]) if vmax is None else vmax
        mesh = _draw(ax, sub, x, y, value, cmap, None, lo, hi)
        if diverging:
            _zero_line(ax, sub, x, y, value)
        _format_axis(ax, x, "x")
        _format_axis(ax, y, "y")
        ax.plot(*_coord(x, [BASE_VALUES[x]]), *_coord(y, [BASE_VALUES[y]]), marker="x",
                color="red", ms=7, mew=2)
        ax.set_title(", ".join(f"{SHORT[k]} = {v:g}" for k, v in fixed.items()) or "all others at base",
                     fontsize=9)
        fig.colorbar(mesh, ax=ax, shrink=0.9)
    fig.suptitle(f"{label(value)} (red x = base point)", fontsize=11)
    return fig


def marginal_effects(df, params, values, ncols=None):
    """For each parameter value: median of each endpoint over the rest of the
    grid (line) with the 10-90% range (band), and the base-slice value (dashed)."""
    fig, axes = plt.subplots(len(values), len(params), figsize=(3.4 * len(params), 2.5 * len(values)),
                             constrained_layout=True, squeeze=False, sharex="col", sharey="row")
    for r, value in enumerate(values):
        for c_, name in enumerate(params):
            ax = axes[r][c_]
            g = df.groupby(name)[value]
            x = _coord(name, g.median().index)
            ax.fill_between(x, g.quantile(0.1), g.quantile(0.9), color=OKABE_ITO["blue"], alpha=0.18, lw=0)
            ax.plot(x, g.median(), color=OKABE_ITO["blue"], lw=2)
            base = select(df, {k: BASE_VALUES[k] for k in params if k != name}).sort_values(name)
            ax.plot(_coord(name, base[name]), base[value], color="black", lw=1.2, ls="--")
            if r == len(values) - 1:
                _format_axis(ax, name, "x")
            if c_ == 0:
                ax.set_ylabel(label(value), fontsize=8)
            ax.grid(alpha=0.2)
    axes[0][0].legend(handles=[plt.Line2D([], [], color=OKABE_ITO["blue"], lw=2, label="median over grid"),
                               Patch(color=OKABE_ITO["blue"], alpha=0.18, label="10-90% of grid"),
                               plt.Line2D([], [], color="black", ls="--", lw=1.2, label="others at base")],
                      fontsize=7, frameon=False)
    return fig


def outcome_map(df, x, y, slice_by, slice_values=None, column="outcome_floor", fixed=None, title=None):
    """Categorical outcome on an (x, y) plane, one panel per `slice_by` value."""
    fixed = fixed or {}
    if slice_values is None:
        vals = np.sort(df[slice_by].unique())
        slice_values = [vals[0], BASE_VALUES[slice_by], vals[-1]]
    slice_values = sorted({nearest(df, slice_by, v) for v in slice_values})
    data = select(df, fixed)
    present = [o for o in OUTCOME_FLOOR_ORDER if o in set(data[np.isin(data[slice_by], slice_values)][column])]
    code = {o: k for k, o in enumerate(present)}
    cmap = ListedColormap([OUTCOME_COLORS[o] for o in present])
    norm = BoundaryNorm(np.arange(-0.5, len(present)), cmap.N)
    fig, axes = plt.subplots(1, len(slice_values), figsize=(3.6 * len(slice_values) + 0.6, 4.3),
                             sharey=True, constrained_layout=True, squeeze=False)
    for ax, sv in zip(axes[0], slice_values):
        sub = data[np.isclose(data[slice_by], sv)].copy()
        sub["_code"] = sub[column].map(code)
        _draw(ax, sub, x, y, "_code", cmap, norm)
        tag = " (base)" if np.isclose(sv, BASE_VALUES[slice_by]) else ""
        ax.set_title(f"{SHORT[slice_by]} = {sv:g}{tag}", fontsize=10)
        _format_axis(ax, x, "x")
    _format_axis(axes[0][0], y, "y")
    fig.legend(handles=[Patch(color=OUTCOME_COLORS[o], label=o) for o in present], loc="outside lower center",
               ncol=min(4, len(present)), fontsize=8, frameon=False)
    default = "Outcome with extinction floor" if column == "outcome_floor" else "Outcome (plain model)"
    fig.suptitle(title or default, fontsize=11)
    return fig


def outcome_counts(df, column="outcome_floor") -> pd.DataFrame:
    counts = df[column].value_counts()
    out = pd.DataFrame({"runs": counts, "% of grid": (counts / len(df) * 100).round(2)})
    return out.reindex([o for o in OUTCOME_FLOOR_ORDER if o in out.index])


def timecourses(cases, t_end=60.0, a0=None, i_dose=None, floor=False, log_y=True, title=None):
    """Plot S, R, antibiotic and resistant fraction for a list of (label, overrides) cases.
    `overrides` may contain Params fields and the doses 'a0' and 'i'."""
    from dataclasses import replace
    from ..model.ma2024_model import DEFAULT_PARAMS, BASE_A0, BASE_I
    from .sweep_engine import simulate, EXTINCTION_FLOOR
    fig, axes = plt.subplots(2, len(cases), figsize=(3.7 * len(cases), 5.2), sharex=True,
                             constrained_layout=True, squeeze=False)
    for j, (name, over) in enumerate(cases):
        over = dict(over)
        a = over.pop("a0", BASE_A0 if a0 is None else a0)
        i = over.pop("i", BASE_I if i_dose is None else i_dose)
        tr = simulate(replace(DEFAULT_PARAMS, **over), a, i, t_end=t_end, dt=t_end / 600,
                      floor=EXTINCTION_FLOOR if floor else None)
        ax = axes[0][j]
        ns, nr = np.where(tr.ns > 0, tr.ns, np.nan), np.where(tr.nr > 0, tr.nr, np.nan)
        ax.plot(tr.t, ns, color=OKABE_ITO["blue"], lw=2, label="S")
        ax.plot(tr.t, nr, color=OKABE_ITO["vermillion"], lw=2, ls="--", label="R")
        ax.plot(tr.t, np.maximum(tr.a, 1e-30), color=OKABE_ITO["grey"], lw=1.3, ls=":", label="antibiotic")
        if log_y:
            ax.set_yscale("log")
            lo = max(min(np.nanmin(ns), np.nanmin(nr)) / 3, 1e-12)
            ax.set_ylim(lo, 200)
            ax.axhline(EXTINCTION_FLOOR, color="black", lw=0.8, ls="-.")
        ax.set_title(name, fontsize=9)
        ax2 = axes[1][j]
        with np.errstate(invalid="ignore"):
            ax2.plot(tr.t, tr.nr / (tr.ns + tr.nr), color=OKABE_ITO["purple"], lw=2)
        ax2.axhline(0.5, color="0.6", lw=0.8, ls="--")
        ax2.set_ylim(-0.02, 1.02)
        ax2.set_xlabel("Time")
    axes[0][0].set_ylabel("Density / concentration")
    axes[1][0].set_ylabel("Resistant fraction f_R")
    axes[0][0].legend(fontsize=7, frameon=False)
    if title:
        fig.suptitle(title, fontsize=11)
    return fig
