"""Expanded Person 3 analysis for the Ma et al. two-population model.

It uses the same five-state model as the shared implementation, but adds selection-
relative metrics, multi-dose response atlases, mechanistic ablations,
initial-composition analysis, robustness, and matched-inhibitor identifiability
checks.

Run locally with:

    conda run -n ar_cp python -m src.analysis.person3_expanded_analysis

Outputs are written under the repository root:

    results/figures/*.png and *.pdf
    results/tables/*.csv and *.json
    results/reports/Person_3_Expanded_Report.pdf
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

MPLCONFIGDIR = Path(tempfile.gettempdir()) / "person3_expanded_mplconfig"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.stats import spearmanr, qmc

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as mpl_font_manager
from matplotlib.colors import Normalize, LogNorm

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image,
    Table,
    TableStyle,
    PageBreak,
)


# DejaVu Sans is distributed with Matplotlib and provides the Greek letters,
# Unicode minus, and subscript characters used in the model primer.
pdfmetrics.registerFont(TTFont("DejaVuSans", mpl_font_manager.findfont("DejaVu Sans")))
pdfmetrics.registerFont(TTFont(
    "DejaVuSans-Bold",
    mpl_font_manager.findfont(mpl_font_manager.FontProperties(family="DejaVu Sans", weight="bold")),
))


# This module lives under src/analysis; write outputs to the repository root.
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
TABLES = RESULTS / "tables"
REPORTS = RESULTS / "reports"
for directory in (FIGURES, TABLES, REPORTS):
    directory.mkdir(parents=True, exist_ok=True)


BASE = {
    "alpha": 0.95,
    "beta_min": 0.90,
    "xi": 0.80,
    "d_a": 0.02,
    "kappa_b": 0.35,
    "d_b": 1.0,
    "gamma": 1.38,
    "h_a": 3.0,
    "i": 1.0,
    "h_i": 2.0,
    "phi_max": 1.0,
    "c": 0.70,
}

NS0 = 0.2
NR0 = 0.2
S0 = 4.0
A0 = 2.15
B0 = 0.0
T_ENDS = {2.15: 300.0, 10.0: 400.0, 100.0: 600.0}
N_POINTS = 301
EXTINCTION_THRESHOLD = 1e-4
RECOVERY_FRACTION = 0.90
SELECTION_TOLERANCE = 1e-3

I_VALUES = np.logspace(-1, 1, 21)
HI_VALUES = np.linspace(1.0, 5.0, 13)
DB_VALUES = np.linspace(1.0, 10.0, 13)
HI_DB_DOSES = [0.1, 1.0, 10.0]
ANCHORS = [2.15, 10.0, 100.0]


def hill(x: float | np.ndarray, h: float) -> float | np.ndarray:
    x = np.maximum(np.asarray(x), 0.0)
    xh = np.power(x, h)
    return xh / (1.0 + xh)


def inhibitor_activity(i: float, h_i: float) -> float:
    return float(hill(i, h_i))


def pathway_values(p: dict) -> tuple[float, float, float, float]:
    iota = inhibitor_activity(p["i"], p["h_i"])
    beta = p["beta_min"] + p["c"] * (1.0 - p["beta_min"]) * iota
    phi = p["phi_max"] * (1.0 - p["c"] * iota)
    return iota, beta, phi, p["d_b"] * iota


def pathway_for_mode(p: dict, mode: str) -> tuple[float, float, float, float]:
    """Return pathway values for a mechanistic ablation.

    The full model is the reference.  In each ablation, one mechanism follows
    the selected inhibitor setting while the remaining mechanisms are held at
    the baseline setting.  This is a causal bookkeeping comparison, not a new
    biological model.
    """
    full = pathway_values(p)
    ref = pathway_values(BASE)
    if mode == "full":
        return full
    if mode == "private_protection_only":
        return ref[0], full[1], ref[2], ref[3]
    if mode == "public_degradation_only":
        return ref[0], ref[1], full[2], ref[3]
    if mode == "free_bla_turnover_only":
        return ref[0], ref[1], ref[2], full[3]
    if mode == "inhibitor_independent_control":
        return ref
    raise ValueError(f"Unknown pathway mode: {mode}")


def rhs_logspace(t: float, z: np.ndarray, p: dict, mode: str = "full") -> np.ndarray:
    """RHS with log-transformed live populations.

    The nutrient, antibiotic, and free-Bla states remain in their original
    coordinates.  Logging n_s and n_r preserves the composition signal during
    high-dose bottlenecks without changing the model equations.
    """
    ns, nr = np.exp(np.clip(z[:2], -745.0, 700.0))
    s, a, b = np.maximum(z[2:], 0.0)
    g = s / (1.0 + s)
    drug_effect = float(hill(a, p["h_a"]))
    l = p["gamma"] * drug_effect * g
    iota, beta, phi, effective_db = pathway_for_mode(p, mode)
    dns = (g - l) * ns
    dnr = (p["alpha"] * g - beta * l) * nr
    ds = (p["xi"] * l - g) * ns + (p["xi"] * beta * l - p["alpha"] * g) * nr
    da = -p["kappa_b"] * b * a - phi * nr * a - p["d_a"] * a
    db = beta * l * nr - effective_db * b
    return np.array([
        g - l,
        p["alpha"] * g - beta * l,
        ds,
        da,
        db,
    ])


def simulate(
    p: dict,
    ns0: float = NS0,
    nr0: float = NR0,
    a0: float = A0,
    t_end: float | None = None,
    n_points: int = N_POINTS,
    mode: str = "full",
) -> dict:
    if t_end is None:
        t_end = T_ENDS.get(float(a0), 300.0)
    y0 = np.array([
        np.log(max(ns0, 1e-300)),
        np.log(max(nr0, 1e-300)),
        S0,
        a0,
        B0,
    ])
    t_eval = np.linspace(0.0, float(t_end), int(n_points))
    sol = solve_ivp(
        lambda t, z: rhs_logspace(t, z, p, mode=mode),
        (0.0, float(t_end)),
        y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=3e-7,
        atol=1e-9,
        max_step=2.0,
    )
    y = np.vstack([
        np.exp(np.clip(sol.y[0], -745.0, 700.0)),
        np.exp(np.clip(sol.y[1], -745.0, 700.0)),
        np.maximum(sol.y[2], 0.0),
        np.maximum(sol.y[3], 0.0),
        np.maximum(sol.y[4], 0.0),
    ])
    return {"t": sol.t, "y": y, "success": bool(sol.success), "message": sol.message}


def shannon(fr: np.ndarray | float) -> np.ndarray | float:
    fr = np.asarray(fr)
    fs = 1.0 - fr
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(fs > 0.0, -fs * np.log(fs), 0.0)
        out += np.where(fr > 0.0, -fr * np.log(fr), 0.0)
    return out


def trapz(y: np.ndarray, x: np.ndarray) -> float:
    fn = getattr(np, "trapezoid", None) or np.trapz
    return float(fn(y, x))


def time_to_steady(t: np.ndarray, y: np.ndarray) -> float:
    if len(t) < 20:
        return np.nan
    dt = np.diff(t)
    # Relative increments for live populations, absolute bounded increments for
    # environmental states. This avoids labeling a small absolute state change
    # as convergence while a population is still changing exponentially.
    relative_live = np.abs(np.diff(np.log(np.maximum(y[:2], 1e-300)), axis=1)) / dt
    environmental = np.abs(np.diff(y[2:], axis=1)) / dt
    changes = np.vstack([relative_live, environmental])
    score = np.max(changes, axis=0)
    window = max(8, len(t) // 25)
    threshold = 2e-4
    for j in range(0, len(score) - window):
        if np.max(score[j:j + window]) < threshold:
            return float(t[j])
    return np.nan


def endpoint_metrics(
    sim: dict,
    p: dict,
    reference: dict | None = None,
    ns0: float = NS0,
    nr0: float = NR0,
    extinction_threshold: float = EXTINCTION_THRESHOLD,
) -> dict:
    t, y = sim["t"], sim["y"]
    ns, nr, s, a, b = y
    total = ns + nr
    with np.errstate(divide="ignore", invalid="ignore"):
        fr = np.divide(nr, total, out=np.zeros_like(nr), where=total > 0)
        log_ratio = np.log(np.maximum(nr, 1e-300) / np.maximum(ns, 1e-300))
    H = shannon(fr)
    initial_total = ns0 + nr0
    initial_fr = nr0 / initial_total if initial_total > 0 else np.nan
    min_idx = int(np.argmin(total))
    ref_ns_final = np.nan
    if reference is not None and len(reference["y"][0]):
        ref_ns_final = float(reference["y"][0, -1])
    target_ns = RECOVERY_FRACTION * ref_ns_final if np.isfinite(ref_ns_final) else np.nan
    recovery_start = min_idx
    t_recovery = np.nan
    if np.isfinite(target_ns):
        candidates = np.flatnonzero(ns[recovery_start:] >= target_ns)
        if len(candidates):
            candidate = int(candidates[0] + recovery_start)
            # Require the threshold to remain satisfied through the final 10%
            # of the sampled trajectory, avoiding one-point false recoveries.
            tail = ns[candidate:]
            if len(tail) and np.all(tail[-max(3, len(tail) // 10):] >= target_ns):
                t_recovery = float(t[candidate])
    final_total = float(total[-1])
    final_ns = float(ns[-1])
    final_nr = float(nr[-1])
    final_fr = float(fr[-1])
    delta_fr = final_fr - initial_fr
    delta_log_ratio = float(log_ratio[-1] - log_ratio[0])
    if abs(delta_fr) <= SELECTION_TOLERANCE:
        selection_direction = "neutral"
    elif delta_fr > 0:
        selection_direction = "resistant-enriched"
    else:
        selection_direction = "susceptible-enriched"
    if final_total < extinction_threshold:
        outcome = "whole-community collapse"
    elif final_nr < extinction_threshold:
        outcome = "resistant extinction"
    elif final_ns < extinction_threshold:
        outcome = "susceptible extinction"
    elif final_fr <= 0.20:
        outcome = "susceptible-dominant"
    elif final_fr >= 0.80:
        outcome = "resistant-dominant"
    else:
        outcome = "mixed"
    iota, beta, phi, effective_db = pathway_values(p)
    return {
        "initial_fR": float(initial_fr),
        "initial_H": float(shannon(initial_fr)),
        "final_ns": final_ns,
        "final_nr": final_nr,
        "final_total": final_total,
        "final_fR": final_fr,
        "delta_fR": float(delta_fr),
        "final_log_ratio": float(log_ratio[-1]),
        "delta_log_ratio": delta_log_ratio,
        "final_H": float(H[-1]),
        "delta_H": float(H[-1] - H[0]),
        "min_total": float(total[min_idx]),
        "peak_loss": float(1.0 - total[min_idx] / initial_total),
        "peak_loss_time": float(t[min_idx]),
        "auc_total": trapz(total, t),
        "auc_a": trapz(a, t),
        "auc_b": trapz(b, t),
        "peak_a": float(np.max(a)),
        "peak_b": float(np.max(b)),
        "time_peak_b": float(t[int(np.argmax(b))]),
        "time_to_ss": float(time_to_steady(t, y)),
        "time_to_recovery": float(t_recovery) if np.isfinite(t_recovery) else np.nan,
        "resistant_extinct": bool(final_nr < extinction_threshold),
        "susceptible_extinct": bool(final_ns < extinction_threshold),
        "community_collapsed": bool(final_total < extinction_threshold),
        "susceptible_recovered": bool(np.isfinite(t_recovery)),
        "selection_direction": selection_direction,
        "outcome": outcome,
        "solver_success": bool(sim["success"]),
        "iota": iota,
        "beta": beta,
        "phi": phi,
        "effective_db_iota": effective_db,
    }


def no_treatment_reference(p: dict, ns0: float = NS0, nr0: float = NR0, a0: float = 0.0, t_end: float | None = None) -> dict:
    return simulate(p, ns0=ns0, nr0=nr0, a0=a0, t_end=t_end)


def record_with_parameters(metrics: dict, p: dict, **extra) -> dict:
    row = dict(metrics)
    row.update({key: float(value) for key, value in p.items()})
    row.update(extra)
    return row


def run_pair_grid(
    x_name: str,
    x_values: Iterable[float],
    y_name: str,
    y_values: Iterable[float],
    fixed: dict,
    a0: float,
    label: str,
    mode: str = "full",
) -> pd.DataFrame:
    records = []
    t_end = T_ENDS[float(a0)]
    reference = no_treatment_reference(fixed, a0=0.0, t_end=t_end)
    x_values = list(x_values)
    y_values = list(y_values)
    total = len(x_values) * len(y_values)
    count = 0
    for y_value in y_values:
        for x_value in x_values:
            p = dict(fixed)
            p[x_name] = float(x_value)
            p[y_name] = float(y_value)
            sim = simulate(p, a0=a0, t_end=t_end, mode=mode)
            metrics = endpoint_metrics(sim, p, reference=reference)
            records.append(record_with_parameters(
                metrics,
                p,
                a0=float(a0),
                plane=label,
                pathway_mode=mode,
                x_name=x_name,
                y_name=y_name,
                x_value=float(x_value),
                y_value=float(y_value),
            ))
            count += 1
        if count == total or count % max(1, len(x_values) * 4) == 0:
            print(f"  {label}, a0={a0:g}: {count}/{total}", flush=True)
    frame = pd.DataFrame(records)
    frame.to_csv(TABLES / f"atlas_a0_{a0:g}_{label}.csv", index=False)
    return frame


def savefig(fig: plt.Figure, name: str) -> None:
    fig.savefig(FIGURES / f"{name}.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def pathway_plot() -> None:
    doses = np.logspace(-1, 1, 300)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for hi in [1.0, 2.0, 5.0]:
        axes[0, 0].plot(doses, [inhibitor_activity(i, hi) for i in doses], label=f"h_i={hi:g}")
    axes[0, 0].axvline(1, color="0.5", ls="--", lw=1)
    axes[0, 0].set_xscale("log")
    axes[0, 0].set(xlabel="inhibitor dose i", ylabel="activity iota", title="Inhibitor Hill response")
    axes[0, 0].legend(frameon=False)
    for hi in [1.0, 2.0, 5.0]:
        betas, phis = [], []
        for i in doses:
            p = dict(BASE, i=float(i), h_i=hi)
            _, beta, phi, _ = pathway_values(p)
            betas.append(beta)
            phis.append(phi)
        axes[0, 1].plot(doses, betas, label=f"h_i={hi:g}")
        axes[1, 0].plot(doses, phis, label=f"h_i={hi:g}")
    for ax in [axes[0, 1], axes[1, 0]]:
        ax.set_xscale("log")
        ax.axvline(1, color="0.5", ls="--", lw=1)
        ax.legend(frameon=False, fontsize=8)
    axes[0, 1].set(xlabel="inhibitor dose i", ylabel="beta", title="Private protection is weakened")
    axes[1, 0].set(xlabel="inhibitor dose i", ylabel="phi", title="Public degradation is weakened")
    for hi in [1.0, 2.0, 5.0]:
        axes[1, 1].plot(doses, [BASE["d_b"] * inhibitor_activity(i, hi) for i in doses], label=f"h_i={hi:g}")
    axes[1, 1].set_xscale("log")
    axes[1, 1].axvline(1, color="0.5", ls="--", lw=1)
    axes[1, 1].set(xlabel="inhibitor dose i", ylabel="effective d_b*iota", title="Free-Bla inactivation rate")
    axes[1, 1].legend(frameon=False, fontsize=8)
    fig.suptitle("Person 3 pharmacodynamic pathway", y=1.02)
    fig.tight_layout()
    savefig(fig, "expanded_01_inhibitor_pathway")


def validation_plot() -> None:
    p = dict(BASE, i=10.0, h_i=2.0, d_b=1.0)
    sim = simulate(p, a0=2.15, t_end=40.0, n_points=401)
    t = sim["t"]
    ns, nr, s, a, b = sim["y"]
    fr = nr / np.maximum(ns + nr, 1e-300)
    iota, beta, phi, effective_db = pathway_values(p)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    axes[0].plot(t, ns, label="n_s", color="tab:green")
    axes[0].plot(t, nr, label="n_r", color="tab:gray")
    axes[0].plot(t, a, label="a", color="tab:red", ls="--")
    axes[0].plot(t, b, label="b", color="black", ls=":")
    axes[0].set(xlabel="time", ylabel="state", title="Ma-style validation trajectory")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(t, fr, color="tab:purple")
    axes[1].axhline(fr[0], color="0.5", ls="--", lw=1)
    axes[1].set(xlabel="time", ylabel="resistant fraction", ylim=(0, 1), title="Selection-relative composition")
    axes[2].axis("off")
    axes[2].text(0.05, 0.90, "Validation setting", fontsize=12, weight="bold")
    axes[2].text(0.05, 0.76, f"a0 = 2.15\ni = 10\nh_i = 2\nd_b = 1\n\nι = {iota:.3f}\nβ = {beta:.3f}\nφ = {phi:.3f}\nd_bι = {effective_db:.3f}\n\nEndpoint: f_R = {fr[-1]:.3f}", va="top")
    fig.tight_layout()
    savefig(fig, "expanded_02_validation_trajectory")


def _grid_pivot(frame: pd.DataFrame, value: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs = np.sort(frame["x_value"].unique())
    ys = np.sort(frame["y_value"].unique())
    matrix = frame.pivot(index="y_value", columns="x_value", values=value).reindex(index=ys, columns=xs).to_numpy()
    return xs, ys, matrix


def phase_panel(frame: pd.DataFrame, a0: float, plane: str) -> None:
    measures = [
        ("delta_fR", "Delta f_R", "coolwarm", None),
        ("final_fR", "final f_R", "viridis", (0, 1)),
        ("delta_log_ratio", "Delta log(n_r/n_s)", "coolwarm", None),
        ("auc_b", "free-Bla AUC", "magma", None),
        ("time_to_recovery", "recovery time", "plasma", None),
        ("peak_loss", "peak population loss", "cividis", (0, 1)),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.2))
    xs, ys, _ = _grid_pivot(frame, "delta_fR")
    for ax, (metric, title, cmap, limits) in zip(axes.flat, measures):
        _, _, matrix = _grid_pivot(frame, metric)
        values = matrix[np.isfinite(matrix)]
        if len(values) == 0:
            vmin, vmax = 0, 1
        elif limits is not None:
            vmin, vmax = limits
        elif metric in {"delta_fR", "delta_log_ratio"}:
            vmax = max(abs(float(np.nanmin(values))), abs(float(np.nanmax(values))))
            vmin = -vmax
        else:
            vmin, vmax = float(np.nanmin(values)), float(np.nanmax(values))
            if vmax == vmin:
                vmax = vmin + 1e-9
        mesh = ax.pcolormesh(xs, ys, matrix, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        if np.all(xs > 0) and np.max(xs) / np.min(xs) > 20:
            ax.set_xscale("log")
        ax.set(xlabel=frame["x_name"].iloc[0], ylabel=frame["y_name"].iloc[0], title=title)
        fig.colorbar(mesh, ax=ax, shrink=0.85)
    fig.suptitle(f"Response atlas: {plane}, a0={a0:g}", y=1.01)
    fig.tight_layout()
    savefig(fig, f"expanded_atlas_a0_{a0:g}_{plane}")


def summarize_anchors(all_atlas: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for a0, group in all_atlas.groupby("a0"):
        direction_counts = group["selection_direction"].value_counts(normalize=True)
        rows.append({
            "a0": a0,
            "n_runs": len(group),
            "min_delta_fR": group["delta_fR"].min(),
            "median_delta_fR": group["delta_fR"].median(),
            "max_delta_fR": group["delta_fR"].max(),
            "median_final_fR": group["final_fR"].median(),
            "median_auc_b": group["auc_b"].median(),
            "median_peak_loss": group["peak_loss"].median(),
            "fraction_resistant_enriched": direction_counts.get("resistant-enriched", 0.0),
            "fraction_susceptible_enriched": direction_counts.get("susceptible-enriched", 0.0),
            "fraction_neutral": direction_counts.get("neutral", 0.0),
            "fraction_recovered": group["susceptible_recovered"].mean(),
            "fraction_solver_success": group["solver_success"].mean(),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES / "anchor_summary.csv", index=False)
    return summary


def anchor_plot(summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    x = summary["a0"].to_numpy()
    axes[0].fill_between(x, summary["min_delta_fR"], summary["max_delta_fR"], alpha=0.2, color="tab:blue", label="range")
    axes[0].plot(x, summary["median_delta_fR"], "o-", color="tab:blue", label="median")
    axes[0].axhline(0, color="0.4", lw=1)
    axes[0].set(xlabel="initial antibiotic a0", ylabel="Delta f_R", title="Selection effect across atlas")
    axes[0].set_xscale("log")
    axes[0].legend(frameon=False)
    axes[1].plot(x, summary["median_auc_b"], "o-", color="tab:orange")
    axes[1].set(xlabel="initial antibiotic a0", ylabel="median free-Bla AUC", title="Public protection exposure")
    axes[1].set_xscale("log")
    axes[2].plot(x, summary["fraction_resistant_enriched"], "o-", label="R-enriched")
    axes[2].plot(x, summary["fraction_susceptible_enriched"], "o-", label="S-enriched")
    axes[2].plot(x, summary["fraction_neutral"], "o-", label="neutral")
    axes[2].set(xlabel="initial antibiotic a0", ylabel="fraction of atlas runs", ylim=(0, 1), title="Selection-direction composition")
    axes[2].set_xscale("log")
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    savefig(fig, "expanded_03_anchor_comparison")


def choose_representatives(all_atlas: pd.DataFrame) -> pd.DataFrame:
    selected = []
    for a0 in ANCHORS:
        group = all_atlas[all_atlas["a0"] == a0].copy()
        if group.empty:
            continue
        for label, index in [
            ("most_susceptible_enriched", group["delta_fR"].idxmin()),
            ("near_neutral", (group["delta_fR"].abs()).idxmin()),
            ("most_resistant_enriched", group["delta_fR"].idxmax()),
            ("deepest_bottleneck", group["peak_loss"].idxmax()),
        ]:
            row = group.loc[index].copy()
            row["representative"] = label
            selected.append(row)
    reps = pd.DataFrame(selected).drop_duplicates(subset=["a0", "representative"])
    reps.to_csv(TABLES / "representative_regimes_expanded.csv", index=False)
    return reps


def representative_trajectory_plot(reps: pd.DataFrame) -> None:
    fig, axes = plt.subplots(len(ANCHORS), 3, figsize=(15, 4.4 * len(ANCHORS)), squeeze=False)
    for row_idx, a0 in enumerate(ANCHORS):
        sub = reps[reps["a0"] == a0]
        for col_idx, label in enumerate(["most_susceptible_enriched", "near_neutral", "most_resistant_enriched"]):
            ax = axes[row_idx, col_idx]
            chosen = sub[sub["representative"] == label]
            if chosen.empty:
                ax.axis("off")
                continue
            rec = chosen.iloc[0]
            p = dict(BASE, i=float(rec["i"]), h_i=float(rec["h_i"]), d_b=float(rec["d_b"]))
            sim = simulate(p, a0=float(a0), t_end=T_ENDS[float(a0)])
            t = sim["t"]
            ns, nr, s, a, b = sim["y"]
            fr = nr / np.maximum(ns + nr, 1e-300)
            ax.plot(t, ns, color="tab:green", label="n_s")
            ax.plot(t, nr, color="tab:gray", label="n_r")
            ax.plot(t, fr, color="tab:purple", ls="--", label="f_R")
            ax.set_ylim(bottom=0)
            ax.set(xlabel="time", ylabel="state / fraction", title=f"a0={a0:g}: {label.replace('_', ' ')}\ni={rec['i']:.3g}, h_i={rec['h_i']:.2g}, d_b={rec['d_b']:.2g}\nDelta f_R={rec['delta_fR']:.3f}")
            if row_idx == 0 and col_idx == 0:
                ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Outcome-selected trajectories (observed selection direction)", y=1.005)
    fig.tight_layout()
    savefig(fig, "expanded_04_representative_trajectories")


def composition_analysis(reps: pd.DataFrame) -> pd.DataFrame:
    # Use the baseline anchor and three regimes, while retaining fixed total
    # biomass at each initial composition.
    anchor = 2.15
    target_labels = ["most_susceptible_enriched", "near_neutral", "most_resistant_enriched"]
    selected = reps[(reps["a0"] == anchor) & reps["representative"].isin(target_labels)].drop_duplicates("representative")
    fractions = np.linspace(0.05, 0.95, 19)
    rows = []
    for _, rec in selected.iterrows():
        p = dict(BASE, i=float(rec["i"]), h_i=float(rec["h_i"]), d_b=float(rec["d_b"]))
        for fr0 in fractions:
            total0 = NS0 + NR0
            nr0 = total0 * fr0
            ns0 = total0 * (1.0 - fr0)
            reference = no_treatment_reference(p, ns0=ns0, nr0=nr0, a0=0.0, t_end=T_ENDS[anchor])
            sim = simulate(p, ns0=ns0, nr0=nr0, a0=anchor, t_end=T_ENDS[anchor])
            metrics = endpoint_metrics(sim, p, reference=reference, ns0=ns0, nr0=nr0)
            rows.append(record_with_parameters(metrics, p, a0=anchor, regime=rec["representative"], initial_fR=fr0))
    frame = pd.DataFrame(rows)
    frame.to_csv(TABLES / "initial_composition_expanded.csv", index=False)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for regime, group in frame.groupby("regime"):
        axes[0, 0].plot(group["initial_fR"], group["final_fR"], "o-", label=regime.replace("_", " "))
        axes[0, 1].plot(group["initial_fR"], group["delta_fR"], "o-", label=regime.replace("_", " "))
        axes[1, 0].plot(group["initial_fR"], group["final_H"], "o-", label=regime.replace("_", " "))
        axes[1, 1].plot(group["initial_fR"], group["time_to_recovery"], "o-", label=regime.replace("_", " "))
    axes[0, 0].set(ylabel="final f_R", title="Initial composition shifts endpoint")
    axes[0, 1].axhline(0, color="0.4", lw=1)
    axes[0, 1].set(ylabel="Delta f_R", title="Direction of selection")
    axes[1, 0].set(xlabel="initial f_R", ylabel="final Shannon H", title="Final diversity")
    axes[1, 1].set(xlabel="initial f_R", ylabel="recovery time", title="Recovery timing")
    for ax in axes.flat:
        ax.set_xlabel(ax.get_xlabel() or "initial f_R")
        ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Fixed-total-biomass initial-composition analysis (a0=2.15)", y=1.01)
    fig.tight_layout()
    savefig(fig, "expanded_05_initial_composition")
    return frame


def ablation_analysis(reps: pd.DataFrame) -> pd.DataFrame:
    target_labels = ["most_susceptible_enriched", "near_neutral", "most_resistant_enriched"]
    selected = reps[(reps["a0"] == 2.15) & reps["representative"].isin(target_labels)].drop_duplicates("representative")
    modes = ["full", "private_protection_only", "public_degradation_only", "free_bla_turnover_only", "inhibitor_independent_control"]
    rows = []
    for _, rec in selected.iterrows():
        p = dict(BASE, i=float(rec["i"]), h_i=float(rec["h_i"]), d_b=float(rec["d_b"]))
        reference = no_treatment_reference(p, a0=0.0, t_end=T_ENDS[2.15])
        for mode in modes:
            sim = simulate(p, a0=2.15, t_end=T_ENDS[2.15], mode=mode)
            metrics = endpoint_metrics(sim, p, reference=reference)
            rows.append(record_with_parameters(metrics, p, a0=2.15, regime=rec["representative"], mode=mode))
    frame = pd.DataFrame(rows)
    frame.to_csv(TABLES / "mechanistic_ablations.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
    labels = frame["mode"].drop_duplicates().tolist()
    positions = np.arange(len(labels))
    width = 0.23
    for offset, regime in enumerate(frame["regime"].drop_duplicates().tolist()):
        group = frame[frame["regime"] == regime].set_index("mode").reindex(labels)
        axes[0].bar(positions + (offset - 1) * width, group["delta_fR"], width, label=regime.replace("_", " "))
        axes[1].bar(positions + (offset - 1) * width, group["auc_b"], width)
        axes[2].bar(positions + (offset - 1) * width, group["peak_loss"], width)
    axes[0].axhline(0, color="0.4", lw=1)
    axes[0].set_ylabel("Delta f_R")
    axes[1].set_ylabel("free-Bla AUC")
    axes[2].set_ylabel("peak population loss")
    for ax, title in zip(axes, ["Selection", "Public protection exposure", "Bottleneck"]):
        ax.set_title(title)
        ax.set_xticks(positions, [x.replace("_", "\n") for x in labels], rotation=0, fontsize=7)
    axes[0].legend(frameon=False, fontsize=7)
    fig.suptitle("Mechanistic ablations at selected regimes", y=1.02)
    fig.tight_layout()
    savefig(fig, "expanded_06_mechanistic_ablations")
    return frame


def ablation_heatmaps() -> pd.DataFrame:
    """Map each ablation over inhibitor dose and Hill steepness.

    Delta-delta-fR is calculated against the inhibitor-independent control at
    the same (i, h_i) coordinate. This distinguishes each run's within-run
    selection shift from the effect attributable to the mechanism being varied.
    """
    modes = [
        "full",
        "private_protection_only",
        "public_degradation_only",
        "free_bla_turnover_only",
        "inhibitor_independent_control",
    ]
    outputs = []
    baseline = dict(BASE)
    for mode in modes:
        if mode == "full":
            frame = pd.read_csv(TABLES / "atlas_a0_2.15_i_hi.csv")
            frame["pathway_mode"] = mode
        elif mode == "inhibitor_independent_control":
            # The control does not depend on i or h_i. One integration suffices;
            # replicate its endpoints over the plotting grid for comparison.
            p = dict(BASE)
            t_end = T_ENDS[2.15]
            reference = no_treatment_reference(p, a0=0.0, t_end=t_end)
            sim = simulate(p, a0=2.15, t_end=t_end, mode=mode)
            metrics = endpoint_metrics(sim, p, reference=reference)
            rows = []
            for hi in HI_VALUES:
                for dose in I_VALUES:
                    q = dict(BASE, i=float(dose), h_i=float(hi))
                    rows.append(record_with_parameters(
                        metrics, q, a0=2.15, plane="i_hi", pathway_mode=mode,
                        x_name="i", y_name="h_i", x_value=float(dose), y_value=float(hi),
                    ))
            frame = pd.DataFrame(rows)
        else:
            frame = run_pair_grid("i", I_VALUES, "h_i", HI_VALUES, baseline, 2.15, f"ablation_{mode}", mode=mode)
        frame["ablation_mode"] = mode
        outputs.append(frame)
        print(f"  completed ablation heatmap grid: {mode}", flush=True)
    result = pd.concat(outputs, ignore_index=True)
    # The control equations do not depend on either grid coordinate, so its
    # single Δf_R is the matched reference everywhere. Avoid floating-point
    # joins on serialized log-grid values.
    control_delta = float(result.loc[result["ablation_mode"] == "inhibitor_independent_control", "delta_fR"].iloc[0])
    result["control_delta_fR"] = control_delta
    result["delta_delta_fR_vs_control"] = result["delta_fR"] - result["control_delta_fR"]
    result.to_csv(TABLES / "ablation_heatmap_i_hi.csv", index=False)

    modes_and_titles = [
        ("full", "Full pathway"),
        ("private_protection_only", "Private protection only"),
        ("public_degradation_only", "Public degradation only"),
        ("free_bla_turnover_only", "Free-Bla turnover only"),
        ("inhibitor_independent_control", "Inhibitor-independent control"),
    ]
    for metric, suptitle, filename in [
        ("delta_fR", "Within-run selection change, Δf_R", "expanded_11_ablation_heatmaps_delta_fR"),
        ("delta_delta_fR_vs_control", "Ablation effect relative to control, ΔΔf_R", "expanded_12_ablation_heatmaps_vs_control"),
    ]:
        fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True, sharey=True)
        arrays = [result.loc[result["ablation_mode"] == mode, metric].to_numpy() for mode, _ in modes_and_titles]
        if metric == "delta_fR":
            bound = max(abs(np.nanmin(np.concatenate(arrays))), abs(np.nanmax(np.concatenate(arrays))))
        else:
            bound = max(abs(np.nanmin(np.concatenate(arrays))), abs(np.nanmax(np.concatenate(arrays))))
        bound = max(float(bound), 1e-8)
        for idx, (mode, title) in enumerate(modes_and_titles):
            ax = axes.flat[idx]
            subset = result[result["ablation_mode"] == mode]
            xs, ys, matrix = _grid_pivot(subset, metric)
            mesh = ax.pcolormesh(xs, ys, matrix, shading="auto", cmap="coolwarm", vmin=-bound, vmax=bound)
            ax.set_xscale("log")
            ax.set_title(title)
            ax.set_xlabel("inhibitor dose i")
            ax.set_ylabel("Hill coefficient h_i")
            fig.colorbar(mesh, ax=ax, shrink=0.85)
        axes.flat[-1].axis("off")
        fig.suptitle(f"Ablation response maps at a0=2.15\n{suptitle}", y=1.01)
        fig.tight_layout()
        savefig(fig, filename)
    return result


def robustness_analysis(seed: int = 574, n_samples: int = 48) -> pd.DataFrame:
    ranges = {
        "alpha": (0.75, 1.0),
        "beta_min": (0.0, 1.0),
        "xi": (0.0, 1.0),
        "kappa_b": (0.0, 1.0),
        "gamma": (1.1, 1.4),
        "h_a": (1.0, 5.0),
        "phi_max": (0.0, 5.0),
        "c": (0.0, 1.0),
    }
    sampler = qmc.LatinHypercube(d=len(ranges), seed=seed)
    unit = sampler.random(n=n_samples)
    names = list(ranges)
    sampled = qmc.scale(unit, [ranges[name][0] for name in names], [ranges[name][1] for name in names])
    rows = []
    for idx, values in enumerate(sampled):
        p = dict(BASE)
        p.update({name: float(value) for name, value in zip(names, values)})
        ref = no_treatment_reference(p, a0=0.0, t_end=T_ENDS[2.15])
        for dose in [0.1, 10.0]:
            p_dose = dict(p, i=dose, h_i=2.0, d_b=1.0)
            sim = simulate(p_dose, a0=2.15, t_end=T_ENDS[2.15])
            metrics = endpoint_metrics(sim, p_dose, reference=ref)
            rows.append(record_with_parameters(metrics, p_dose, sample=idx, dose=dose))
        if (idx + 1) % 8 == 0 or idx + 1 == n_samples:
            print(f"  robustness: {idx + 1}/{n_samples}", flush=True)
    frame = pd.DataFrame(rows)
    wide = frame.pivot(index="sample", columns="dose", values="delta_fR")
    wide.columns = [f"delta_fR_dose_{str(c).replace('.', 'p')}" for c in wide.columns]
    frame = frame.merge(wide, left_on="sample", right_index=True, how="left")
    high_col = "delta_fR_dose_10p0"
    low_col = "delta_fR_dose_0p1"
    frame["high_minus_low_delta_fR"] = frame[high_col] - frame[low_col]
    frame.to_csv(TABLES / "robustness_lhs.csv", index=False)
    high = frame[frame["dose"] == 10.0].copy()
    corr_rows = []
    for name in names:
        rho, pval = spearmanr(high[name], high["high_minus_low_delta_fR"], nan_policy="omit")
        corr_rows.append({"parameter": name, "spearman_rho": rho, "p_value": pval})
    pd.DataFrame(corr_rows).to_csv(TABLES / "robustness_rank_correlations.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    high_effect = high["high_minus_low_delta_fR"].dropna()
    axes[0].hist(high_effect, bins=12, color="tab:blue", alpha=0.8)
    axes[0].axvline(0, color="0.4", lw=1)
    axes[0].set(xlabel="high-dose minus low-dose Delta f_R", ylabel="sample count", title="Robustness effect distribution")
    strongest = pd.DataFrame(corr_rows).sort_values("spearman_rho", key=np.abs, ascending=False).head(2)["parameter"].tolist()
    for name in strongest:
        axes[1].scatter(high[name], high["high_minus_low_delta_fR"], s=18, alpha=0.75, label=name)
    axes[1].axhline(0, color="0.4", lw=1)
    axes[1].set(xlabel="sampled parameter", ylabel="high-minus-low effect", title="Strongest rank associations")
    axes[1].legend(frameon=False, fontsize=8)
    directions = high["selection_direction"].value_counts(normalize=True)
    axes[2].bar(directions.index, directions.values, color=["tab:orange", "tab:green", "tab:gray"][:len(directions)])
    axes[2].set(ylim=(0, 1), ylabel="fraction", title="High-dose endpoint directions")
    axes[2].tick_params(axis="x", rotation=25)
    fig.tight_layout()
    savefig(fig, "expanded_07_robustness")
    return frame


def identifiability_analysis() -> tuple[pd.DataFrame, pd.DataFrame]:
    target_iotas = [0.1, 0.5, 0.9]
    hi_choices = [1.0, 2.0, 5.0]
    rows = []
    for target in target_iotas:
        for hi in hi_choices:
            i_value = (target / (1.0 - target)) ** (1.0 / hi)
            p = dict(BASE, i=i_value, h_i=hi, d_b=1.0)
            ref = no_treatment_reference(p, a0=0.0, t_end=T_ENDS[2.15])
            sim = simulate(p, a0=2.15, t_end=T_ENDS[2.15])
            metrics = endpoint_metrics(sim, p, reference=ref)
            rows.append(record_with_parameters(metrics, p, target_iota=target, matched_group=f"iota_{target:g}"))
    frame = pd.DataFrame(rows)
    frame.to_csv(TABLES / "identifiability_matched_iota.csv", index=False)
    db_rows = []
    for db in DB_VALUES:
        p = dict(BASE, i=1.0, h_i=2.0, d_b=float(db))
        ref = no_treatment_reference(p, a0=0.0, t_end=T_ENDS[2.15])
        sim = simulate(p, a0=2.15, t_end=T_ENDS[2.15])
        metrics = endpoint_metrics(sim, p, reference=ref)
        db_rows.append(record_with_parameters(metrics, p, matched_iota=0.5))
    db_frame = pd.DataFrame(db_rows)
    db_frame.to_csv(TABLES / "identifiability_db_at_fixed_iota.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for target, group in frame.groupby("target_iota"):
        axes[0].plot(group["h_i"], group["i"], "o-", label=f"iota={target:g}")
        axes[1].plot(group["h_i"], group["final_fR"], "o-", label=f"iota={target:g}")
        axes[2].plot(group["h_i"], group["auc_b"], "o-", label=f"iota={target:g}")
    axes[0].set_yscale("log")
    axes[0].set(xlabel="h_i", ylabel="matched i", title="Many (i, h_i) pairs share iota")
    axes[1].set(xlabel="h_i", ylabel="final f_R", title="Endpoint under matched iota")
    axes[2].set(xlabel="h_i", ylabel="free-Bla AUC", title="Trajectory exposure under matched iota")
    for ax in axes:
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    savefig(fig, "expanded_08_matched_iota_identifiability")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    axes[0].plot(db_frame["d_b"], db_frame["final_fR"], "o-")
    axes[1].plot(db_frame["d_b"], db_frame["auc_b"], "o-", color="tab:orange")
    axes[2].plot(db_frame["d_b"], db_frame["time_to_recovery"], "o-", color="tab:green")
    axes[0].set(xlabel="d_b", ylabel="final f_R", title="Endpoint composition")
    axes[1].set(xlabel="d_b", ylabel="free-Bla AUC", title="Free-Bla exposure")
    axes[2].set(xlabel="d_b", ylabel="recovery time", title="Recovery timing")
    fig.suptitle("d_b remains informative at fixed inhibitor activity", y=1.02)
    fig.tight_layout()
    savefig(fig, "expanded_09_db_identifiability")
    return frame, db_frame


def selection_criterion_plot() -> None:
    doses = np.logspace(-1, 1, 100)
    rows = []
    for hi in [1.0, 2.0, 5.0]:
        for i in doses:
            p = dict(BASE, i=float(i), h_i=hi)
            iota, beta, _, _ = pathway_values(p)
            lhs = 1.0 - beta
            rhs = (1.0 - p["alpha"]) / p["gamma"]
            rows.append({"i": i, "h_i": hi, "lhs": lhs, "rhs": rhs, "iota": iota})
    frame = pd.DataFrame(rows)
    frame.to_csv(TABLES / "selection_criterion_scan.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for hi, group in frame.groupby("h_i"):
        axes[0].plot(group["i"], group["lhs"], label=f"1-beta, h_i={hi:g}")
    axes[0].axhline(frame["rhs"].iloc[0], color="black", ls="--", label="(1-alpha)/gamma")
    axes[0].set_xscale("log")
    axes[0].set(xlabel="i", ylabel="selection-criterion terms", title="Saturating selection criterion")
    axes[0].legend(frameon=False, fontsize=7)
    for hi, group in frame.groupby("h_i"):
        axes[1].plot(group["i"], group["lhs"] - group["rhs"], label=f"h_i={hi:g}")
    axes[1].axhline(0, color="black", ls="--")
    axes[1].set_xscale("log")
    axes[1].set(xlabel="i", ylabel="criterion margin", title="Positive margin favors resistant enrichment")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    savefig(fig, "expanded_10_selection_criterion")


def write_config(all_atlas: pd.DataFrame) -> None:
    config = {
        "model": "Ma et al. 2024 five-state dimensionless two-population model",
        "state_order": ["n_s", "n_r", "s", "a", "b"],
        "base_parameters": BASE,
        "initial_conditions": {"n_s": NS0, "n_r": NR0, "s": S0, "a": A0, "b": B0},
        "antibiotic_anchors": ANCHORS,
        "time_horizons": T_ENDS,
        "grid": {"i": I_VALUES.tolist(), "h_i": HI_VALUES.tolist(), "d_b": DB_VALUES.tolist(), "h_i_d_b_doses": HI_DB_DOSES},
        "solver": {"method": "LSODA", "rtol": 3e-7, "atol": 1e-9, "max_step": 2.0, "n_points": N_POINTS},
        "thresholds": {"extinction": EXTINCTION_THRESHOLD, "recovery_fraction": RECOVERY_FRACTION, "selection_tolerance": SELECTION_TOLERANCE},
        "n_atlas_rows": int(len(all_atlas)),
        "local_only": False,
        "notes": "Expanded Person 3 analysis; interpretations remain conditional on group review of baselines and endpoints.",
    }
    (TABLES / "analysis_config.json").write_text(json.dumps(config, indent=2))


def para(text: str, style) -> Paragraph:
    return Paragraph(text.replace("&", "&amp;"), style)


def report_table(data: list[list[str]], widths: list[float] | None = None) -> Table:
    table = Table(data, colWidths=widths, repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#d9e8f5")),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.grey),
        ("FONTNAME", (0, 0), (-1, 0), "DejaVuSans-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7.5),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f8fb")]),
    ]))
    return table


def report_image(path: Path, max_width: float = 7.1 * inch, max_height: float = 5.25 * inch) -> Image:
    """Embed a plot without distorting its native aspect ratio."""
    reader = ImageReader(str(path))
    pixel_width, pixel_height = reader.getSize()
    scale = min(max_width / pixel_width, max_height / pixel_height)
    return Image(str(path), width=pixel_width * scale, height=pixel_height * scale)


def wrapped_table(data: list[list[str]], styles, widths: list[float] | None = None) -> Table:
    """Wrap long Unicode table cells so parameter meanings remain readable."""
    wrapped = []
    for row in data:
        wrapped.append([para(str(cell), styles["TableCell"]) for cell in row])
    return report_table(wrapped, widths=widths)


def build_report(all_atlas: pd.DataFrame, summary: pd.DataFrame, reps: pd.DataFrame, composition: pd.DataFrame, ablations: pd.DataFrame, robustness: pd.DataFrame, matched: pd.DataFrame, db_frame: pd.DataFrame, ablation_maps: pd.DataFrame | None = None, comprehensive: dict | None = None) -> None:
    report_path = REPORTS / "Person_3_Expanded_Report.pdf"
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small", parent=styles["BodyText"], fontSize=8.4, leading=10.5, spaceAfter=5))
    styles.add(ParagraphStyle(name="Tiny", parent=styles["BodyText"], fontSize=7.2, leading=8.7, spaceAfter=3))
    styles.add(ParagraphStyle(name="UnicodeSmall", parent=styles["BodyText"], fontName="DejaVuSans", fontSize=8.4, leading=11.2, spaceAfter=5))
    styles.add(ParagraphStyle(name="Equation", parent=styles["BodyText"], fontName="DejaVuSans", fontSize=10.5, leading=15, leftIndent=14, spaceAfter=3))
    styles.add(ParagraphStyle(name="TableCell", parent=styles["BodyText"], fontName="DejaVuSans", fontSize=7.1, leading=8.6, spaceAfter=0))
    styles["Title"].fontSize = 20
    doc = SimpleDocTemplate(str(report_path), pagesize=letter, rightMargin=0.55 * inch, leftMargin=0.55 * inch, topMargin=0.55 * inch, bottomMargin=0.55 * inch)
    story = []
    story.append(Paragraph("Person 3 Expanded Analysis Report", styles["Title"]))
    story.append(Paragraph("Ma two-population beta-lactam/Bla inhibitor model; Person 3 analysis", styles["Heading2"]))
    story.append(Paragraph("This report presents the exhaustive Person 3 analysis. The primer below defines the model before presenting the analysis, so the plots can be interpreted mechanistically rather than as purely mathematical parameter sweeps.", styles["Small"]))
    story.append(Paragraph("Model equations and mechanistic interpretation", styles["Heading1"]))
    story.append(Paragraph("The model is dimensionless. τ is scaled time; nₛ and nᵣ are susceptible and resistant cell abundances; s is the shared nutrient/resource pool; a is active antibiotic; and b is extracellular free beta-lactamase (Bla). The equations below are the Ma paper/supplement base model used for this local analysis.", styles["UnicodeSmall"]))
    equations = [
        "dnₛ/dτ = (g − ℓ)nₛ",
        "dnᵣ/dτ = (αg − βℓ)nᵣ",
        "ds/dτ = (ξℓ − g)nₛ + (ξβℓ − αg)nᵣ",
        "da/dτ = −κᵦba − φnᵣa − dₐa",
        "db/dτ = βℓnᵣ − dᵦιb",
        "g = s/(1+s)",
        "ℓ = γ [aʰᵃ/(1+aʰᵃ)] g",
        "ι = iʰⁱ/(1+iʰⁱ)",
        "β = βₘᵢₙ + c(1−βₘᵢₙ)ι",
        "φ = φₘₐₓ(1−cι)",
    ]
    for equation in equations:
        story.append(Paragraph(equation, styles["Equation"]))
    story.append(Paragraph("Mechanistic reading: g is nutrient-limited bacterial growth; ℓ is antibiotic-mediated lysis; ι is the fractional inhibitor activity; β is the resistant-cell lysis factor after inhibitor action; and φ is the rate at which living resistant cells degrade antibiotic. Free Bla is produced when resistant cells lyse and is removed at the inhibitor-dependent rate dᵦι. Thus, inhibitor activity simultaneously weakens private resistance protection through β and public antibiotic degradation through φ.", styles["UnicodeSmall"]))
    state_rows = [
        ["State/control", "Where it appears", "Mechanistic real-world equivalent"],
        ["nₛ", "First ODE, nutrient and lysis terms", "Number or density of antibiotic-susceptible bacteria; the population that experiences the full lysis rate ℓ."],
        ["nᵣ", "Second, third, fourth, and fifth ODEs", "Number or density of resistant bacteria; the population carrying Bla/resistance and paying the growth burden α."],
        ["s", "g and the nutrient-balance equation", "Shared nutrient or limiting resource available to both populations; lysis can recycle part of it through ξ."],
        ["a", "ℓ and antibiotic-loss equation", "Active antibiotic concentration or effective antibiotic exposure in the environment."],
        ["b", "Antibiotic-loss and Bla-balance equations", "Free extracellular Bla enzyme; a public-good pool that degrades antibiotic but is inactivated by the inhibitor."],
        ["τ", "All differential equations", "Dimensionless time; map to experimental time only after specifying the model time scale."],
    ]
    story.append(wrapped_table(state_rows, styles, widths=[0.85 * inch, 1.75 * inch, 4.5 * inch]))
    story.append(Spacer(1, 0.12 * inch))
    parameter_rows = [
        ["Parameter", "Equation role", "Mechanistic real-world equivalent", "Effect of increasing it"],
        ["α", "αg in resistant growth; αg in nutrient use", "Fitness of resistant bacteria relative to susceptible bacteria when antibiotic pressure is absent; α<1 represents a resistance cost.", "Raises resistant growth and nutrient consumption; usually makes resistant enrichment more likely."],
        ["βₘᵢₙ", "βℓ without inhibitor", "Baseline resistant-cell lysis sensitivity. 1−βₘᵢₙ is the private protection supplied by resistance/Bla inside or near the cell.", "Raises resistant lysis at low inhibitor activity; reduces the private benefit of resistance."],
        ["ξ", "ξℓ and ξβℓ in ds/dτ", "Fraction of lysed biomass converted back into shared usable nutrient; a proxy for recycling efficiency.", "Increases resource regeneration after lysis and can support later regrowth."],
        ["dₐ", "−dₐa", "Background antibiotic clearance or decay independent of Bla and resistant cells; e.g., chemical degradation, dilution, or washout.", "Shortens antibiotic persistence and generally reduces treatment pressure duration."],
        ["κᵦ", "−κᵦba", "Catalytic efficiency of free Bla in degrading active antibiotic; the strength of extracellular public protection.", "Makes a given free-Bla concentration remove antibiotic faster."],
        ["γ", "ℓ = γ·Hill(a)·g", "Maximum antibiotic-induced lysis strength relative to nutrient-limited growth; an antibiotic potency/lysis-scale parameter.", "Increases lysis at a given antibiotic exposure and usually deepens the initial population bottleneck."],
        ["hₐ", "Hill exponent in ℓ", "Steepness of the antibiotic concentration-response curve; how switch-like lysis is around the antibiotic transition.", "Makes the antibiotic response more threshold-like, with effects depending on where a lies relative to the transition."],
        ["i", "Hill exponent input for ι", "Dimensionless inhibitor dose, scaled by the inhibitor half-maximal effective concentration (IC₅₀-like reference).", "Increases inhibitor activity and therefore strengthens inhibitor effects when hᵢ is fixed."],
        ["hᵢ", "Hill exponent in ι", "Steepness/cooperativity of inhibitor action; how sharply inhibitor activity changes around i=1.", "Increases ι for i>1 but decreases ι for i<1, so it is not globally equivalent to a stronger dose."],
        ["dᵦ", "−dᵦιb", "Scale of inhibitor-dependent free-Bla inactivation/turnover; a proxy for enzyme neutralization or loss rate.", "Shortens free-Bla persistence and reduces the duration of extracellular antibiotic protection."],
        ["φₘₐₓ", "φnᵣa", "Maximum antibiotic degradation rate by living resistant cells; the cell-associated/public degradation pathway.", "Makes resistant cells remove antibiotic more effectively when inhibitor activity is low."],
        ["c", "β and φ inhibitor coupling", "Fractional efficacy of inhibitor action on the resistance/protection pathway; it couples intracellular/private and public effects in this model.", "Makes inhibitor activity more strongly increase β and decrease φ."],
        ["ι", "β, φ, and dᵦι", "Effective fraction of inhibitor target activity, from zero to one; the direct pharmacodynamic variable generated by i and hᵢ.", "Increasing ι weakens private protection, weakens public degradation, and increases free-Bla turnover."],
    ]
    story.append(wrapped_table(parameter_rows, styles, widths=[0.65 * inch, 1.25 * inch, 3.25 * inch, 1.95 * inch]))
    story.append(Spacer(1, 0.12 * inch))
    story.append(Paragraph("Initial conditions are nₛ(0)=nᵣ(0)=0.2, s(0)=4, b(0)=0, and a(0)=a₀. The analysis uses a₀=2.15, 10, and 100 as antibiotic anchors. The inhibitor dose i and its response steepness hᵢ are Person 3’s main controls; the other parameters are frozen at the shared baseline unless a coordination/robustness run is explicitly identified.", styles["UnicodeSmall"]))
    story.append(Paragraph("Executive summary", styles["Heading1"]))
    anchor_text = []
    for _, row in summary.iterrows():
        anchor_text.append(f"a0={row['a0']:g}: median Delta f_R={row['median_delta_fR']:.4f}, range [{row['min_delta_fR']:.4f}, {row['max_delta_fR']:.4f}], resistant-enriched fraction={row['fraction_resistant_enriched']:.2f}, susceptible-enriched fraction={row['fraction_susceptible_enriched']:.2f}, recovered fraction={row['fraction_recovered']:.2f}")
    story.append(Paragraph("The expanded atlas contains " + f"{len(all_atlas):,}" + " simulations across three antibiotic anchors and three inhibitor parameter planes. " + " ".join(anchor_text), styles["Small"]))
    story.append(Paragraph("The primary interpretation is selection-relative: the change from the initial resistant fraction is reported alongside final composition. The inhibitor activity ι changes beta and phi simultaneously in the full model, and d_b controls the duration of free-Bla protection through the effective rate d_bι. Therefore, matched inhibitor activity does not imply matched trajectories when d_b differs.", styles["Small"]))
    story.append(Paragraph("Methods and computational contract", styles["Heading1"]))
    story.append(Paragraph("The equations are integrated with LSODA using a log transform for n_s and n_r to prevent high-dose underflow from corrupting the composition signal. The environmental states remain in their original coordinates and are clipped to nonnegative values for reporting. The baseline is alpha=0.95, beta_min=0.90, xi=0.80, d_a=0.02, kappa_b=0.35, gamma=1.38, h_a=3, h_i=2, d_b=1, phi_max=1, and c=0.70. Initial conditions are n_s=n_r=0.2, s=4, b=0, with antibiotic anchors a0=2.15, 10, and 100. The full configuration is in results/tables/analysis_config.json.", styles["Small"]))
    story.append(Paragraph("The atlas uses i from 0.1 to 10 on a logarithmic grid, h_i from 1 to 5, and d_b from 1 to 10. A run is classified as neutral when |Delta f_R| <= 0.001. Recovery means n_s reaches 90 percent of the no-treatment reference endpoint and remains there in the final trajectory segment. Peak loss is one minus the minimum total live biomass divided by the initial total biomass.", styles["Small"]))
    story.append(Paragraph("Atlas results", styles["Heading1"]))
    table_data = [["a0", "runs", "Delta f_R median", "Delta f_R range", "R-enriched", "S-enriched", "recovered"]]
    for _, row in summary.iterrows():
        table_data.append([
            f"{row['a0']:g}", f"{int(row['n_runs'])}", f"{row['median_delta_fR']:.4f}", f"[{row['min_delta_fR']:.3f}, {row['max_delta_fR']:.3f}]",
            f"{row['fraction_resistant_enriched']:.2f}", f"{row['fraction_susceptible_enriched']:.2f}", f"{row['fraction_recovered']:.2f}",
        ])
    story.append(report_table(table_data, widths=[0.5 * inch, 0.55 * inch, 0.85 * inch, 1.1 * inch, 0.75 * inch, 0.75 * inch, 0.65 * inch]))
    story.append(Spacer(1, 0.15 * inch))
    story.append(Paragraph("The corresponding response-atlas figures show whether the primary selection effect is accompanied by a biomass bottleneck, a change in free-Bla exposure, or loss of susceptible recovery. Companion maps should be read jointly: final f_R alone can be unfavorable even when total biomass is recovered, and a high f_R can also arise because the initial f_R was high.", styles["Small"]))
    for name, caption in [
        ("expanded_01_inhibitor_pathway", "Figure 1. Effective inhibitor pathway: iota, beta, phi, and d_b iota."),
        ("expanded_02_validation_trajectory", "Figure 2. Representative Ma-style trajectory and resistant-fraction change."),
        ("expanded_03_anchor_comparison", "Figure 3. Cross-anchor comparison of selection and exposure."),
        ("expanded_04_representative_trajectories", "Figure 4. Algorithmically selected trajectories for susceptible-enriched, near-neutral, and resistant-enriched outcomes."),
    ]:
        path = FIGURES / f"{name}.png"
        if path.exists():
            story.append(report_image(path))
            story.append(Paragraph(caption, styles["Tiny"]))
    story.append(PageBreak())
    story.append(Paragraph("Initial composition", styles["Heading1"]))
    story.append(Paragraph("At fixed initial total biomass, changing f_R(0) generally shifts the final composition as expected from the initial condition, but the slope and recovery behavior depend on the inhibitor regime. The selection-relative curves are the relevant test: if Delta f_R keeps the same sign over the initial-composition grid, the inhibitor setting has a robust direction of selection; if it crosses zero, the initial mixture is part of the endpoint mechanism.", styles["Small"]))
    story.append(report_image(FIGURES / "expanded_05_initial_composition.png"))
    story.append(Paragraph("Figure 5. Initial resistant fraction versus endpoint composition, selection change, diversity, and recovery.", styles["Tiny"]))
    story.append(Paragraph("Mechanistic ablations", styles["Heading1"]))
    story.append(Paragraph("The ablations isolate the three routes by which inhibitor activity enters the model. Private-protection-only changes beta; public-degradation-only changes phi; free-Bla-turnover-only changes d_b iota. Comparing these with the full pathway shows which mechanism explains the endpoint and which primarily changes the transient antibiotic/Bla exposure.", styles["Small"]))
    story.append(report_image(FIGURES / "expanded_06_mechanistic_ablations.png"))
    story.append(Paragraph("Figure 6. Mechanistic ablations across the selected low, neutral, and high selection regimes.", styles["Tiny"]))
    ab_table = [["regime", "mode", "Delta f_R", "free-Bla AUC", "peak loss"]]
    for _, row in ablations.iterrows():
        ab_table.append([str(row["regime"]).replace("_", " "), str(row["mode"]).replace("_", " "), f"{row['delta_fR']:.4f}", f"{row['auc_b']:.3g}", f"{row['peak_loss']:.3f}"])
    story.append(report_table(ab_table, widths=[1.35 * inch, 1.35 * inch, 0.8 * inch, 0.85 * inch, 0.7 * inch]))
    if ablation_maps is not None:
        story.append(PageBreak())
        story.append(Paragraph("Ablation heatmaps across inhibitor settings", styles["Heading1"]))
        story.append(Paragraph("Each heatmap varies inhibitor dose i and Hill steepness h_i at a0=2.15, with d_b fixed at its baseline. Figure 7 shows Δf_R within each simulation, measured from its own initial composition (f_R(0)=0.5). Figure 8 shows ΔΔf_R: each mode’s Δf_R minus the inhibitor-independent control’s Δf_R at the same grid point. Thus the second set directly compares the ablation effect against the control. The full pathway and control are included as reference maps.", styles["UnicodeSmall"]))
        story.append(report_image(FIGURES / "expanded_11_ablation_heatmaps_delta_fR.png", max_height=6.6 * inch))
        story.append(Paragraph("Figure 7. Selection change within each ablation setting.", styles["Tiny"]))
        story.append(PageBreak())
        story.append(report_image(FIGURES / "expanded_12_ablation_heatmaps_vs_control.png", max_height=6.6 * inch))
        story.append(Paragraph("Figure 8. Difference in selection change relative to the matched inhibitor-independent control. Positive values mean more resistant enrichment (or less susceptible enrichment) than the control.", styles["Tiny"]))
        story.append(PageBreak())
    story.append(PageBreak())
    story.append(Paragraph("Robustness and identifiability", styles["Heading1"]))
    story.append(Paragraph("The robustness run samples the non-Person 3 parameters over the supplementary ranges. It is a coordination analysis: it tests whether the Person 3 interpretation survives plausible uncertainty, but it does not replace the group’s agreed baseline. Rank correlations and the full sampled table are saved under results/tables.", styles["Small"]))
    story.append(report_image(FIGURES / "expanded_07_robustness.png"))
    story.append(Paragraph("Figure 9. Latin-hypercube robustness analysis.", styles["Tiny"]))
    story.append(Paragraph("Matched-ι analysis", styles["Heading2"]))
    story.append(Paragraph("Because iota = i^h_i/(1+i^h_i), multiple pairs of i and h_i produce the same instantaneous inhibitor activity. Endpoint similarity within a matched group supports treating iota as the direct pathway variable, while remaining differences in free-Bla AUC, antibiotic exposure, or recovery indicate information carried by other rates or by numerical/transient effects. Varying d_b at fixed iota provides a direct test of the free-Bla time scale.", styles["Small"]))
    story.append(report_image(FIGURES / "expanded_08_matched_iota_identifiability.png"))
    story.append(report_image(FIGURES / "expanded_09_db_identifiability.png"))
    story.append(Paragraph("Figures 10-11. Matched inhibitor activity and fixed-activity d_b identifiability checks.", styles["Tiny"]))
    story.append(Paragraph("Interpretation, implications, and limitations", styles["Heading1"]))
    story.append(Paragraph("1. Inhibitor activity has two linked ecological effects in the full model: beta moves toward one, weakening private lysis protection, and phi decreases, weakening resistant-mediated antibiotic degradation. A dose setting can therefore change both selection and the depth/duration of the population bottleneck.\n2. h_i is not a monotonic 'strength' parameter over the entire dose range. Below i=1 it lowers iota relative to a lower h_i, while above i=1 it raises iota; the response atlas and pathway plot should be interpreted with this sign change in mind.\n3. d_b is a time-scale parameter. It can leave the instantaneous inhibitor activity unchanged while changing free-Bla persistence, antibiotic exposure, recovery timing, and possibly the final composition.\n4. Initial composition is not a nuisance detail. The same dynamic selection effect can produce different final f_R values depending on f_R(0), so Delta f_R and Delta log(n_r/n_s) should be reported with any endpoint claim.\n5. All conclusions remain conditional on the dimensionless base model, no spatial structure, no horizontal gene transfer, constant treatment parameters, and the group’s choices for antibiotic anchors and endpoint thresholds.", styles["Small"]))
    story.append(Paragraph("Hand-off to the other roles", styles["Heading1"]))
    story.append(Paragraph("Person 1 should confirm alpha, beta_min, c, the selection criterion, and the final robustness ranges. Person 2 should confirm the a0, gamma, and h_a anchors and whether the high-dose horizon is appropriate. Person 4 should confirm kappa_b, d_a, phi_max, xi uncertainty, direct a/b observables, and sampling timing. The group should then agree on the state order, solver tolerances, extinction/recovery thresholds, and which coordination runs belong in the final shared report. Until those decisions are made, the tables and figures in this report should be treated as a reproducible local analysis package rather than final group claims.", styles["Small"]))
    if comprehensive is not None:
        story.append(PageBreak())
        story.append(Paragraph("Comprehensive inhibitor analysis: full parameter sweep", styles["Heading1"]))
        story.append(Paragraph(f"To match the breadth of the other workstreams, this addendum evaluates a structured 21 × 21 × 21 grid ({comprehensive['n_runs']:,} simulations) over inhibitor dose i, inhibitor Hill coefficient hᵢ, and free-Bla inactivation rate dᵦ, at the shared Ma-style antibiotic anchor a₀=2.15. It reports the same family of endpoint maps, marginal effects, sensitivity models, tradeoffs, trajectory contrasts, and outcome regimes used in the other notebooks. The previously generated multi-anchor response atlas, composition analysis, robustness runs, identifiability study, ablations, and heatmaps remain included above.", styles["UnicodeSmall"]))
        csum = comprehensive["summary"]
        summary_rows = [["Endpoint", "Minimum", "Median", "Maximum"]]
        for _, row in csum.iterrows():
            summary_rows.append([str(row["endpoint"]), f"{row['min']:.4g}", f"{row['median']:.4g}", f"{row['max']:.4g}"])
        story.append(wrapped_table(summary_rows, styles, widths=[2.7 * inch, 1.1 * inch, 1.1 * inch, 1.1 * inch]))
        story.append(Spacer(1, 0.12 * inch))
        story.append(Paragraph(comprehensive["regime_text"], styles["UnicodeSmall"]))
        story.append(Paragraph("Validation and selection criterion", styles["Heading2"]))
        story.append(Paragraph(comprehensive["validation_text"], styles["UnicodeSmall"]))
        story.append(Paragraph("Global sensitivity and marginal response", styles["Heading2"]))
        story.append(Paragraph("The standardized regression includes standardized main effects and all pairwise interactions among i, hᵢ, and dᵦ. Coefficients are conditional associations over this structured grid, not causal effects outside the stated model. Binned curves show the mean and ±1 standard deviation across the other two varied parameters.", styles["UnicodeSmall"]))
        if "sensitivity_table" in comprehensive:
            story.append(wrapped_table(comprehensive["sensitivity_table"], styles, widths=[1.5 * inch, 1.65 * inch, 1.65 * inch, 1.0 * inch, 0.8 * inch]))
        for fig_name, caption in comprehensive["figures"]:
            story.append(PageBreak())
            story.append(Paragraph(caption, styles["Heading2"]))
            fig_path = FIGURES / f"{fig_name}.png"
            if fig_path.exists():
                story.append(report_image(fig_path, max_height=6.6 * inch))
        story.append(PageBreak())
        story.append(Paragraph("Sensitivity tables and reproducibility", styles["Heading2"]))
        story.append(Paragraph("Machine-readable full-grid results, endpoint summaries, binned response summaries, Spearman associations, standardized regression coefficients, criterion comparisons, convergence classifications, and trajectory-case selections are saved in results/tables/ with the `comprehensive_` prefix. The runner is person3_comprehensive_analysis.py. All comprehensive figures have both PNG and PDF versions under results/figures/ and are mirrored to the top-level figures/ folder.", styles["UnicodeSmall"]))
    doc.build(story)
    print(f"Saved report: {report_path}")


def main() -> None:
    print("Running expanded Person 3 analysis in local biobase-compatible script", flush=True)
    pathway_plot()
    validation_plot()
    all_frames = []
    for a0 in ANCHORS:
        fixed = dict(BASE)
        frame_i_hi = run_pair_grid("i", I_VALUES, "h_i", HI_VALUES, fixed, a0, "i_hi")
        frame_i_db = run_pair_grid("i", I_VALUES, "d_b", DB_VALUES, fixed, a0, "i_db")
        all_frames.extend([frame_i_hi, frame_i_db])
        for dose in HI_DB_DOSES:
            fixed_dose = dict(BASE, i=dose)
            frame_hi_db = run_pair_grid("h_i", HI_VALUES, "d_b", DB_VALUES, fixed_dose, a0, f"hi_db_i_{dose:g}")
            all_frames.append(frame_hi_db)
    all_atlas = pd.concat(all_frames, ignore_index=True)
    all_atlas.to_csv(TABLES / "all_atlas_results.csv", index=False)
    for a0 in ANCHORS:
        anchor = all_atlas[all_atlas["a0"] == a0]
        for plane in anchor["plane"].drop_duplicates():
            phase_panel(anchor[anchor["plane"] == plane], a0, plane)
    summary = summarize_anchors(all_atlas)
    anchor_plot(summary)
    reps = choose_representatives(all_atlas)
    representative_trajectory_plot(reps)
    composition = composition_analysis(reps)
    ablations = ablation_analysis(reps)
    ablation_maps = ablation_heatmaps()
    robustness = robustness_analysis()
    matched, db_frame = identifiability_analysis()
    selection_criterion_plot()
    write_config(all_atlas)
    build_report(all_atlas, summary, reps, composition, ablations, robustness, matched, db_frame, ablation_maps)
    print("Expanded analysis complete", flush=True)


if __name__ == "__main__":
    main()
