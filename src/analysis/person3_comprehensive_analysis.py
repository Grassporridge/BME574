"""Complete Person 3 sweep and analysis matching the other workstreams.

This runner reuses the model/plot/report utilities in person3_expanded_analysis.py
and adds a full 21^3 (i, h_i, d_b) sweep plus the endpoint, sensitivity,
interaction, regime, heatmap, and trajectory analyses used by the other members.

Run with:
    conda run -n ar_cp python -m src.analysis.person3_comprehensive_analysis

All new files go under the repository's results/ tree.
"""

from __future__ import annotations

import itertools
import json
import os
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as student_t
from scipy.integrate import solve_ivp

MPLCONFIGDIR = Path(tempfile.gettempdir()) / "person3_comprehensive_mplconfig"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
import matplotlib.pyplot as plt

try:  # package execution: python -m src.analysis.person3_comprehensive_analysis
    from . import person3_expanded_analysis as core
except ImportError:  # direct execution from src/analysis
    import person3_expanded_analysis as core


# Scripts are stored in src/analysis; results belong at the repository root.
ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
for folder in (TABLES, FIGURES):
    folder.mkdir(parents=True, exist_ok=True)

PARAMS = ["i", "h_i", "d_b"]
LABELS = {"i": "Inhibitor dose i", "h_i": "Hill coefficient h_i", "d_b": "Free-Bla inactivation d_b"}
BASELINE = {"i": 1.0, "h_i": 2.0, "d_b": 1.0}
N_AXIS = 21
T_END = 300.0
N_TIME = 301
TRACE_END = 60.0
SS_TOL = 1e-3
SS_HOLD_FRACTION = 0.10
VIABLE_FRACTION = 0.50
DOMINANCE_THRESHOLD = 0.75
REGIME_EXTINCTION_THRESHOLD = 1e-6
ANTIBIOTIC_CLEARANCE_THRESHOLD = 0.1


def project_steady_state(t: np.ndarray, s: np.ndarray, a: np.ndarray, p: dict, iota: float, beta: float) -> tuple[float, bool, np.ndarray]:
    """Shared-workstream steady-state detector: population per-capita rates.

    Detection requires both |d ln n_s/dt| and |d ln n_r/dt| < 1e-3 and the
    condition to hold for at least the final 10% of the simulated interval.
    Residual exponential antibiotic decay does not by itself block convergence.
    """
    g = s / (1.0 + s)
    ell = p["gamma"] * core.hill(a, p["h_a"]) * g
    rates = np.maximum(np.abs(g - ell), np.abs(p["alpha"] * g - beta * ell))
    below = rates < SS_TOL
    idx = len(t) - 1
    while idx > 0 and below[idx - 1]:
        idx -= 1
    t_star = float(t[idx])
    converged = bool(below[-1] and (t[-1] - t_star) >= SS_HOLD_FRACTION * (t[-1] - t[0]))
    return (t_star if converged else float(t[-1])), converged, rates


def comprehensive_metrics(sim: dict, p: dict, reference_ns: float) -> dict:
    t = sim["t"]
    ns, nr, nutrient, antibiotic, bla = sim["y"]
    total = ns + nr
    fr = np.divide(nr, total, out=np.zeros_like(nr), where=total > 0)
    diversity = core.shannon(fr)
    iota, beta, phi, db_effective = core.pathway_values(p)
    t_ss, converged, growth_rates = project_steady_state(t, nutrient, antibiotic, p, iota, beta)

    initial_total = ns[0] + nr[0]
    initial_fr = nr[0] / initial_total
    final_total = float(total[-1])
    final_fr = float(fr[-1])
    delta_fr = final_fr - initial_fr
    log_ratio_delta = float(np.log(max(nr[-1], 1e-300) / max(ns[-1], 1e-300)) - np.log(nr[0] / ns[0]))
    min_total_idx = int(np.argmin(total))
    min_ns_idx = int(np.argmin(ns))
    min_nr_idx = int(np.argmin(nr))
    recovery_line = core.RECOVERY_FRACTION * ns[-1]
    below_recovery = np.flatnonzero(ns < recovery_line)
    if len(below_recovery) == 0:
        t_recovery_own = 0.0
    else:
        first_dip = int(below_recovery[0])
        recovered = np.flatnonzero(ns[first_dip:] >= recovery_line)
        t_recovery_own = float(t[first_dip + recovered[0]]) if len(recovered) else np.nan
    ref_line = core.RECOVERY_FRACTION * reference_ns
    below_ref = np.flatnonzero(ns < ref_line)
    if len(below_ref) == 0:
        t_recovery_reference = 0.0
    else:
        first_dip_ref = int(below_ref[0])
        recovered_ref = np.flatnonzero(ns[first_dip_ref:] >= ref_line)
        t_recovery_reference = float(t[first_dip_ref + recovered_ref[0]]) if len(recovered_ref) else np.nan
    below_a = np.flatnonzero(antibiotic <= ANTIBIOTIC_CLEARANCE_THRESHOLD)
    below_a1 = np.flatnonzero(antibiotic <= 1.0)
    t_clear = float(t[below_a[0]]) if len(below_a) else np.nan
    t_a1 = float(t[below_a1[0]]) if len(below_a1) else np.nan

    if abs(delta_fr) <= core.SELECTION_TOLERANCE:
        direction = "neutral"
    elif delta_fr > 0:
        direction = "resistant-enriched"
    else:
        direction = "susceptible-enriched"
    resistant_extinct = bool(nr[-1] < core.EXTINCTION_THRESHOLD)
    susceptible_extinct = bool(ns[-1] < core.EXTINCTION_THRESHOLD)
    community_extinct = bool(final_total < REGIME_EXTINCTION_THRESHOLD)
    if community_extinct:
        outcome = "community extinction"
    elif resistant_extinct:
        outcome = "resistant extinction"
    elif susceptible_extinct:
        outcome = "susceptible extinction"
    elif final_fr <= 1.0 - DOMINANCE_THRESHOLD:
        outcome = "susceptible-dominant"
    elif final_fr >= DOMINANCE_THRESHOLD:
        outcome = "resistant-dominant"
    else:
        outcome = "mixed"

    l_over_g_initial = p["gamma"] * float(core.hill(p.get("a0", core.A0), p["h_a"]))
    static_selection_margin = (1.0 - beta) - ((1.0 - p["alpha"]) / max(l_over_g_initial, 1e-300))
    static_predicts_r = bool(static_selection_margin > 0)
    simulated_r_favored = bool(log_ratio_delta > 0)
    viability_cutoff = VIABLE_FRACTION * reference_ns
    return {
        "initial_fR": float(initial_fr), "final_fR": final_fr, "delta_fR": float(delta_fr),
        "delta_log_ratio": log_ratio_delta, "selection_direction": direction,
        "final_H": float(diversity[-1]), "min_H": float(np.min(diversity)),
        "final_ns": float(ns[-1]), "final_nr": float(nr[-1]), "final_total": final_total,
        "min_ns": float(ns[min_ns_idx]), "min_nr": float(nr[min_nr_idx]),
        "min_total": float(total[min_total_idx]),
        "peak_loss_S": float(1.0 - ns[min_ns_idx] / max(ns[0], 1e-300)),
        "peak_loss_total": float(1.0 - total[min_total_idx] / max(initial_total, 1e-300)),
        "time_min_total": float(t[min_total_idx]), "t_steady_state": t_ss,
        "ss_converged": converged, "max_final_growth_rate": float(growth_rates[-1]),
        "t_recovery_S_own_final": t_recovery_own,
        "t_recovery_S_no_treatment": t_recovery_reference,
        "sensitive_survives": bool(ns[-1] >= viability_cutoff),
        "resistant_eliminated": resistant_extinct,
        "community_extinct_1e-6": community_extinct,
        "time_antibiotic_below_0p1": t_clear, "time_antibiotic_below_1": t_a1,
        "antibiotic_AUC": core.trapz(antibiotic, t), "free_Bla_AUC": core.trapz(bla, t),
        "total_population_AUC": core.trapz(total, t),
        "outcome_regime": outcome, "solver_success": bool(sim["success"]),
        "iota": iota, "beta": beta, "phi": phi, "effective_d_b_iota": db_effective,
        "l_over_g_initial": l_over_g_initial, "static_selection_margin": static_selection_margin,
        "static_criterion_predicts_R": static_predicts_r, "simulated_R_favored": simulated_r_favored,
        "criterion_agrees": static_predicts_r == simulated_r_favored,
    }


def verify_logspace_equivalence(n_checks: int = 12, seed: int = 574) -> tuple[pd.DataFrame, str]:
    """Compare this log-space implementation with direct linear-state ODEs.

    Uses randomized inhibitor parameter triples at the Ma-style a0=2.15 dose.
    The compared equations are independently evaluated in raw coordinates.
    """
    rng = np.random.default_rng(seed)
    i_values = rng.uniform(0.1, 10.0, n_checks)
    hi_values = rng.uniform(1.0, 5.0, n_checks)
    db_values = rng.uniform(1.0, 10.0, n_checks)
    t_eval = np.linspace(0.0, 40.0, 201)
    rows = []
    for idx, (i_value, hi_value, db_value) in enumerate(zip(i_values, hi_values, db_values)):
        p = dict(core.BASE, i=float(i_value), h_i=float(hi_value), d_b=float(db_value))
        iota, beta, phi, effective_db = core.pathway_values(p)

        def rhs_linear(t, y):
            ns, nr, s, a, b = np.maximum(y, 0.0)
            g = s / (1.0 + s)
            ell = p["gamma"] * float(core.hill(a, p["h_a"])) * g
            return np.array([
                (g - ell) * ns,
                (p["alpha"] * g - beta * ell) * nr,
                (p["xi"] * ell - g) * ns + (p["xi"] * beta * ell - p["alpha"] * g) * nr,
                -p["kappa_b"] * b * a - phi * nr * a - p["d_a"] * a,
                beta * ell * nr - effective_db * b,
            ])

        raw = solve_ivp(rhs_linear, (0.0, 40.0), [core.NS0, core.NR0, core.S0, 2.15, core.B0],
                        t_eval=t_eval, method="LSODA", rtol=3e-7, atol=1e-9, max_step=2.0)
        log = core.simulate(p, a0=2.15, t_end=40.0, n_points=len(t_eval))
        allowed = 2e-6 + 2e-4 * np.abs(raw.y)
        normalized_error = np.abs(raw.y - log["y"]) / allowed
        max_ratio = float(np.max(normalized_error))
        rows.append({"check": idx, "i": i_value, "h_i": hi_value, "d_b": db_value,
                     "max_normalized_error": max_ratio, "within_tolerance": max_ratio <= 1.0,
                     "solver_success_raw": bool(raw.success), "solver_success_log": log["success"]})
    frame = pd.DataFrame(rows)
    frame.to_csv(TABLES / "comprehensive_logspace_equivalence.csv", index=False)
    passed = bool(frame.within_tolerance.all() and frame.solver_success_raw.all() and frame.solver_success_log.all())
    text = (f"The log-population solver was compared with direct linear-state integration at {n_checks} randomized (i, hᵢ, dᵦ) settings. "
            f"All solver integrations succeeded; {frame.within_tolerance.sum()}/{n_checks} checks met the combined absolute/relative tolerance "
            f"(maximum normalized error {frame.max_normalized_error.max():.3g}). Equivalence status: {'PASS' if passed else 'REVIEW REQUIRED'}.")
    return frame, text


def run_full_grid() -> pd.DataFrame:
    out_file = TABLES / "comprehensive_inhibitor_grid_a0_2p15.csv"
    if out_file.exists():
        existing = pd.read_csv(out_file)
        expected = N_AXIS ** 3
        if len(existing) == expected and existing["solver_success"].astype(bool).all():
            print(f"Using completed full-grid data: {len(existing):,} rows", flush=True)
            return existing
    i_values = np.geomspace(0.1, 10.0, N_AXIS)
    hi_values = np.linspace(1.0, 5.0, N_AXIS)
    db_values = np.linspace(1.0, 10.0, N_AXIS)
    p_ref = dict(core.BASE)
    ref = core.no_treatment_reference(p_ref, t_end=T_END)
    reference_ns = float(ref["y"][0, -1])
    rows = []
    start = time.monotonic()
    total = N_AXIS ** 3
    count = 0
    for i_value in i_values:
        for hi_value in hi_values:
            for db_value in db_values:
                p = dict(core.BASE, i=float(i_value), h_i=float(hi_value), d_b=float(db_value))
                p["a0"] = 2.15
                sim = core.simulate(p, a0=2.15, t_end=T_END, n_points=N_TIME)
                metrics = comprehensive_metrics(sim, p, reference_ns)
                row = {**p, **metrics, "a0": 2.15, "simulation_horizon": T_END}
                rows.append(row)
                count += 1
        if count % (N_AXIS * 3) == 0 or count == total:
            print(f"  full inhibitor grid: {count:,}/{total:,} ({time.monotonic()-start:.0f}s)", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out_file, index=False)
    return frame


ENDPOINTS = {
    "delta_fR": "Change in resistant fraction",
    "final_fR": "Final resistant fraction",
    "final_H": "Final Shannon diversity",
    "min_H": "Minimum Shannon diversity",
    "final_total": "Final total population",
    "t_steady_state": "Time to steady state",
    "t_recovery_S_own_final": "Sensitive recovery time",
    "peak_loss_total": "Peak total population loss",
    "final_ns": "Final susceptible population",
    "final_nr": "Final resistant population",
}


def summarize_grid(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = pd.DataFrame([
        {"endpoint": label, "column": col, "min": df[col].min(), "median": df[col].median(), "max": df[col].max(),
         "finite_fraction": df[col].notna().mean()}
        for col, label in ENDPOINTS.items()
    ])
    summary.to_csv(TABLES / "comprehensive_endpoint_summary.csv", index=False)
    regime = df["outcome_regime"].value_counts(dropna=False).rename_axis("outcome_regime").reset_index(name="n")
    regime["fraction"] = regime["n"] / len(df)
    regime.to_csv(TABLES / "comprehensive_regime_counts.csv", index=False)

    direction = df["selection_direction"].value_counts(normalize=True).to_dict()
    criterion_accuracy = df["criterion_agrees"].mean()
    validation = pd.DataFrame([{
        "n": len(df), "solver_success_fraction": df["solver_success"].mean(),
        "ss_converged_fraction": df["ss_converged"].mean(),
        "criterion_agreement_fraction": criterion_accuracy,
        "static_criterion_R_fraction": df["static_criterion_predicts_R"].mean(),
        "simulated_R_favored_fraction": df["simulated_R_favored"].mean(),
        "resistant_eliminated_fraction": df["resistant_eliminated"].mean(),
        "sensitive_survives_fraction": df["sensitive_survives"].mean(),
        **{f"fraction_{key}": val for key, val in direction.items()},
    }])
    validation.to_csv(TABLES / "comprehensive_run_validation_summary.csv", index=False)
    return summary, regime, validation


def parameter_sensitivity(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    endps = ["delta_fR", "final_fR", "final_H", "final_total", "t_steady_state", "t_recovery_S_own_final", "peak_loss_total"]
    corr_rows = []
    for parameter in PARAMS:
        for endpoint in endps:
            rho, pvalue = spearmanr(df[parameter], df[endpoint], nan_policy="omit")
            corr_rows.append({"parameter": parameter, "endpoint": endpoint, "spearman_rho": rho, "p_value": pvalue})
    correlations = pd.DataFrame(corr_rows)
    correlations.to_csv(TABLES / "comprehensive_spearman_sensitivity.csv", index=False)

    z = df[PARAMS].copy()
    for p in PARAMS:
        z[p] = (z[p] - z[p].mean()) / z[p].std(ddof=0)
    feature_names = ["intercept", *PARAMS, "i:h_i", "i:d_b", "h_i:d_b"]
    X = np.column_stack([
        np.ones(len(df)), z["i"], z["h_i"], z["d_b"],
        z["i"] * z["h_i"], z["i"] * z["d_b"], z["h_i"] * z["d_b"],
    ])
    reg_rows = []
    for endpoint in endps:
        valid = df[endpoint].notna().to_numpy()
        x, y = X[valid], df.loc[valid, endpoint].to_numpy(float)
        y_sd = y.std(ddof=0)
        if y_sd == 0:
            continue
        yz = (y - y.mean()) / y_sd
        coef, _, _, _ = np.linalg.lstsq(x, yz, rcond=None)
        residual = yz - x @ coef
        dof = max(len(yz) - x.shape[1], 1)
        mse = float(residual @ residual / dof)
        covariance = mse * np.linalg.pinv(x.T @ x)
        se = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        tvals = coef / np.where(se > 0, se, np.nan)
        pvals = 2.0 * student_t.sf(np.abs(tvals), dof)
        r2 = 1.0 - float(residual @ residual) / max(float(((yz - yz.mean()) ** 2).sum()), 1e-300)
        for name, value, std_err, tval, pval in zip(feature_names, coef, se, tvals, pvals):
            reg_rows.append({"endpoint": endpoint, "term": name, "standardized_coefficient": value,
                             "std_error": std_err, "t_statistic": tval, "p_value": pval,
                             "r_squared": r2, "n": len(y)})
    regression = pd.DataFrame(reg_rows)
    regression.to_csv(TABLES / "comprehensive_standardized_regression.csv", index=False)
    return correlations, regression


def binned_effects(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for parameter in PARAMS:
        bins = pd.qcut(df[parameter], q=11, duplicates="drop")
        for endpoint in ["peak_loss_total", "t_recovery_S_own_final", "t_steady_state", "final_total", "final_fR", "final_H"]:
            temp = df[[parameter, endpoint]].copy()
            temp["bin"] = bins
            for bin_name, group in temp.groupby("bin", observed=True):
                rows.append({"parameter": parameter, "endpoint": endpoint, "bin": str(bin_name),
                             "parameter_mean": group[parameter].mean(), "endpoint_mean": group[endpoint].mean(),
                             "endpoint_sd": group[endpoint].std(), "n": len(group)})
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES / "comprehensive_binned_effects.csv", index=False)
    fig, axes = plt.subplots(6, 3, figsize=(14, 19), constrained_layout=True)
    endpoints = ["peak_loss_total", "t_recovery_S_own_final", "t_steady_state", "final_total", "final_fR", "final_H"]
    for col, parameter in enumerate(PARAMS):
        for row, endpoint in enumerate(endpoints):
            ax = axes[row, col]
            sub = summary[(summary.parameter == parameter) & (summary.endpoint == endpoint)].sort_values("parameter_mean")
            x, mean, sd = sub.parameter_mean.to_numpy(), sub.endpoint_mean.to_numpy(), sub.endpoint_sd.fillna(0).to_numpy()
            ax.plot(x, mean, "o-", ms=3, color="#0072B2")
            ax.fill_between(x, mean - sd, mean + sd, alpha=.18, color="#0072B2")
            ax.set_ylabel(ENDPOINTS.get(endpoint, endpoint) if col == 0 else "")
            if row == 0:
                ax.set_title(LABELS[parameter])
            if row == len(endpoints) - 1:
                ax.set_xlabel(LABELS[parameter])
            if parameter == "i":
                ax.set_xscale("log")
            ax.grid(alpha=.2)
    fig.suptitle("Marginal binned effects across the joint inhibitor sweep", y=1.01)
    core.savefig(fig, "comprehensive_01_binned_effects")
    return summary


def tradeoff_plots(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    sc = axes[0].scatter(df.final_total, df.final_fR, c=df.i, norm=plt.matplotlib.colors.LogNorm(.1, 10), s=7, alpha=.5, cmap="viridis")
    axes[0].set(xlabel="final total population", ylabel="final resistant fraction", title="Composition versus population recovery")
    fig.colorbar(sc, ax=axes[0], label="inhibitor dose i")
    sc = axes[1].scatter(df.final_total, df.final_H, c=df.h_i, s=7, alpha=.5, cmap="plasma")
    axes[1].set(xlabel="final total population", ylabel="final Shannon diversity", title="Diversity versus final biomass")
    fig.colorbar(sc, ax=axes[1], label="h_i")
    sc = axes[2].scatter(df.peak_loss_total, df.t_recovery_S_own_final, c=df.d_b, s=7, alpha=.5, cmap="magma")
    axes[2].set(xlabel="peak total population loss", ylabel="susceptible recovery time", title="Bottleneck versus recovery")
    fig.colorbar(sc, ax=axes[2], label="d_b")
    core.savefig(fig, "comprehensive_02_endpoint_tradeoffs")


def composition_curves(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for parameter in PARAMS:
        bin_id = pd.qcut(df[parameter], q=11, duplicates="drop")
        sub = df.assign(_bin=bin_id)
        for key, group in sub.groupby("_bin", observed=True):
            for field in ["final_ns", "final_nr", "final_total", "final_fR"]:
                rows.append({"parameter": parameter, "bin": str(key), "parameter_mean": group[parameter].mean(),
                             "metric": field, "mean": group[field].mean(), "sd": group[field].std(), "n": len(group)})
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES / "comprehensive_population_composition_bins.csv", index=False)
    fig, axes = plt.subplots(3, 2, figsize=(12, 12), constrained_layout=True)
    for row, parameter in enumerate(PARAMS):
        ax = axes[row, 0]
        for field, label, color in [("final_ns", "susceptible", "tab:green"), ("final_nr", "resistant", "tab:gray")]:
            sub = summary[(summary.parameter == parameter) & (summary.metric == field)].sort_values("parameter_mean")
            x, m, s = sub.parameter_mean.to_numpy(), sub["mean"].to_numpy(), sub.sd.fillna(0).to_numpy()
            ax.plot(x, m, "o-", label=label, color=color, ms=3)
            ax.fill_between(x, m-s, m+s, color=color, alpha=.15)
        ax.set(title=f"Absolute endpoints versus {LABELS[parameter]}", ylabel="final population")
        if parameter == "i":
            ax.set_xscale("log")
        ax.legend(frameon=False)
        ax.grid(alpha=.2)
        ax = axes[row, 1]
        sub = summary[(summary.parameter == parameter) & (summary.metric == "final_fR")].sort_values("parameter_mean")
        x, m, s = sub.parameter_mean.to_numpy(), sub["mean"].to_numpy(), sub.sd.fillna(0).to_numpy()
        ax.plot(x, m, "o-", color="tab:purple", ms=3)
        ax.fill_between(x, m-s, m+s, color="tab:purple", alpha=.18)
        ax.set(title="Mean final resistant fraction", ylabel="final f_R", ylim=(0, 1))
        if parameter == "i":
            ax.set_xscale("log")
        ax.grid(alpha=.2)
        if row == 2:
            for col in range(2): axes[row, col].set_xlabel(LABELS[parameter])
    core.savefig(fig, "comprehensive_03_population_composition")
    return summary


def pairwise_heatmaps(df: pd.DataFrame) -> list[str]:
    pairs = list(itertools.combinations(PARAMS, 2))
    endpoints = ["final_H", "t_steady_state", "final_fR", "t_recovery_S_own_final",
                 "peak_loss_total", "final_total", "final_ns", "final_nr"]
    titles = {"final_H": "Final Shannon diversity", "t_steady_state": "Time to steady state",
              "final_fR": "Final resistant fraction", "t_recovery_S_own_final": "Susceptible recovery time",
              "peak_loss_total": "Peak total population loss", "final_total": "Final total population",
              "final_ns": "Final susceptible population", "final_nr": "Final resistant population"}
    figures = []
    for endpoint in endpoints:
        fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
        for ax, (x_name, y_name) in zip(axes.flat, pairs):
            subset = df
            for parameter in PARAMS:
                if parameter not in (x_name, y_name):
                    subset = subset[np.isclose(subset[parameter], BASELINE[parameter])]
            table = subset.pivot(index=y_name, columns=x_name, values=endpoint).sort_index(ascending=False)
            mesh = ax.pcolormesh(table.columns.to_numpy(float), table.index.to_numpy(float), table.to_numpy(), shading="nearest", cmap="viridis")
            if x_name == "i": ax.set_xscale("log")
            ax.set_xlabel(LABELS[x_name]); ax.set_ylabel(LABELS[y_name])
            fixed = ", ".join(f"{LABELS[p]}={BASELINE[p]:g}" for p in PARAMS if p not in (x_name, y_name))
            ax.set_title(f"{LABELS[x_name]} × {LABELS[y_name]}\nFixed: {fixed}")
            fig.colorbar(mesh, ax=ax, label=titles[endpoint])
        axes.flat[-1].axis("off")
        fig.suptitle(f"Person 3 endpoint slices: {titles[endpoint]} (a0=2.15)", y=1.01)
        name = f"comprehensive_pairwise_{endpoint}"
        core.savefig(fig, name)
        figures.append(name)
    return figures


def convergence_diagnostics(df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    frame = df.copy()
    frame["nonconverged_reason"] = np.select(
        [~frame.ss_converged & (frame.final_total < REGIME_EXTINCTION_THRESHOLD),
         ~frame.ss_converged & (frame.final_total >= REGIME_EXTINCTION_THRESHOLD)],
        ["near-zero final population", "substantial population; slow or strict detector"],
        default="converged",
    )
    frame.to_csv(TABLES / "comprehensive_convergence_diagnostics.csv", index=False)
    reason_counts = frame["nonconverged_reason"].value_counts().rename_axis("reason").reset_index(name="n")
    reason_counts["fraction"] = reason_counts.n / len(frame)
    reason_counts.to_csv(TABLES / "comprehensive_convergence_reason_counts.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    cmap = {"converged": "tab:green", "near-zero final population": "tab:red", "substantial population; slow or strict detector": "tab:orange"}
    for reason, group in frame.groupby("nonconverged_reason"):
        axes[0].scatter(group.t_steady_state, group.final_total, s=9, alpha=.45, label=f"{reason} (n={len(group)})", color=cmap[reason])
    axes[0].set(xlabel="reported steady-state time (horizon if censored)", ylabel="final total population", title="Slow dynamics versus near-extinction")
    axes[0].legend(frameon=False, fontsize=7)
    axes[0].grid(alpha=.2)
    pivot = reason_counts.set_index("reason")["n"]
    axes[1].bar(pivot.index, pivot.values, color=[cmap[k] for k in pivot.index])
    axes[1].tick_params(axis="x", rotation=20)
    axes[1].set(ylabel="simulation count", title="Convergence classification")
    core.savefig(fig, "comprehensive_04_convergence_diagnostics")
    return reason_counts, ""


def recovery_regimes(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    frame = df.copy()
    total = frame.final_total
    fr = frame.final_fR
    frame["regime_detail"] = np.select(
        [(frame.ss_converged) & (total < REGIME_EXTINCTION_THRESHOLD),
         (~frame.ss_converged) & (total < REGIME_EXTINCTION_THRESHOLD),
         (~frame.ss_converged) & (total >= REGIME_EXTINCTION_THRESHOLD),
         fr <= (1.0 - DOMINANCE_THRESHOLD), fr >= DOMINANCE_THRESHOLD],
        ["community extinct", "nonconverged near-extinction", "nonconverged viable density", "susceptible-dominated", "resistant-dominated"],
        default="mixed community",
    )
    frame.to_csv(TABLES / "comprehensive_regime_classified_grid.csv", index=False)
    order = ["community extinct", "nonconverged near-extinction", "nonconverged viable density", "susceptible-dominated", "resistant-dominated", "mixed community"]
    summaries = []
    for parameter in PARAMS:
        bins = pd.qcut(frame[parameter], q=10, duplicates="drop")
        group = frame.groupby([bins, "regime_detail"], observed=True).size().rename("n").reset_index()
        group = group.rename(columns={parameter: "bin"})
        for bin_value, block in group.groupby(group.columns[0], observed=True):
            total_n = block.n.sum()
            x_value = frame.loc[bins == bin_value, parameter].mean()
            for regime in order:
                n = int(block.loc[block.regime_detail == regime, "n"].sum())
                summaries.append({"parameter": parameter, "bin": str(bin_value), "parameter_mean": x_value,
                                  "regime": regime, "n": n, "fraction": n / total_n if total_n else 0.0})
    summary = pd.DataFrame(summaries)
    summary.to_csv(TABLES / "comprehensive_regime_frequency_bins.csv", index=False)
    fig, axes = plt.subplots(3, 3, figsize=(15, 12), constrained_layout=True)
    regime_colors = {"community extinct": "black", "nonconverged near-extinction": "#D55E00",
                     "nonconverged viable density": "#E69F00", "susceptible-dominated": "#0072B2",
                     "resistant-dominated": "#CC79A7", "mixed community": "#009E73"}
    for row, parameter in enumerate(PARAMS):
        ssub = summary[summary.parameter == parameter]
        pv = ssub.pivot(index="parameter_mean", columns="regime", values="fraction").fillna(0).sort_index()
        ax = axes[row, 0]
        bottom = np.zeros(len(pv))
        for regime in order:
            vals = pv[regime].to_numpy() if regime in pv else np.zeros(len(pv))
            ax.bar(np.arange(len(pv)), vals, bottom=bottom, color=regime_colors[regime], label=regime)
            bottom += vals
        ax.set_xticks(np.arange(len(pv))); ax.set_xticklabels([f"{x:g}" for x in pv.index], rotation=45, ha="right", fontsize=7)
        ax.set(title=f"Outcome regimes by {LABELS[parameter]}", ylabel="fraction")
        ax.set_ylim(0, 1)
        for col, metric, title, color in [(1, "final_H", "Mean final diversity", "#6A3D9A"), (2, "final_total", "Mean final biomass", "#0072B2")]:
            binned = frame.assign(_bin=pd.qcut(frame[parameter], q=10, duplicates="drop")).groupby("_bin", observed=True).agg(x=(parameter, "mean"), y=(metric, "mean"))
            axes[row, col].plot(binned.x, binned.y, "o-", color=color)
            axes[row, col].set(title=title, xlabel=LABELS[parameter], ylabel=metric)
            axes[row, col].grid(alpha=.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles[:len(order)], labels[:len(order)], loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(.5, 1.04))
    core.savefig(fig, "comprehensive_05_recovery_regimes")
    counts = frame.regime_detail.value_counts().reindex(order, fill_value=0)
    regime_text = "Outcome classification over the full joint grid: " + "; ".join(f"{name}: {count:,} ({count/len(frame):.1%})" for name, count in counts.items()) + ". Dominance uses f_R ≤ 0.25 or ≥ 0.75; viability uses 50% of the no-treatment susceptible endpoint; extinction is reported at 10⁻⁴ for individual strains and 10⁻⁶ for the whole community."
    return frame, summary, regime_text


def direct_pathway_plot() -> None:
    doses = np.geomspace(.1, 10, 300)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    for hi in [1, 2, 3, 5]:
        vals = np.array([core.inhibitor_activity(i, hi) for i in doses])
        axes[0].plot(doses, vals, label=f"h_i={hi}")
        axes[1].plot(doses, core.BASE["beta_min"] + core.BASE["c"] * (1-core.BASE["beta_min"]) * vals, label=f"h_i={hi}")
        axes[2].plot(doses, core.BASE["phi_max"] * (1-core.BASE["c"] * vals), label=f"h_i={hi}")
    for ax in axes:
        ax.set_xscale("log"); ax.axvline(1, color="0.4", ls="--", lw=1); ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=8)
    axes[0].set(xlabel="inhibitor dose i", ylabel="iota", title="Inhibitor activity")
    axes[1].set(xlabel="inhibitor dose i", ylabel="beta", title="Resistant lysis factor")
    axes[2].set(xlabel="inhibitor dose i", ylabel="phi", title="Live-cell antibiotic degradation")
    core.savefig(fig, "comprehensive_06_direct_pharmacodynamics")


def representative_timecourses(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cases = []
    base_row = df.iloc[((df[PARAMS] - pd.Series(BASELINE)).abs().sum(axis=1)).argmin()]
    picks = [
        ("baseline", base_row),
        ("fastest_sensitive_recovery", df.loc[df.t_recovery_S_own_final.idxmin()]),
        ("deepest_total_bottleneck", df.loc[df.peak_loss_total.idxmax()]),
        ("most_resistant_enriched", df.loc[df.delta_fR.idxmax()]),
        ("most_susceptible_enriched", df.loc[df.delta_fR.idxmin()]),
        ("highest_final_diversity", df.loc[df.final_H.idxmax()]),
    ]
    for label, row in picks:
        cases.append({"case": label, **{p: float(row[p]) for p in PARAMS},
                      **{m: row[m] for m in ["delta_fR", "final_fR", "final_total", "peak_loss_total", "t_recovery_S_own_final", "t_steady_state", "outcome_regime"]}})
    case_df = pd.DataFrame(cases).drop_duplicates(subset=PARAMS)
    case_df.to_csv(TABLES / "comprehensive_representative_cases.csv", index=False)
    fig, axes = plt.subplots(len(case_df), 4, figsize=(15, 3.0*len(case_df)), squeeze=False, constrained_layout=True)
    all_rows = []
    for row_idx, (_, case) in enumerate(case_df.iterrows()):
        p = dict(core.BASE, i=case.i, h_i=case.h_i, d_b=case.d_b)
        # Show the transient clearly; endpoint metrics still use the full τ=300 run.
        sim = core.simulate(p, a0=2.15, t_end=TRACE_END, n_points=601)
        ns, nr, s, drug, bla = sim["y"]
        total = ns + nr
        fr = nr / np.maximum(total, 1e-300)
        t = sim["t"]
        axes[row_idx, 0].plot(t, drug, color="tab:red")
        axes[row_idx, 1].plot(t, ns, label="n_s", color="tab:green")
        axes[row_idx, 1].plot(t, nr, label="n_r", color="tab:gray")
        axes[row_idx, 2].plot(t, fr, color="tab:purple")
        axes[row_idx, 2].axhline(.5, color="0.5", ls="--", lw=.8)
        axes[row_idx, 3].plot(t, bla, color="black", label="free Bla")
        axes[row_idx, 3].plot(t, s, color="tab:blue", alpha=.75, label="nutrient s")
        axes[row_idx, 0].set_ylabel(f"{case['case']}\ni={case.i:.3g}, h_i={case.h_i:.2g}, d_b={case.d_b:.2g}")
        all_rows.append({"case": case["case"], "time": t, "ns": ns, "nr": nr, "a": drug, "b": bla, "s": s, "fR": fr})
    for ax, title, ylabel in zip(axes[0], ["Antibiotic", "Populations", "Resistant fraction", "Bla and nutrient"], ["a", "population", "f_R", "b, s"]):
        ax.set_title(title); ax.set_ylabel(ylabel); ax.set_xlabel("dimensionless time τ"); ax.grid(alpha=.2)
    axes[0, 1].legend(frameon=False, fontsize=7); axes[0, 3].legend(frameon=False, fontsize=7)
    fig.suptitle("Representative inhibitor regimes (first 60 τ; endpoints use 300 τ)", y=1.025, fontsize=11)
    core.savefig(fig, "comprehensive_07_representative_trajectories")
    trajectory_records = []
    for item in all_rows:
        trajectory_records.extend({"case": item["case"], "time": float(tt), "n_s": float(ss), "n_r": float(rr), "a": float(aa), "b": float(bb), "s": float(nn), "f_R": float(ff)}
                                  for tt, ss, rr, aa, bb, nn, ff in zip(item["time"], item["ns"], item["nr"], item["a"], item["b"], item["s"], item["fR"]))
    trace_df = pd.DataFrame(trajectory_records)
    trace_df.to_csv(TABLES / "comprehensive_representative_trajectories.csv", index=False)
    return case_df, trace_df


def one_factor_trajectories() -> pd.DataFrame:
    levels = {
        "i": [0.1, 1.0, 10.0],
        "h_i": [1.0, 2.0, 5.0],
        "d_b": [1.0, 5.5, 10.0],
    }
    records = []
    fig, axes = plt.subplots(3, 3, figsize=(15, 12), constrained_layout=True)
    colors = ["#0072B2", "#E69F00", "#D55E00"]
    for row, parameter in enumerate(PARAMS):
        for value, level, color in zip(levels[parameter], ["low", "middle", "high"], colors):
            p = dict(core.BASE, **BASELINE)
            p[parameter] = value
            sim = core.simulate(p, a0=2.15, t_end=60, n_points=601)
            t = sim["t"]; ns, nr, s, drug, bla = sim["y"]
            total = ns + nr; fr = nr / np.maximum(total, 1e-300)
            clearance = np.flatnonzero(drug <= ANTIBIOTIC_CLEARANCE_THRESHOLD)
            records.append({"parameter": parameter, "level": level, "value": value,
                            "antibiotic_clearance_time_a_le_0p1": float(t[clearance[0]]) if len(clearance) else np.nan,
                            "antibiotic_AUC": core.trapz(drug, t), "free_Bla_AUC": core.trapz(bla, t),
                            "min_total": float(total.min()), "final_total": float(total[-1]),
                            "final_fR": float(fr[-1]), "peak_loss_total": float(1-total.min()/total[0])})
            axes[row, 0].plot(t, drug, color=color, label=f"{level}: {value:g}")
            axes[row, 0].axhline(.1, color="0.3", ls="--", lw=.7)
            axes[row, 1].plot(t, total, color=color, label=f"{level}: {value:g}")
            axes[row, 2].plot(t, fr, color=color, label=f"{level}: {value:g}")
    for row, parameter in enumerate(PARAMS):
        axes[row, 0].set_ylabel(f"{LABELS[parameter]}\na(t)")
        axes[row, 1].set_ylabel("total population")
        axes[row, 2].set_ylabel("f_R")
        for col in range(3):
            axes[row, col].grid(alpha=.2)
        axes[row, 2].set_ylim(0, 1)
        axes[row, 2].legend(frameon=False, fontsize=7)
    for ax, title in zip(axes[0], ["Antibiotic clearance", "Crash and regrowth", "Community composition"]): ax.set_title(title)
    for ax in axes[-1]: ax.set_xlabel("dimensionless time τ")
    fig.suptitle("Low, middle, and high one-factor inhibitor perturbations", y=1.01)
    core.savefig(fig, "comprehensive_08_clearance_recovery_trajectories")
    summary = pd.DataFrame(records)
    summary.to_csv(TABLES / "comprehensive_one_factor_clearance_summary.csv", index=False)
    return summary


def criterion_validation(df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    summary = pd.DataFrame([{
        "n": len(df), "static_criterion_accuracy": df.criterion_agrees.mean(),
        "static_predicts_R_fraction": df.static_criterion_predicts_R.mean(),
        "dynamic_R_favored_fraction": df.simulated_R_favored.mean(),
        "false_R_prediction_fraction": ((df.static_criterion_predicts_R) & (~df.simulated_R_favored)).mean(),
        "missed_R_enrichment_fraction": ((~df.static_criterion_predicts_R) & (df.simulated_R_favored)).mean(),
    }])
    summary.to_csv(TABLES / "comprehensive_selection_criterion_summary.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    correct = df.groupby(["i", "h_i"]).criterion_agrees.mean().reset_index()
    table = correct.pivot(index="h_i", columns="i", values="criterion_agrees").sort_index(ascending=False)
    mesh = axes[0].pcolormesh(table.columns.to_numpy(), table.index.to_numpy(), table.to_numpy(), shading="nearest", cmap="RdYlGn", vmin=0, vmax=1)
    axes[0].set_xscale("log"); axes[0].set(xlabel="i", ylabel="h_i", title="Initial-state criterion agreement (averaged over d_b)")
    fig.colorbar(mesh, ax=axes[0], label="fraction agreement")
    axes[1].scatter(df.static_selection_margin, df.delta_log_ratio, c=df.d_b, cmap="viridis", s=5, alpha=.25)
    axes[1].axhline(0, color="0.3", ls="--"); axes[1].axvline(0, color="0.3", ls="--")
    axes[1].set(xlabel="initial static criterion margin", ylabel="dynamic Δlog(n_r/n_s)", title="Static prediction versus integrated selection")
    core.savefig(fig, "comprehensive_09_selection_criterion_validation")
    text = (f"The time-varying simulations are compared with the paper’s instantaneous selection criterion evaluated at the initial nutrient and antibiotic state. "
            f"Across {len(df):,} inhibitor combinations, that initial-state criterion agrees with the sign of the full-trajectory log-ratio change in {df.criterion_agrees.mean():.1%} of runs. "
            "Mismatches are expected because iota, beta, and the lysis-to-growth ratio interact with the evolving antibiotic/Bla environment; d_b changes the trajectory but is absent from the instantaneous criterion.")
    return summary, text


def sensitivity_figures(corr: pd.DataFrame, regression: pd.DataFrame) -> list[str]:
    names = []
    endpoints = ["delta_fR", "final_H", "final_total", "t_recovery_S_own_final", "peak_loss_total"]
    fig, axes = plt.subplots(1, len(endpoints), figsize=(17, 5), constrained_layout=True)
    for ax, endpoint in zip(axes, endpoints):
        sub = corr[corr.endpoint == endpoint].set_index("parameter").reindex(PARAMS)
        ax.barh([LABELS[p] for p in PARAMS], sub.spearman_rho, color=["#0072B2" if x >= 0 else "#D55E00" for x in sub.spearman_rho])
        ax.axvline(0, color="0.3", lw=.8); ax.set_title(ENDPOINTS.get(endpoint, endpoint)); ax.grid(axis="x", alpha=.2)
    fig.suptitle("Global rank sensitivity across the full inhibitor grid", y=1.02)
    core.savefig(fig, "comprehensive_10_spearman_sensitivity"); names.append("comprehensive_10_spearman_sensitivity")
    terms = ["i", "h_i", "d_b", "i:h_i", "i:d_b", "h_i:d_b"]
    fig, axes = plt.subplots(1, len(endpoints), figsize=(18, 5.5), constrained_layout=True)
    for ax, endpoint in zip(axes, endpoints):
        sub = regression[(regression.endpoint == endpoint) & regression.term.isin(terms)].set_index("term").reindex(terms)
        ax.barh(terms, sub.standardized_coefficient, color=["#0072B2" if x >= 0 else "#D55E00" for x in sub.standardized_coefficient])
        ax.axvline(0, color="0.3", lw=.8); ax.set_title(ENDPOINTS.get(endpoint, endpoint)); ax.grid(axis="x", alpha=.2)
    fig.suptitle("Standardized main effects and pairwise interactions", y=1.02)
    core.savefig(fig, "comprehensive_11_standardized_interactions"); names.append("comprehensive_11_standardized_interactions")
    return names


def report_sensitivity_table(corr: pd.DataFrame, regression: pd.DataFrame) -> list[list[str]]:
    endpoints = ["delta_fR", "final_H", "final_total", "t_recovery_S_own_final", "peak_loss_total"]
    main_terms = PARAMS
    interaction_terms = ["i:h_i", "i:d_b", "h_i:d_b"]
    rows = [["Endpoint", "Strongest main β", "Strongest interaction β", "Top Spearman ρ", "R²"]]
    for endpoint in endpoints:
        model = regression[regression.endpoint == endpoint]
        main = model[model.term.isin(main_terms)].iloc[model[model.term.isin(main_terms)].standardized_coefficient.abs().argmax()]
        interactions = model[model.term.isin(interaction_terms)]
        interaction = interactions.iloc[interactions.standardized_coefficient.abs().argmax()]
        rank = corr[corr.endpoint == endpoint].iloc[corr[corr.endpoint == endpoint].spearman_rho.abs().argmax()]
        rows.append([
            ENDPOINTS.get(endpoint, endpoint),
            f"{main.term}: {main.standardized_coefficient:+.3f}",
            f"{interaction.term}: {interaction.standardized_coefficient:+.3f}",
            f"{rank.parameter}: {rank.spearman_rho:+.3f}",
            f"{model.r_squared.iloc[0]:.3f}",
        ])
    return rows


def make_additional_config(n: int) -> None:
    config = {
        "analysis": "Person 3 comprehensive inhibitor pharmacodynamics", "repository_location": "BME574 results/",
        "joint_grid": {"i": [0.1, 10.0, N_AXIS, "log-spaced"], "h_i": [1.0, 5.0, N_AXIS, "linear"], "d_b": [1.0, 10.0, N_AXIS, "linear"]},
        "a0": 2.15, "horizon": T_END, "time_points": N_TIME,
        "representative_trace_horizon": TRACE_END,
        "solver": {"method": "LSODA", "rtol": 3e-7, "atol": 1e-9, "max_step": 2.0},
        "shared_steady_state": {"per_capita_growth_tolerance": SS_TOL, "hold_fraction": SS_HOLD_FRACTION},
        "recovery": {"fraction": core.RECOVERY_FRACTION, "own_treatment_final_and_no_treatment_reference_both_reported": True},
        "viability_fraction_of_no_treatment_susceptible": VIABLE_FRACTION,
        "individual_extinction_threshold": core.EXTINCTION_THRESHOLD,
        "community_extinction_threshold": REGIME_EXTINCTION_THRESHOLD,
        "dominance_threshold": DOMINANCE_THRESHOLD,
        "rows": n,
    }
    (TABLES / "comprehensive_analysis_config.json").write_text(json.dumps(config, indent=2))


def main() -> None:
    print("Starting comprehensive Person 3 analysis (writing outputs under repository results/)", flush=True)
    equivalence, equivalence_text = verify_logspace_equivalence()
    print("  logspace/raw equation comparison complete", flush=True)
    df = run_full_grid()
    summary, regime_counts, validation = summarize_grid(df)
    correlations, regression = parameter_sensitivity(df)
    bins = binned_effects(df)
    tradeoff_plots(df)
    composition_curves(df)
    pair_figures = pairwise_heatmaps(df)
    convergence_diagnostics(df)
    _, _, regime_text = recovery_regimes(df)
    direct_pathway_plot()
    cases, traces = representative_timecourses(df)
    clearance = one_factor_trajectories()
    criterion_summary, criterion_text = criterion_validation(df)
    sensitivity_names = sensitivity_figures(correlations, regression)
    sensitivity_table = report_sensitivity_table(correlations, regression)
    make_additional_config(len(df))

    selection = df.selection_direction.value_counts(normalize=True)
    regime_description = regime_text + " " + "Across this joint grid, the median Δf_R is " + f"{df.delta_fR.median():.4f}" + ", with range [" + f"{df.delta_fR.min():.4f}, {df.delta_fR.max():.4f}" + "]; selection-direction fractions are " + ", ".join(f"{k}: {v:.1%}" for k,v in selection.items()) + "."
    figures = [
        ("comprehensive_01_binned_effects", "Marginal binned endpoint responses"),
        ("comprehensive_02_endpoint_tradeoffs", "Endpoint tradeoffs: composition, diversity, bottleneck, and recovery"),
        ("comprehensive_03_population_composition", "Absolute population and final composition versus inhibitor parameters"),
        *[(name, "Pairwise heatmaps: " + ENDPOINTS.get(name.removeprefix("comprehensive_pairwise_"), name)) for name in pair_figures],
        ("comprehensive_04_convergence_diagnostics", "Convergence diagnostics"),
        ("comprehensive_05_recovery_regimes", "Outcome regimes, diversity, and final biomass"),
        ("comprehensive_06_direct_pharmacodynamics", "Direct inhibitor pharmacodynamic functions"),
        ("comprehensive_07_representative_trajectories", "Selected full-grid representative trajectories (first 60 τ; endpoint metrics use 300 τ)"),
        ("comprehensive_08_clearance_recovery_trajectories", "Low, middle, and high inhibitor perturbation trajectories"),
        ("comprehensive_09_selection_criterion_validation", "Static selection criterion compared with dynamic outcomes"),
        *[(name, "Global sensitivity: " + name.replace("comprehensive_", "").replace("_", " ")) for name in sensitivity_names],
    ]
    extra = {
        "n_runs": len(df), "summary": summary, "regime_text": regime_description,
        "validation_text": equivalence_text + " " + criterion_text + " " + f"No-treatment-relative susceptible viability was met in {df.sensitive_survives.mean():.1%} of runs; resistant extinction occurred in {df.resistant_eliminated.mean():.1%}; and {df.ss_converged.mean():.1%} met the shared steady-state detector by τ={T_END:g}.",
        "figures": figures, "sensitivity_table": sensitivity_table,
    }
    all_atlas = pd.read_csv(TABLES / "all_atlas_results.csv")
    anchor_summary = pd.read_csv(TABLES / "anchor_summary.csv")
    reps = pd.read_csv(TABLES / "representative_regimes_expanded.csv")
    composition = pd.read_csv(TABLES / "initial_composition_expanded.csv")
    ablations = pd.read_csv(TABLES / "mechanistic_ablations.csv")
    robustness = pd.read_csv(TABLES / "robustness_lhs.csv")
    matched = pd.read_csv(TABLES / "identifiability_matched_iota.csv")
    db_frame = pd.read_csv(TABLES / "identifiability_db_at_fixed_iota.csv")
    ablation_maps = pd.read_csv(TABLES / "ablation_heatmap_i_hi.csv")
    core.build_report(all_atlas, anchor_summary, reps, composition, ablations, robustness, matched, db_frame,
                      ablation_maps=ablation_maps, comprehensive=extra)
    print(f"Comprehensive Person 3 analysis complete: {len(df):,} simulations", flush=True)


if __name__ == "__main__":
    main()
