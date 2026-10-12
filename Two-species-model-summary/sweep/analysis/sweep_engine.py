"""
sweep_engine.py -- one simulation + endpoint engine shared by every sweep.

Every parameter group (alpha/beta_min/c, gamma/h_a/a0, i/h_i/d_b,
kappa_b/d_a/phi_max/xi) calls the same `compute_endpoints()`, so all result
tables have the same columns, the same definitions and the same numerical
settings.

NUMERICS
--------
* Same equations as src/model/ma2024_model.py (SI Eqs. S1-S10), integrated in
  (ln n_s, ln n_r, s, a, b) so populations stay accurate through deep crashes.
  `verify_equivalence()` checks this against the linear-coordinate model.
* Adaptive horizon. A run is integrated to t = 150 (sampled every 0.5). If the
  populations have not settled it is continued to t = 600 (every 1.0) and then
  to t = 2400 (every 2.0). This replaces the per-notebook fixed horizons
  (150, 300, 400) used before.
* Exact minima. S stops declining at the moment the antibiotic falls to
  a*_S = (gamma - 1)^(-1/h_a) (where lysis equals growth), and R at
  a*_R = (alpha/(beta gamma - alpha))^(1/h_a). The solver locates those two
  moments exactly, so lowest densities and peak losses do not depend on the
  sampling step.
* Steady state: both populations' per-capita growth rates |d ln n/dt| stay
  below TOL_SS from some time t* to the end of the run, and that quiet stretch
  covers at least SS_HOLD_FRACTION of the horizon.

ENDPOINT DEFINITIONS (identical for every sweep)
------------------------------------------------
  ns_final, nr_final, total_final   densities at the end of the run
  f_R_final                          n_r/(n_s+n_r) at the end
  diversity_final, diversity_min     Shannon H over {S, R}: final / lowest
  ns_min, nr_min, total_min          lowest densities reached
  peak_loss_S, peak_loss_R, peak_loss_total
                                     largest fractional drop below the START
  t_S_min                            time of the lowest S
  t_regain_S                         time after the dose at which S is back at
                                     its starting density (0 if S never dipped;
                                     NaN if it never comes back). This is the
                                     recovery time from the crash.
  t_90_final_S                       time for S to reach 90% of its final
                                     density (growth included). This is what
                                     earlier notebooks called "t_recovery_S".
  t_steady_state, ss_converged       see above (NaN if not converged)
  t_a_below_1                        time for the antibiotic to fall below a = 1
  delta_log_ratio                    change in ln(n_r/n_s); > 0 means R gained
  sensitive_survives                 ns_final >= VIABLE_FRACTION x drug-free S
  resistant_eliminated_abs           legacy flag: nr_final < 1e-4
  outcome                            end-state class of the plain model

EXTINCTION FLOOR (the "_floor" columns)
---------------------------------------
The plain ODE lets a population fall to 1e-24 and regrow. Real populations
cannot: below one cell they are gone. `EXTINCTION_FLOOR` is one cell in a
culture that saturates at ~1e9 cells (1e-9 of the carrying capacity). When a
run's lowest density crosses the floor it is re-run with that population
removed for good at the moment it crosses (its consequences included: e.g. no
more antibiotic degradation by R). Columns: S_extinct_floor, R_extinct_floor,
t_extinct_S, t_extinct_R, ns_final_floor, nr_final_floor, f_R_final_floor,
outcome_floor, goal_met_floor (R gone for good and S viable at the end). Runs
that never cross the floor are identical to the plain model.
"""

from __future__ import annotations
import math
from dataclasses import replace
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp

from ..model.ma2024_model import (
    Params, DEFAULT_PARAMS, PARAM_NAMES, SI_INITIAL_CONDITIONS, SI_RANDOMIZED_RANGES,
    GROUP_CHOSEN_RANGES, SI_DOSE_RANGES, BASE_A0, BASE_I, default_y0,
    simulate as simulate_linear,
)

# ---------------------------------------------------------------------------
# Settings shared by every sweep
# ---------------------------------------------------------------------------
CARRYING_CAPACITY = (SI_INITIAL_CONDITIONS["ns0"] + SI_INITIAL_CONDITIONS["nr0"]
                     + SI_INITIAL_CONDITIONS["s0"])          # 4.4: drug-free final total
CELLS_AT_SATURATION = 1e9                                     # assumption (cells per mL)
EXTINCTION_FLOOR = CARRYING_CAPACITY / CELLS_AT_SATURATION    # density of one cell
F_R_EXTINCT = 1e-9          # R "extinct" at the end if f_R below this
F_R_SUPPRESSED = 1e-6       # R "suppressed" if f_R below this
LEGACY_ABS_THRESHOLD = 1e-4  # the earlier notebooks' n_r < 1e-4 flag
VIABLE_FRACTION = 0.5
RECOVERY_FRACTION = 0.9
TOL_SS = 1e-3
SS_HOLD_FRACTION = 0.1
HORIZON_STAGES = ((150.0, 0.5), (600.0, 1.0), (2400.0, 2.0))   # (end time, sampling step)
RTOL, ATOL = 1e-8, 1e-10

ABLATIONS = ("full", "private_only", "degradation_only", "turnover_only", "none")


def pathway(p: Params, i_dose: float, ablation: str = "full"):
    """(beta, phi, d_b*iota) for the dose. The inhibitor acts on three things:
    R's private protection (beta), antibiotic degradation by living R (phi) and
    inactivation of free beta-lactamase (d_b*iota). `ablation` lets the inhibitor
    act on only one of them (the others see no inhibitor), to separate the three.
    """
    i_dose = max(float(i_dose), 0.0)
    ih = i_dose ** p.h_i
    iota = ih / (1.0 + ih)
    on = {"full": (1, 1, 1), "private_only": (1, 0, 0), "degradation_only": (0, 1, 0),
          "turnover_only": (0, 0, 1), "none": (0, 0, 0)}[ablation]
    beta = p.beta_min + p.c * (1.0 - p.beta_min) * iota * on[0]
    phi = p.phi_max * (1.0 - p.c * iota * on[1])
    return beta, phi, p.d_b * iota * on[2]


class Trajectory:
    """Result of `integrate`. ns/nr are 0 after a population is removed."""
    __slots__ = ("t", "ln_ns", "ln_nr", "ns", "nr", "s", "a", "b", "converged", "t_ss",
                 "t_extinct_S", "t_extinct_R", "p", "a0", "i_dose", "beta", "phi", "y",
                 "turning_points")


def _steady_state(t, ln_ns, ln_nr, s, a, p: Params, beta: float):
    """(t*, converged): see module docstring."""
    g = s / (1.0 + s)
    ah = np.power(np.maximum(a, 0.0), p.h_a)
    l = p.gamma * ah / (1.0 + ah) * g
    rate_s = np.where(np.isfinite(ln_ns), np.abs(g - l), 0.0)
    rate_r = np.where(np.isfinite(ln_nr), np.abs(p.alpha * g - beta * l), 0.0)
    below = np.maximum(rate_s, rate_r) < TOL_SS
    idx = len(t) - 1
    if not below[idx]:
        return float("nan"), False
    while idx > 0 and below[idx - 1]:
        idx -= 1
    converged = (t[-1] - t[idx]) >= SS_HOLD_FRACTION * t[-1]
    return (float(t[idx]) if converged else float("nan")), bool(converged)


def integrate(p: Params, a0: float, i_dose: float, y0_lin=None,
              stages=HORIZON_STAGES, adaptive: bool = True, floor: float | None = None,
              ablation: str = "full", method: str = "LSODA",
              rtol: float = RTOL, atol: float = ATOL) -> Trajectory:
    """Integrate the model. With `adaptive`, stop after the first stage in which
    the populations have settled. With `floor`, a population that falls below
    `floor` is removed permanently at that moment."""
    beta, phi, dbi = pathway(p, i_dose, ablation)
    alpha, gamma, h_a, xi, kb, da = p.alpha, p.gamma, p.h_a, p.xi, p.kappa_b, p.d_a
    alive = [1.0, 1.0]

    def fun(t, y):
        s, a, b = y[2], y[3], y[4]
        g = s / (1.0 + s)
        if a > 0.0:
            ah = a ** h_a
            l = gamma * ah / (1.0 + ah) * g
        else:
            l = 0.0
        ns = math.exp(y[0]) * alive[0]
        nr = math.exp(y[1]) * alive[1]
        return [(g - l) * alive[0],
                (alpha * g - beta * l) * alive[1],
                (xi * l - g) * ns + (xi * beta * l - alpha * g) * nr,
                -kb * b * a - phi * nr * a - da * a,
                beta * l * nr - dbi * b]

    y0 = default_y0(a0) if y0_lin is None else np.array(y0_lin, dtype=float).copy()
    y0[3] = a0
    y = np.array([math.log(y0[0]), math.log(y0[1]), y0[2], y0[3], y0[4]])
    ln_floor = math.log(floor) if floor else None
    t_ext = [float("nan"), float("nan")]

    # Antibiotic levels at which each population stops declining (lysis = growth).
    a_star = [(gamma - 1.0) ** (-1.0 / h_a) if gamma > 1.0 else None,
              (alpha / (beta * gamma - alpha)) ** (1.0 / h_a) if beta * gamma > alpha else None]
    turning = []                                   # (time, ln n_s, ln n_r) at those moments

    T = [np.array([0.0])]
    Y = [y.reshape(5, 1).copy()]
    t_cur, converged, t_ss = 0.0, False, float("nan")

    def stored(block):
        block = block.copy()
        for k in (0, 1):
            if alive[k] == 0.0:
                block[k] = -np.inf
        return block

    for t_end, dt in stages:
        n = int(round((t_end - t_cur) / dt))
        grid = t_cur + dt * np.arange(1, n + 1)
        grid[-1] = t_end
        while t_cur < t_end and (alive[0] or alive[1]):
            t_eval = grid[grid > t_cur + 1e-12]
            events = []
            for k in (0, 1):
                if alive[k] and a_star[k] is not None and y[3] > a_star[k] * (1.0 + 1e-9):
                    ev = (lambda t, y, k=k: y[3] - a_star[k])
                    ev.terminal, ev.direction, ev.pop = False, -1, k
                    events.append(ev)
                if alive[k] and ln_floor is not None:
                    ev = (lambda t, y, k=k: y[k] - ln_floor)
                    ev.terminal, ev.direction, ev.pop = True, -1, k
                    events.append(ev)
            try:
                sol = solve_ivp(fun, (t_cur, t_end), y, t_eval=t_eval, events=events or None,
                                method=method, rtol=rtol, atol=atol)
            except ValueError:
                # Root-finding for a turning-point event can fail when the event sits
                # exactly on a step boundary. Those events only refine the minima, so
                # retry with the extinction events alone.
                events = [ev for ev in events if ev.terminal]
                sol = solve_ivp(fun, (t_cur, t_end), y, t_eval=t_eval, events=events or None,
                                method=method, rtol=rtol, atol=atol)
            if sol.status < 0:
                raise RuntimeError(f"Integration failed: {sol.message}")
            if len(sol.t):
                T.append(sol.t)
                Y.append(stored(sol.y))
            for ev, te, ye in zip(events, sol.t_events or [], sol.y_events or []):
                if not ev.terminal and len(te):
                    turning.append((float(te[0]),
                                    float(ye[0][0]) if alive[0] else -np.inf,
                                    float(ye[0][1]) if alive[1] else -np.inf))
            if sol.status == 1:                         # a population hit the floor
                for ev, te, ye in zip(events, sol.t_events, sol.y_events):
                    if ev.terminal and len(te):
                        t_cur, y = float(te[0]), np.array(ye[0])
                        break
                for k in (0, 1):                        # both can cross at the same moment
                    if alive[k] and y[k] <= ln_floor + 1e-9:
                        alive[k] = 0.0
                        t_ext[k] = t_cur
                T.append(np.array([t_cur]))
                Y.append(stored(y.reshape(5, 1)))
            else:
                t_cur, y = t_end, sol.y[:, -1].copy()
        if not (alive[0] or alive[1]):
            converged, t_ss = True, t_cur
            break
        t_all, y_all = np.concatenate(T), np.concatenate(Y, axis=1)
        t_ss, converged = _steady_state(t_all, y_all[0], y_all[1], y_all[2], y_all[3], p, beta)
        if adaptive and converged:
            break

    tr = Trajectory()
    tr.t, tr.y = np.concatenate(T), np.concatenate(Y, axis=1)
    tr.ln_ns, tr.ln_nr, tr.s, tr.a, tr.b = tr.y
    with np.errstate(over="ignore"):
        tr.ns, tr.nr = np.exp(tr.ln_ns), np.exp(tr.ln_nr)
    tr.converged, tr.t_ss = converged, t_ss
    tr.turning_points = turning
    tr.t_extinct_S, tr.t_extinct_R = t_ext
    tr.p, tr.a0, tr.i_dose, tr.beta, tr.phi = p, float(a0), float(i_dose), beta, phi
    return tr


def simulate(p: Params = DEFAULT_PARAMS, a0: float = BASE_A0, i_dose: float = BASE_I,
             t_end: float = 60.0, dt: float = 0.1, **kwargs) -> Trajectory:
    """Fixed-horizon run for plotting time courses."""
    return integrate(p, a0, i_dose, stages=((float(t_end), float(dt)),), adaptive=False, **kwargs)


def shannon(ns, nr):
    """Shannon diversity over {S, R}; NaN where both are zero."""
    ns, nr = np.asarray(ns, float), np.asarray(nr, float)
    total = ns + nr
    with np.errstate(divide="ignore", invalid="ignore"):
        ps, pr = ns / total, nr / total
        h = -(np.where(ps > 0, ps * np.log(ps), 0.0) + np.where(pr > 0, pr * np.log(pr), 0.0))
    return np.where(total > 0, h, np.nan)


@lru_cache(maxsize=None)
def no_drug_reference() -> float:
    """Final S of the base mixture with no antibiotic and no inhibitor: the one
    fixed yardstick for 'S is viable' in every sweep."""
    return float(integrate(DEFAULT_PARAMS, 0.0, 0.0).ns[-1])


def _first_time(t, mask):
    idx = np.flatnonzero(mask)
    return float(t[idx[0]]) if len(idx) else float("nan")


def classify(f_R, total_final, converged, sensitive_survives) -> str:
    """End-state class of the plain model."""
    if total_final < EXTINCTION_FLOOR:
        return "community collapsed"
    if not converged:
        return "not converged"
    if f_R < F_R_EXTINCT:
        return "R cleared"
    if f_R < F_R_SUPPRESSED:
        return "R suppressed"
    return "S viable, R persists" if sensitive_survives else "S not viable, R persists"


OUTCOME_ORDER = ["R cleared", "R suppressed", "S viable, R persists", "S not viable, R persists",
                 "community collapsed", "not converged"]
OUTCOME_FLOOR_ORDER = ["R eradicated, S survives", "R eradicated, S not viable",
                       "community eradicated", "S eradicated, R survives"] + OUTCOME_ORDER


def endpoints(tr: Trajectory, ns_ref: float | None = None) -> dict:
    """All endpoints of one plain-model trajectory (see module docstring)."""
    ns_ref = no_drug_reference() if ns_ref is None else ns_ref
    t, ns, nr = tr.t, tr.ns, tr.nr
    total = ns + nr
    ns0, nr0 = ns[0], nr[0]
    nsf, nrf, totf = float(ns[-1]), float(nr[-1]), float(total[-1])
    collapsed = totf < EXTINCTION_FLOOR
    f_R = nrf / totf if totf > 0 else float("nan")
    H = shannon(ns, nr)
    # Lowest densities: sampled values plus the exactly located turning points.
    i_min = int(np.argmin(ns))
    ns_min, t_S_min = float(ns[i_min]), float(t[i_min])
    nr_min, total_min = float(nr.min()), float(total.min())
    for tt, ln_s, ln_r in tr.turning_points:
        es, er = math.exp(ln_s), math.exp(ln_r)
        if es < ns_min:
            ns_min, t_S_min = es, tt
        nr_min = min(nr_min, er)
        total_min = min(total_min, es + er)
    peak_loss_S = float(max(0.0, 1.0 - ns_min / ns0))

    if peak_loss_S == 0.0:
        t_regain = 0.0
    else:
        after = t >= t_S_min
        t_regain = _first_time(t[after], ns[after] >= ns0)
    t_90 = float("nan") if collapsed else _first_time(t, ns >= RECOVERY_FRACTION * nsf)
    survives = bool(nsf >= VIABLE_FRACTION * ns_ref)

    return {
        "ns_final": nsf, "nr_final": nrf, "total_final": totf,
        "f_R_final": float("nan") if collapsed else f_R,
        "diversity_final": float("nan") if collapsed else float(H[-1]),
        "diversity_min": float("nan") if collapsed else float(np.nanmin(H)),
        "ns_min": ns_min, "nr_min": nr_min, "total_min": total_min,
        "peak_loss_S": peak_loss_S,
        "peak_loss_R": float(max(0.0, 1.0 - nr_min / nr0)),
        "peak_loss_total": float(max(0.0, 1.0 - total_min / (ns0 + nr0))),
        "t_S_min": t_S_min,
        "t_regain_S": t_regain,
        "t_90_final_S": t_90,
        "t_steady_state": tr.t_ss, "ss_converged": tr.converged, "t_end": float(t[-1]),
        "t_a_below_1": 0.0 if tr.a0 <= 1.0 else _first_time(t, tr.a <= 1.0),
        "delta_log_ratio": float((tr.ln_nr[-1] - tr.ln_ns[-1]) - (tr.ln_nr[0] - tr.ln_ns[0])),
        "sensitive_survives": survives,
        "resistant_eliminated_abs": bool(nrf < LEGACY_ABS_THRESHOLD),
        "outcome": classify(f_R, totf, tr.converged, survives),
    }


def compute_endpoints(overrides: dict | None = None, a0: float = BASE_A0, i_dose: float = BASE_I,
                      ns_ref: float | None = None, with_floor: bool = True,
                      ablation: str = "full") -> dict:
    """Run one parameter combination and return its endpoint row.
    `overrides` maps Params field names to values; everything else stays at the
    SI base value."""
    p = replace(DEFAULT_PARAMS, **(overrides or {}))
    ns_ref = no_drug_reference() if ns_ref is None else ns_ref
    tr = integrate(p, a0, i_dose, ablation=ablation)
    row = endpoints(tr, ns_ref)
    if not with_floor:
        return row

    crossed = row["ns_min"] < EXTINCTION_FLOOR or row["nr_min"] < EXTINCTION_FLOOR
    if crossed:
        trf = integrate(p, a0, i_dose, floor=EXTINCTION_FLOOR, ablation=ablation)
        s_ext, r_ext = not math.isnan(trf.t_extinct_S), not math.isnan(trf.t_extinct_R)
        nsf, nrf = float(trf.ns[-1]), float(trf.nr[-1])
        if s_ext and r_ext:
            out = "community eradicated"
        elif r_ext:
            out = ("R eradicated, S survives" if nsf >= VIABLE_FRACTION * ns_ref
                   else "R eradicated, S not viable")
        elif s_ext:
            out = "S eradicated, R survives"
        else:
            out = row["outcome"]
        row.update({
            "S_extinct_floor": s_ext, "R_extinct_floor": r_ext,
            "t_extinct_S": trf.t_extinct_S, "t_extinct_R": trf.t_extinct_R,
            "ns_final_floor": nsf, "nr_final_floor": nrf,
            "f_R_final_floor": nrf / (nsf + nrf) if (nsf + nrf) > 0 else float("nan"),
            "outcome_floor": out,
        })
    else:
        row.update({
            "S_extinct_floor": False, "R_extinct_floor": False,
            "t_extinct_S": float("nan"), "t_extinct_R": float("nan"),
            "ns_final_floor": row["ns_final"], "nr_final_floor": row["nr_final"],
            "f_R_final_floor": row["f_R_final"], "outcome_floor": row["outcome"],
        })
    row["goal_met_floor"] = row["outcome_floor"] == "R eradicated, S survives"
    return row


# ---------------------------------------------------------------------------
# Verification helper (used by notebook 00)
# ---------------------------------------------------------------------------
ALL_RANGES = {**SI_RANDOMIZED_RANGES, **GROUP_CHOSEN_RANGES}


def random_parameter_sets(n: int, seed: int = 0, names=None) -> list[dict]:
    """n parameter sets drawn uniformly from the sweep ranges (all parameters,
    or only `names`)."""
    rng = np.random.default_rng(seed)
    names = list(ALL_RANGES) if names is None else list(names)
    return [{k: float(rng.uniform(*ALL_RANGES[k])) for k in names} for _ in range(n)]


def verify_equivalence(n_checks: int = 40, seed: int = 0, rtol: float = 1e-4,
                       atol: float = 1e-6, t_end: float = 40.0) -> dict:
    """Compare the log-coordinate engine with the linear-coordinate reference
    model at random points over ALL parameters and both doses. Returns the worst
    |difference| / (atol + rtol*|reference|) per state variable (<= 1 passes)."""
    rng = np.random.default_rng(seed + 1)
    worst = {k: 0.0 for k in ("ns", "nr", "s", "a", "b")}
    for over in random_parameter_sets(n_checks, seed):
        a0 = float(np.exp(rng.uniform(*np.log(SI_DOSE_RANGES["a0"]))))
        i_dose = float(np.exp(rng.uniform(*np.log(SI_DOSE_RANGES["i"]))))
        p = replace(DEFAULT_PARAMS, **over)
        n_pts = int(t_end / 0.25) + 1
        ref = simulate_linear(p, a0, i_dose, t_span=(0.0, t_end), n_points=n_pts,
                              rtol=1e-10, atol=1e-12)
        tr = integrate(p, a0, i_dose, stages=((t_end, 0.25),), adaptive=False)
        for key, r, v in zip(worst, ref.y, (tr.ns, tr.nr, tr.s, tr.a, tr.b)):
            worst[key] = max(worst[key], float(np.max(np.abs(r - v) / (atol + rtol * np.abs(r)))))
    return worst
