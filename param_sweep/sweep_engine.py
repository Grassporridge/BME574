"""
sweep_engine.py

Sweep infrastructure built on top of ma2024_model.py for the alpha x beta_min x c
parameter sweep (growth burden, private benefit, inhibitor efficacy -- the three
parameters behind the paper's selection criterion, Eqs. 1-2).

WHY A SEPARATE ENGINE, NOT JUST ma2024_model.simulate()
--------------------------------------------------------
ma2024_model.py integrates the state (n_s, n_r, s, a, b) directly. That is fine
for a single illustrative run, but across a wide sweep of alpha/beta_min/c some
corners of the grid drive n_s or n_r down by 4-6 orders of magnitude before they
recover (a real "bottleneck"). Integrating the raw densities through that squeeze
loses precision exactly where the sweep's endpoints (peak loss, extinction flags)
matter most.

This module integrates ln(n_s) and ln(n_r) instead of n_s and n_r themselves.
The right-hand side is algebraically identical (see the docstring on
`rhs_logspace` below) -- it is the same five equations from ma2024_model.py,
just carried in a coordinate system that stays well-conditioned as a population
gets small. `verify_equivalence()` at the bottom checks this reformulation
against ma2024_model.simulate() point by point before it is used for anything.

ENDPOINTS COMPUTED PER PARAMETER COMBINATION
---------------------------------------------
  - f_R_final          : resistant fraction n_r/(n_s+n_r) once the system has
                          reached (numerical) steady state
  - resistant_eliminated : bool, n_r at steady state below EXTINCTION_THRESHOLD
  - sensitive_survives  : bool, n_s at steady state above a "viable" threshold,
                          defined as a fraction of the drug-free control's own
                          steady-state n_s (see `no_drug_reference()`)
  - t_recovery_S        : first time n_s climbs back above 90% of its own
                          eventual steady-state value (NaN if it never gets
                          there within the horizon)
  - t_steady_state       : first time the whole system's relative rate of
                          change stays below TOL_SS for the remainder of the
                          run (right-censored at t_span[1] if never reached)
  - peak_loss_S, peak_loss_total : deepest fractional drop from the starting
                          density for S alone and for the whole community,
                          i.e. how close the run passed to a bottleneck
  - diversity_final      : Shannon diversity H = -sum(p_i ln p_i) over
                          (n_s, n_r) at steady state
  - diversity_min        : the lowest Shannon diversity reached during the
                          run (how close the community came to a single-taxon
                          state, even if it recovered)

Two treatment scenarios are swept:
  "treatment" -- a0=2.15, i=10, the single dose of the paper's Fig. 2C. This is
                 where the requested endpoints (diversity, time-to-SS, f_R,
                 extinction, survival, recovery time, peak loss) are computed.
  "criterion" -- a0=100, i=10, the highest antibiotic and inhibitor doses in
                 the paper's dose range (SI: 1 < a0 < 100, i = 0.1-10), i.e. the
                 closest the paper's own dose grid gets to the saturating limit
                 its simplified criterion (Eq. 2) assumes. Used to validate that
                 criterion across a full structured alpha x beta_min x c grid
                 (the paper's Fig. 3C used 10,000 randomized strains).

PARAMETERS: every model parameter and initial condition is the Supplementary
Information base value (ma2024_model.DEFAULT_PARAMS / default_y0). alpha,
beta_min and c are swept over the SI's randomized ranges (alpha 0.75-1,
beta_min 0-1, c 0-1); gamma, h_a, h_i, kappa_b, d_a, d_b, xi, phi_max are held
at their SI base values (they are other group members' sweep variables).
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import replace
from typing import NamedTuple

from ma2024_model import (
    Params, DEFAULT_PARAMS, default_y0, hill, growth_rate, inhibitor_effect, beta_of, phi_of,
    simulate as simulate_raw,
)

# ---------------------------------------------------------------------------
# Thresholds (documented explicitly, used consistently across every endpoint)
# ---------------------------------------------------------------------------

EXTINCTION_THRESHOLD = 1e-4   # absolute density; below this, a population counts as
                              # practically extinct. The model is a deterministic
                              # ODE, so n_r can approach but never exactly reach 0
                              # in finite time (see module docstring in sweep
                              # notebooks) -- this is the operational cutoff.
VIABLE_FRACTION = 0.5         # "S survives" means its steady-state density is at
                              # least this fraction of the drug-free control's own
                              # steady-state density.
RECOVERY_FRACTION = 0.9       # "recovered" means back above this fraction of its
                              # own eventual steady-state value.
TOL_SS = 1e-3                 # tolerance on both populations' per-capita growth
                              # rates (|d ln n_s/dt|, |d ln n_r/dt|, units 1/time)
                              # used to detect steady state.
SS_HOLD_FRACTION = 0.1        # the relative rate must stay below TOL_SS for the
                              # last SS_HOLD_FRACTION of the run to be accepted
                              # as converged, not just a single quiet instant.


# ---------------------------------------------------------------------------
# Log-space reformulation
# ---------------------------------------------------------------------------

def rhs_logspace(tau, y, p: Params, i_dose: float):
    """Same five equations as ma2024_model.rhs, carried as
    (ln n_s, ln n_r, s, a, b) instead of (n_s, n_r, s, a, b).

    d(ln n_s)/dt = (1/n_s) dn_s/dt = g - l                [same as ma2024_model.rhs]
    d(ln n_r)/dt = (1/n_r) dn_r/dt = alpha*g - beta*l      [same]
    ds/dt, da/dt, db/dt : identical to ma2024_model.rhs, with n_s = exp(ln_n_s)
    and n_r = exp(ln_n_r) substituted in wherever they appear.
    """
    ln_ns, ln_nr, s, a, b = y
    ns = np.exp(ln_ns)
    nr = np.exp(ln_nr)
    g = growth_rate(s)
    l = p.gamma * hill(a, p.h_a) * g
    iota = inhibitor_effect(i_dose, p)
    beta = beta_of(iota, p)
    phi = phi_of(iota, p)

    dln_ns = g - l
    dln_nr = p.alpha * g - beta * l
    ds = (p.xi * l - g) * ns + (p.xi * beta * l - p.alpha * g) * nr
    da = -p.kappa_b * b * a - phi * nr * a - p.d_a * a
    db = beta * l * nr - p.d_b * iota * b

    return np.array([dln_ns, dln_nr, ds, da, db])


def simulate_logspace(
    p: Params,
    a0: float,
    i_dose: float,
    y0_lin: np.ndarray | None = None,
    t_span: tuple[float, float] = (0.0, 150.0),
    n_points: int = 300,
):
    """Integrate the log-space system. y0_lin is the initial condition in the
    ORIGINAL (linear) coordinates [n_s0, n_r0, s0, a0, b0]; a0 in y0_lin is
    overwritten by the `a0` argument, matching ma2024_model.simulate()'s
    convention.

    Returns a solve_ivp result whose .y rows are still (ln n_s, ln n_r, s, a, b)
    -- use `.ns`, `.nr` attributes (attached below) for the linear-space values.
    """
    if y0_lin is None:
        y0_lin = default_y0(a0)   # SI initial conditions
    else:
        y0_lin = np.array(y0_lin, dtype=float).copy()
        y0_lin[3] = a0

    y0 = y0_lin.copy()
    y0[0] = np.log(y0_lin[0])
    y0[1] = np.log(y0_lin[1])

    t_eval = np.linspace(*t_span, n_points)
    sol = solve_ivp(
        rhs_logspace, t_span, y0, args=(p, i_dose), t_eval=t_eval,
        method="LSODA", rtol=1e-8, atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")

    sol.ns = np.exp(sol.y[0])
    sol.nr = np.exp(sol.y[1])
    sol.s = sol.y[2]
    sol.a = sol.y[3]
    sol.b = sol.y[4]
    return sol


# ---------------------------------------------------------------------------
# Equivalence check against ma2024_model.py (must be run before trusting this
# module for anything -- see notebook 1)
# ---------------------------------------------------------------------------

def verify_equivalence(n_checks: int = 12, seed: int = 0, rtol: float = 1e-4,
                        atol: float = 1e-6,
                        t_span=(0.0, 40.0), a0=2.15, i_dose=10.0) -> dict:
    """Compares simulate_logspace() against ma2024_model.simulate() (the
    already-delivered, reviewed implementation) at randomized alpha/beta_min/c
    points. Returns a dict of per-check max |raw-log| combined absolute/relative
    differences in n_s(t), n_r(t), s(t), a(t), b(t); raises AssertionError if
    any exceeds its tolerance.

    Uses an absolute-plus-relative closeness test (|raw-log| <= atol +
    rtol*|raw|), the same pattern as numpy.isclose, rather than a pure
    relative error. A pure relative test blows up harmlessly once a state
    variable (s in particular, once nutrient is depleted) decays past both
    solvers' own atol=1e-10 floor: two numbers like 2.5e-14 and -1.2e-10 are
    both "zero" for any practical purpose, but their relative difference is
    formally huge. This only re-derives the SAME equations in a different
    coordinate system, so any real bug (a sign error, a mis-transcribed term)
    would still show up here as a large ABSOLUTE difference at a state where
    the raw trajectory is not itself near the solver's noise floor.
    """
    rng = np.random.default_rng(seed)
    alphas = rng.uniform(0.75, 1.0, n_checks)
    beta_mins = rng.uniform(0.0, 1.0, n_checks)
    cs = rng.uniform(0.0, 1.0, n_checks)

    max_diffs = {"ns": 0.0, "nr": 0.0, "s": 0.0, "a": 0.0, "b": 0.0}
    for alpha, beta_min, c in zip(alphas, beta_mins, cs):
        p = replace(DEFAULT_PARAMS, alpha=alpha, beta_min=beta_min, c=c)
        sol_raw = simulate_raw(p, a0=a0, i_dose=i_dose, t_span=t_span, n_points=150)
        sol_log = simulate_logspace(p, a0=a0, i_dose=i_dose, t_span=t_span, n_points=150)

        for key, raw_row, log_val in [
            ("ns", sol_raw.y[0], sol_log.ns),
            ("nr", sol_raw.y[1], sol_log.nr),
            ("s", sol_raw.y[2], sol_log.s),
            ("a", sol_raw.y[3], sol_log.a),
            ("b", sol_raw.y[4], sol_log.b),
        ]:
            allowed = atol + rtol * np.abs(raw_row)
            excess = np.abs(raw_row - log_val) - allowed
            # Report how far past its own tolerance the worst point is, as a
            # fraction of that tolerance (0 or negative = within tolerance).
            worst_ratio = float(np.max(excess / allowed))
            max_diffs[key] = max(max_diffs[key], worst_ratio)

    for key, val in max_diffs.items():
        assert val <= 0.0, (
            f"Log-space reformulation disagrees with ma2024_model.py on '{key}': "
            f"worst point exceeds its atol+rtol*|raw| bound by a factor of {val:.2f}"
        )
    return max_diffs


# ---------------------------------------------------------------------------
# Reference (drug-free) control, used to define "viable" for S
# ---------------------------------------------------------------------------

def no_drug_reference(t_span=(0.0, 150.0)) -> float:
    """Steady-state n_s of the BASE strain (SI base parameters, alpha = 0.95)
    with no antibiotic at all. With a = 0 the lysis term is zero, so beta_min, c
    and gamma play no role; alpha still does (R competes with S for nutrient),
    so this is the base mixture's drug-free S density, used as ONE fixed
    yardstick for "S is viable" across the whole grid.
    """
    sol = simulate_logspace(DEFAULT_PARAMS, a0=0.0, i_dose=0.0, t_span=t_span, n_points=50)
    return float(sol.ns[-1])


# ---------------------------------------------------------------------------
# Steady-state detection
# ---------------------------------------------------------------------------

def _population_growth_rates(sol, p: Params, i_dose: float) -> np.ndarray:
    """max(|d ln n_s/dt|, |d ln n_r/dt|) at each saved time point -- i.e. the
    larger of the two populations' own per-capita growth rates, evaluated
    from the analytic RHS.

    These two quantities are already unitless (1/time) per-capita rates, so
    no further normalization is needed, and -- unlike a plain relative rate
    on s, a or b -- they correctly go to zero once the population has
    actually settled. (An earlier version of this function also required
    |da/dt|/a and |db/dt|/b to be small; that criterion never triggers,
    because once resistant cells are essentially gone the antibiotic simply
    decays as a*exp(-d_a*t), which has a CONSTANT relative rate d_a forever
    by definition of exponential decay -- so it looks perpetually
    "unconverged" even though the population dynamics have long since
    stopped changing. What matters for the sweep's endpoints is whether
    n_s and n_r have stopped moving, not whether every trace amount of
    residual antibiotic has fully vanished.)"""
    n = len(sol.t)
    rate = np.zeros(n)
    for k in range(n):
        y = np.array([sol.y[0, k], sol.y[1, k], sol.s[k], sol.a[k], sol.b[k]])
        dydt = rhs_logspace(sol.t[k], y, p, i_dose)
        rate[k] = max(abs(dydt[0]), abs(dydt[1]))
    return rate


def time_to_steady_state(sol, p: Params, i_dose: float, tol: float = TOL_SS,
                          hold_fraction: float = SS_HOLD_FRACTION) -> tuple[float, bool]:
    """First time t* such that both populations' per-capita growth rates
    (|d ln n_s/dt|, |d ln n_r/dt|) stay below `tol` for the rest of the run
    (i.e. from t* to the end covers at least hold_fraction of the total
    horizon and never exceeds tol). Returns (t*, converged) -- converged=
    False means the system had not settled by the end of the horizon, and
    t* is reported as the horizon's end (right-censored)."""
    rate = _population_growth_rates(sol, p, i_dose)
    t = sol.t
    below = rate < tol
    # Walk backward from the end; find the earliest index from which "below"
    # holds all the way to the end.
    n = len(t)
    idx = n - 1
    while idx > 0 and below[idx - 1]:
        idx -= 1
    t_star = t[idx]
    span = t[-1] - t[0]
    converged = below[-1] and (t[-1] - t_star) >= hold_fraction * span
    if not converged:
        return float(t[-1]), False
    return float(t_star), True


# ---------------------------------------------------------------------------
# Endpoint extraction for a single (alpha, beta_min, c) combination
# ---------------------------------------------------------------------------

def shannon_diversity_arr(ns: np.ndarray, nr: np.ndarray) -> np.ndarray:
    total = ns + nr
    with np.errstate(divide="ignore", invalid="ignore"):
        p_s = np.where(total > 0, ns / total, 0.0)
        p_r = np.where(total > 0, nr / total, 0.0)
        h_s = np.where(p_s > 0, -p_s * np.log(p_s), 0.0)
        h_r = np.where(p_r > 0, -p_r * np.log(p_r), 0.0)
    return h_s + h_r


class Endpoints(NamedTuple):
    alpha: float
    beta_min: float
    c: float
    f_R_final: float
    resistant_eliminated: bool
    sensitive_survives: bool
    t_recovery_S: float
    t_steady_state: float
    ss_converged: bool
    peak_loss_S: float
    peak_loss_total: float
    diversity_final: float
    diversity_min: float
    ns_final: float
    nr_final: float


def compute_endpoints(
    alpha: float, beta_min: float, c: float,
    a0: float, i_dose: float,
    ns0_no_drug: float,
    t_span: tuple[float, float] = (0.0, 150.0),
    n_points: int = 300,
) -> Endpoints:
    p = replace(DEFAULT_PARAMS, alpha=alpha, beta_min=beta_min, c=c)
    sol = simulate_logspace(p, a0=a0, i_dose=i_dose, t_span=t_span, n_points=n_points)
    ns, nr, t = sol.ns, sol.nr, sol.t
    total = ns + nr

    ns0, nr0 = ns[0], nr[0]
    ns_final, nr_final = ns[-1], nr[-1]

    f_R_final = nr_final / (ns_final + nr_final) if (ns_final + nr_final) > 0 else np.nan
    resistant_eliminated = bool(nr_final < EXTINCTION_THRESHOLD)
    sensitive_survives = bool(ns_final > VIABLE_FRACTION * ns0_no_drug)

    # Time-to-recovery for S: first time n_s climbs back above
    # RECOVERY_FRACTION * its own eventual (steady-state) value, only counted
    # if it actually dipped below that line at some point (otherwise there
    # was no crash to recover from, and t_recovery_S = 0 by convention).
    recovery_line = RECOVERY_FRACTION * ns_final
    dipped = np.any(ns < recovery_line)
    if not dipped:
        t_recovery_S = 0.0
    else:
        above_again = np.where(ns >= recovery_line)[0]
        above_again = above_again[above_again > np.argmax(ns < recovery_line)] \
            if len(above_again) else above_again
        # first index, after the dip begins, where ns is back above the line
        dip_start = np.argmax(ns < recovery_line)
        after_dip = np.where(ns[dip_start:] >= recovery_line)[0]
        t_recovery_S = float(t[dip_start + after_dip[0]]) if len(after_dip) else np.nan

    t_ss, ss_converged = time_to_steady_state(sol, p, i_dose)

    peak_loss_S = float((ns0 - np.min(ns)) / ns0) if ns0 > 0 else np.nan
    peak_loss_total = float(((ns0 + nr0) - np.min(total)) / (ns0 + nr0))

    div = shannon_diversity_arr(ns, nr)
    diversity_final = float(div[-1])
    diversity_min = float(np.min(div))

    return Endpoints(
        alpha=alpha, beta_min=beta_min, c=c,
        f_R_final=float(f_R_final),
        resistant_eliminated=resistant_eliminated,
        sensitive_survives=sensitive_survives,
        t_recovery_S=t_recovery_S,
        t_steady_state=t_ss,
        ss_converged=ss_converged,
        peak_loss_S=peak_loss_S,
        peak_loss_total=peak_loss_total,
        diversity_final=diversity_final,
        diversity_min=diversity_min,
        ns_final=float(ns_final),
        nr_final=float(nr_final),
    )


# ---------------------------------------------------------------------------
# Criterion check for a single combination (saturating-dose scenario)
# ---------------------------------------------------------------------------

def resistant_advantage_logspace(sol, p: Params, i_dose: float) -> float:
    """Net change in ln(n_r/n_s) over the run. Since
    d ln(n_r/n_s)/dt = (alpha*g - beta*l) - (g - l) = -(1-alpha)*g + (1-beta)*l,
    this is the integral of the resistant cells' relative growth advantage.
    The log-space solution carries ln n_s and ln n_r directly, so the change is
    read off exactly from the solution's end points (no quadrature error, no
    division of two underflowed densities). Positive => resistant fraction grew."""
    return float((sol.y[1, -1] - sol.y[0, -1]) - (sol.y[1, 0] - sol.y[0, 0]))


def criterion_check(alpha: float, beta_min: float, c: float,
                     a0: float = 100.0, i_dose: float = 10.0,
                     t_span=(0.0, 600.0), n_points: int = 50) -> dict:
    p = replace(DEFAULT_PARAMS, alpha=alpha, beta_min=beta_min, c=c)
    sol = simulate_logspace(p, a0=a0, i_dose=i_dose, t_span=t_span, n_points=n_points)
    advantage = resistant_advantage_logspace(sol, p, i_dose)
    x = (1.0 - alpha) / p.gamma
    y = (1.0 - c) * (1.0 - beta_min)
    return {
        "alpha": alpha, "beta_min": beta_min, "c": c,
        "x_val": x, "y_val": y,
        "criterion_predicts_R": bool(y > x),
        "simulated_R_favored": bool(advantage > 0),
        "advantage": advantage,
    }
