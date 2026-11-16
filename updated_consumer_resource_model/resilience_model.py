"""
Pre-existing resistance and post-antibiotic community recovery
==============================================================

A minimal model of the barcoded E. coli community experiment: grow the community,
treat it with carbenicillin, dilute into fresh drug-free medium, let it recover, and
measure how much of the community came back.

Everything needed is in this one file:

    BASELINE          every parameter, with units
    run(**changes)    one matched experiment -> all readouts
    sweep(...)        one parameter over a list of values
    grid(...)         two parameters over a rectangle
    check()           sanity checks; run once per session, not on import

A notebook needs nothing but:

    from resilience_model import run, sweep, grid, BASELINE, check
    r = run(growth_mean=0.2)          # anything not named keeps its baseline value
    r.metrics                         # one row of readouts
    df = sweep("n_resistant", [0, 1, 2, 4, 8, 16])

See README.md for the assumptions, their justification, and the tested ranges.

Units: time in hours, biomass in OD, drug in ug/ml of carbenicillin.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

__all__ = ["BASELINE", "Params", "run", "sweep", "grid", "check", "Result",
           "TESTED_RANGES", "SAMPLE_HOURS"]

SAMPLE_HOURS = (4.0, 6.0, 24.0, 48.0)   # when the experiment samples


# --------------------------------------------------------------------------- parameters
@dataclass(frozen=True)
class Params:
    """Every parameter of the model. `BASELINE` is the default instance.

    Drug concentrations are in ug/ml of carbenicillin, so `dose` is the number the
    protocol specifies and `km` is TEM-1's measured Michaelis constant in the same unit.
    """
    # --- community -------------------------------------------------------------
    seed: int = 42                  # fixes every strain trait; change to test robustness
    n_strains: int = 57             # barcoded library size
    focal_index: int = 29           # k30, the strain that carries the resistance gene

    # --- growth ----------------------------------------------------------------
    growth_mean: float = 0.55       # 1/h, mean maximum growth rate across strains
    growth_cv: float = 0.08         # strain-to-strain spread in maximum growth rate
    food: float = 1.0               # OD-equivalents of nutrient in fresh medium
    half_saturation: float = 0.05   # OD-equivalents; Monod constant for nutrient uptake
    death_rate: float = 0.004       # 1/h, constant death rate (fed or starving)

    # --- resistance ------------------------------------------------------------
    n_resistant: int = 1            # how many strains carry the gene; 0 = sensitive community
    focal_multiplier: float = 1.0   # starting biomass of k30 relative to the others (1x or 10x)
    resistance_cost: float = 0.05   # fractional loss of maximum growth rate for a carrier
    carrier_lysis: float = 0.0      # carrier killing as a fraction of a susceptible cell's
                                    # (Ma et al.'s beta): 0 = immune, 1 = no private protection

    # --- antibiotic ------------------------------------------------------------
    dose: float = 50.0              # ug/ml carbenicillin added at the start of treatment
    gamma: float = 2.5              # max killing rate, as a multiple of the cell's current growth
    hill: float = 2.0               # Hill coefficient of the killing response
    susceptible_ic50: float = 5.0   # ug/ml, median IC50 of the susceptible strains
    ic50_cv: float = 0.70           # strain-to-strain spread in IC50 (lognormal sd)
    vmax: float = 10000.0           # ug/ml/h per OD of resistant biomass, at saturation
    km: float = 5.3                 # ug/ml; TEM-1 Km for carbenicillin (14 uM)
    free_enzyme_activity: float = 0.35  # activity of enzyme released by lysed carriers, per OD
                                        # lysed, relative to vmax (Ma et al.'s kappa_b/phi_max)
    enzyme_decay: float = 0.0       # 1/h, inactivation of free enzyme (Ma et al.: none without
                                    # inhibitor)

    # --- protocol --------------------------------------------------------------
    initial_biomass: float = 0.005  # total OD at inoculation, split equally across strains
    pregrowth_h: float = 4.0        # drug-free growth before treatment
    lag_h: float = 0.0              # lag after inoculation and after transfer into fresh
                                    # medium (Baranyi-Roberts); 0 = no lag
    treatment_h: float = 20.0       # with drug
    dilution: float = 100.0         # fold dilution at transfer into fresh medium
    recovery_h: float = 24.0        # drug-free growth after transfer

    # --- measurement -----------------------------------------------------------
    cells_per_OD: float = 1e9       # sets the extinction floor at one cell
    detection_rel: float = 1e-4     # relative abundance below which a barcode is not counted

    # --- numerics --------------------------------------------------------------
    dt_out: float = 0.1             # output spacing, h (does not affect the solution)
    rtol: float = 1e-8
    atol: float = 1e-10

    @property
    def floor(self) -> float:
        """Extinction threshold: one cell."""
        return 1.0 / self.cells_per_OD


BASELINE = Params()


# --------------------------------------------------------------------------- community
@lru_cache(maxsize=256)
def _traits(seed: int, n: int, growth_mean: float, growth_cv: float,
            susceptible_ic50: float, ic50_cv: float):
    """Per-strain maximum growth rate and IC50.

    Strains differ in these two properties and nothing else: they are near-isogenic
    single-gene deletions of one parent, competing for the same nutrient. Both draws
    depend only on the seed, so changing any other parameter compares the same strains.
    """
    rng = np.random.default_rng(seed)
    s = np.sqrt(np.log1p(growth_cv ** 2))
    mu = growth_mean * np.exp(s * rng.standard_normal(n) - s ** 2 / 2)
    ic50 = susceptible_ic50 * np.exp(ic50_cv * rng.standard_normal(n) - ic50_cv ** 2 / 2)
    return mu, ic50


@lru_cache(maxsize=256)
def _carriers(seed: int, n: int, focal_index: int, k: int):
    """Which strains carry the resistance gene.

    k30 is always first, then others in a fixed order from a dedicated random stream,
    so the sets are nested: raising k adds carriers without reshuffling the existing
    ones, and the draw never disturbs the trait sampling.
    """
    mask = np.zeros(n, bool)
    if k <= 0:
        return mask
    others = np.array([i for i in range(n) if i != focal_index])
    order = np.random.default_rng(seed + 991).permutation(others)
    mask[np.r_[focal_index, order[:k - 1]].astype(int)] = True
    return mask


def community(p: Params):
    """Strain traits, carrier mask and starting biomass for one parameter set."""
    mu, ic50 = _traits(p.seed, p.n_strains, p.growth_mean, p.growth_cv,
                       p.susceptible_ic50, p.ic50_cv)
    carrier = _carriers(p.seed, p.n_strains, p.focal_index, p.n_resistant)
    mu = mu * (1 - p.resistance_cost * carrier)   # the enzyme is expressed constantly, so
                                                  # its cost is a constant growth penalty
    n0 = np.full(p.n_strains, p.initial_biomass / p.n_strains)
    n0[p.focal_index] *= p.focal_multiplier
    return mu, ic50, carrier, n0


# --------------------------------------------------------------------------- dynamics
def _lag_factor(t: float, t_fresh: float, p: Params) -> float:
    """Fraction of full growth capacity, t hours after entering fresh medium at t_fresh.

    Baranyi-Roberts adjustment alpha = q/(1+q) with dq/dt = growth_mean * q, written in
    closed form. q(t_fresh) is set so that log biomass lags its unlagged course by `lag_h`
    hours. With lag_h = 0 the factor is exactly 1 and the model is unchanged.
    """
    if p.lag_h <= 0.0:
        return 1.0
    q0_inv = np.expm1(p.growth_mean * p.lag_h)          # 1/q0 = exp(mu*lag) - 1
    return 1.0 / (1.0 + q0_inv * np.exp(-p.growth_mean * (t - t_fresh)))


def _rates(N, R, a, mu, ic50, carrier, p: Params, E: float = 0.0, lag: float = 1.0):
    """Per-strain growth and killing, nutrient draw-down, drug degradation, free enzyme.

    `E` is free enzyme released by lysed carriers, in OD of lysed carrier biomass; `lag` is
    the lag-phase growth factor from `_lag_factor`.
    """
    R = max(R, 0.0)
    a = max(a, 0.0)
    E = max(E, 0.0)
    g = lag * mu * R / (p.half_saturation + R)        # Monod growth, slowed during lag
    if a > 0.0:
        x = (a / ic50) ** p.hill
        kill = p.gamma * (x / (1.0 + x)) * g          # beta-lactam killing needs growth
        kill = np.where(carrier, p.carrier_lysis * kill, kill)   # 0 = immune carriers
    else:
        kill = np.zeros_like(g)
    dR = -float(np.sum(N * g))                        # 1 OD of biomass consumes 1 unit of food
    resistant_biomass = float(np.sum(N[carrier]))
    enzyme = resistant_biomass + p.free_enzyme_activity * E
    da = -p.vmax * enzyme * a / (p.km + a)            # saturating degradation, cells + free
    dE = float(np.sum((N * kill)[carrier])) - p.enzyme_decay * E   # released on lysis
    return g, kill, dR, da, dE


def _segment(N0, R0, a0, E0, t0, t1, t_fresh, mu, ic50, carrier, p: Params, alive):
    """Integrate one phase, stopping whenever a strain crosses the one-cell floor.

    Biomass is integrated in log space so that abundances spanning many orders of
    magnitude stay accurate; a strain that reaches one cell is removed for good.
    `t_fresh` is when the culture last entered fresh medium, which sets the lag phase.
    """
    floor_log = np.log(p.floor)
    grid_t = np.arange(t0, t1 + p.dt_out / 2, p.dt_out)
    ts, Ns, Rs, as_, Es, deaths = [], [], [], [], [], []
    death_total = 0.0
    N, R, a, E = N0.copy(), R0, a0, E0

    while t0 < t1 - 1e-12:
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            te = grid_t[(grid_t >= t0) & (grid_t <= t1)]
            for t in te:
                ts.append(t); Ns.append(N.copy()); Rs.append(R); as_.append(a)
                Es.append(E); deaths.append(death_total)
            break

        def rhs(t, y):
            Nl = np.zeros_like(N)
            Nl[idx] = np.exp(y[:idx.size])
            g, kill, dR, da, dE = _rates(Nl, y[-4], y[-3], mu, ic50, carrier, p,
                                         E=y[-2], lag=_lag_factor(t, t_fresh, p))
            dlog = g[idx] - p.death_rate - kill[idx]
            dloss = float(np.sum(Nl * (p.death_rate + kill)))
            return np.r_[dlog, dR, da, dE, dloss]

        def make_event(j):
            def ev(t, y):
                return y[j] - floor_log
            ev.terminal = True
            ev.direction = -1
            return ev

        y0 = np.r_[np.log(np.maximum(N[idx], p.floor * 1e-3)), R, a, E, death_total]
        te = np.minimum(grid_t[(grid_t >= t0) & (grid_t <= t1)], t1)
        sol = solve_ivp(rhs, (t0, t1), y0, method="LSODA", t_eval=te,
                        events=[make_event(j) for j in range(idx.size)],
                        rtol=p.rtol, atol=p.atol, dense_output=False)
        if not sol.success:
            raise RuntimeError(sol.message)

        for k, t in enumerate(sol.t):
            Nk = np.zeros_like(N)
            Nk[idx] = np.exp(sol.y[:idx.size, k])
            ts.append(float(t)); Ns.append(Nk); Rs.append(float(sol.y[-4, k]))
            as_.append(float(sol.y[-3, k])); Es.append(float(sol.y[-2, k]))
            deaths.append(float(sol.y[-1, k]))

        if sol.status == 1:                                  # a strain hit the floor
            t0 = float(sol.t_events[[len(e) > 0 for e in sol.t_events].index(True)][0])
            yend = np.concatenate([e for e in sol.y_events if len(e)])[0]
            N = np.zeros_like(N); N[idx] = np.exp(yend[:idx.size])
            R, a, E = float(yend[-4]), float(yend[-3]), float(yend[-2])
            death_total = float(yend[-1])
            hit = idx[np.isclose(yend[:idx.size], floor_log, atol=1e-9)]
            if hit.size == 0:
                hit = idx[[int(np.argmin(yend[:idx.size]))]]
            death_total += float(np.sum(N[hit]))             # the last cell is accounted for
            alive[hit] = False
            N[hit] = 0.0
        else:
            t0 = t1
            N = Ns[-1].copy() if Ns else N
            R, a, E, death_total = Rs[-1], as_[-1], Es[-1], deaths[-1]

    return (np.array(ts), np.array(Ns), np.array(Rs), np.array(as_), np.array(Es),
            np.array(deaths), N, R, a, E, death_total)


def simulate(p: Params) -> dict:
    """Pre-growth -> treatment -> 1:dilution transfer -> recovery."""
    mu, ic50, carrier, n0 = community(p)
    alive = n0 > p.floor
    t_treat = p.pregrowth_h
    t_transfer = p.pregrowth_h + p.treatment_h
    t_end = t_transfer + p.recovery_h

    T, NN, RR, AA, EE, DD = [], [], [], [], [], []
    removed = 0.0

    # phase 1: drug-free pre-growth, from inoculation into fresh medium at t = 0
    ts, Ns, Rs, as_, es, dd, N, R, a, E, death = _segment(
        n0, p.food, 0.0, 0.0, 0.0, t_treat, 0.0, mu, ic50, carrier, p, alive)
    T += list(ts); NN += list(Ns); RR += list(Rs); AA += list(as_); EE += list(es)
    DD += list(dd)

    # phase 2: drug added to the same culture (no new lag)
    a = p.dose
    ts, Ns, Rs, as_, es, dd, N, R, a, E, death2 = _segment(
        N, R, a, E, t_treat, t_transfer, 0.0, mu, ic50, carrier, p, alive)
    T += list(ts); NN += list(Ns); RR += list(Rs); AA += list(as_); EE += list(es)
    DD += list(np.asarray(dd) + death)
    death += death2

    # transfer: 1:dilution into fresh drug-free medium
    pre_transfer = N.copy()
    removed += R          # the spent medium, with whatever nutrient is left in it, is discarded
    N = N / p.dilution
    E = E / p.dilution    # free enzyme is carried over in the transferred volume
    removed += float(np.sum(pre_transfer) - np.sum(N))
    gone = (N < p.floor) & alive
    removed += float(np.sum(N[gone]))
    N[gone] = 0.0
    alive = alive & ~gone
    post_transfer = N.copy()

    # phase 3: recovery in fresh medium, with a new lag from the transfer
    ts, Ns, Rs, as_, es, dd, N, R, a, E, death3 = _segment(
        N, p.food, 0.0, E, t_transfer, t_end, t_transfer, mu, ic50, carrier, p, alive)
    T += list(ts); NN += list(Ns); RR += list(Rs); AA += list(as_); EE += list(es)
    DD += list(np.asarray(dd) + death)
    death += death3

    T = np.array(T); NN = np.array(NN)
    # mass balance: food supplied = biomass standing + biomass dead + biomass discarded
    supplied = p.food * 2 + float(np.sum(n0))      # fresh medium twice, plus the inoculum
    standing = float(np.sum(NN[-1])) + float(RR[-1])
    budget_error = standing + death + removed - supplied

    return dict(params=p, t=T, N=NN, R=np.array(RR), a=np.array(AA), E=np.array(EE),
                carrier=carrier, mu=mu, ic50=ic50, t_transfer=t_transfer,
                pre_transfer=pre_transfer, post_transfer=post_transfer,
                budget_error=budget_error)


# --------------------------------------------------------------------------- readouts
def snapshot(sim: dict, hour: float) -> np.ndarray:
    """Absolute abundance of every strain at a sampling time.

    The stored point at the transfer time is the end-of-treatment culture, i.e. the
    24 h sample taken before the community is diluted.
    """
    t = sim["t"]
    hit = np.flatnonzero(np.isclose(t, hour, atol=1e-7))
    return sim["N"][hit[0] if hit.size else int(np.argmin(np.abs(t - hour)))]


def detected(x: np.ndarray, p: Params) -> np.ndarray:
    """Abundances with strains below the sequencing detection limit set to zero."""
    tot = x.sum()
    if tot <= 0:
        return np.zeros_like(x)
    return np.where(x / tot >= p.detection_rel, x, 0.0)


def shannon(x: np.ndarray) -> float:
    """H = -sum p ln p, as reported in the manuscript."""
    s = x.sum()
    if s <= 0:
        return 0.0
    q = x[x > 0] / s
    return float(-(q * np.log(q)).sum())


def richness(x: np.ndarray) -> int:
    return int((x > 0).sum())


def bray_curtis(x: np.ndarray, y: np.ndarray) -> float:
    px = x / x.sum() if x.sum() > 0 else x
    py = y / y.sum() if y.sum() > 0 else y
    return float(np.abs(px - py).sum() / max((px + py).sum(), 1e-300))


@dataclass
class Result:
    """One matched experiment: resistant and sensitive communities, treated and not."""
    params: Params
    runs: dict[tuple[str, bool], dict]      # (community, treated) -> simulation
    metrics: pd.Series                      # the readouts, one row
    timeseries: pd.DataFrame                # per-strain abundance at the sampling hours
    protection: pd.Series                   # per-strain protection score at 24 h

    @property
    def treated(self) -> dict:
        return self.runs[("resistant", True)]

    @property
    def control(self) -> dict:
        return self.runs[("resistant", False)]


def _metrics(runs, p: Params) -> pd.Series:
    carrier = runs[("resistant", True)]["carrier"]
    prot = ~carrier                                  # members without the gene

    def at(key, hour):
        return detected(snapshot(runs[key], hour), p)

    t48 = {k: at(k, 48.0) for k in runs}
    xt, xc = t48[("resistant", True)], t48[("resistant", False)]
    st, sc = t48[("sensitive", True)], t48[("sensitive", False)]

    def ratio(num, den):
        return float(num / den) if den > 0 else np.nan

    # drug clearance: hours from drug addition until the drug first falls below the median
    # susceptible IC50 (where killing of a typical strain is at half its maximum)
    trt = runs[("resistant", True)]
    # the treatment phase starts at the first stored point carrying drug; the point stored
    # at the same time just before it is the end of drug-free pre-growth
    start = int(np.argmax((trt["t"] >= p.pregrowth_h - 1e-9) & (trt["a"] > 0)))
    stop = int(np.argmax(trt["t"] >= trt["t_transfer"] - 1e-9))     # end of treatment
    below = np.flatnonzero(trt["a"][start:stop + 1] < p.susceptible_ic50)
    clearance = (float(trt["t"][start + below[0]] - p.pregrowth_h)
                 if p.dose > p.susceptible_ic50 and below.size else np.nan)
    i24 = int(np.argmin(np.abs(trt["t"] - 24.0)))
    carriers_24 = snapshot(trt, 24.0)[carrier].sum()
    # early response of the non-carrier members, as in the lab's 4 h -> 6 h biomass panel;
    # 4 h is the last drug-free point, so this is the first two hours of treatment
    early = ratio(snapshot(trt, 6.0)[prot].sum(), snapshot(trt, 4.0)[prot].sum())

    return pd.Series(dict(
        # --- manuscript readouts, whole community --------------------------------
        biomass_treated=float(xt.sum()), biomass_control=float(xc.sum()),
        biomass_resilience=ratio(xt.sum(), xc.sum()),
        shannon_treated=shannon(xt), shannon_control=shannon(xc),
        shannon_resilience=ratio(shannon(xt), shannon(xc)),
        richness_treated=richness(xt), richness_control=richness(xc),
        compositional_resilience=1 - bray_curtis(xt, xc),
        # --- the resistant member ------------------------------------------------
        resistant_share=ratio(xt[carrier].sum(), xt.sum()),
        resistant_biomass=float(xt[carrier].sum()),
        n_resistant=int(carrier.sum()),
        resistant_fraction=float(carrier.sum() / p.n_strains),
        # --- the members the resistance is protecting ----------------------------
        protected_shannon=shannon(xt[prot]), protected_shannon_control=shannon(xc[prot]),
        protected_shannon_ratio=ratio(shannon(xt[prot]), shannon(xc[prot])),
        protected_richness=richness(xt[prot]), protected_richness_control=richness(xc[prot]),
        protected_biomass_ratio=ratio(xt[prot].sum(), xc[prot].sum()),
        protected_change_4_6h=early,
        # truly gone (below one cell), as opposed to merely below the detection limit
        protected_extinct=int((snapshot(runs[("resistant", True)], 48.0)[prot] <= 0).sum()),
        # --- the matched all-sensitive community ---------------------------------
        sensitive_biomass_resilience=ratio(st.sum(), sc.sum()),
        sensitive_shannon_treated=shannon(st), sensitive_shannon_control=shannon(sc),
        sensitive_richness_treated=richness(st),
        # --- mechanism -----------------------------------------------------------
        drug_clearance_h=clearance,
        free_enzyme_at_24h=float(trt["E"][i24]),
        carrier_biomass_24h=float(carriers_24),
        carrier_survived=bool((detected(snapshot(trt, 48.0), p)[carrier] > 0).any()),
        # --- diagnostics ---------------------------------------------------------
        drug_at_24h=float(trt["a"][i24]),
        budget_error=max(abs(r["budget_error"]) for r in runs.values()),
    ))


def _timeseries(runs, p: Params) -> pd.DataFrame:
    rows = []
    for (comm, trt), sim in runs.items():
        for h in SAMPLE_HOURS:
            x = snapshot(sim, h)
            tot = x.sum()
            for i in range(p.n_strains):
                rows.append(dict(community=comm, treated=trt, hour=h, strain=f"k{i + 1:02d}",
                                 carrier=bool(sim["carrier"][i]), absolute=float(x[i]),
                                 relative=float(x[i] / tot) if tot > 0 else 0.0))
    return pd.DataFrame(rows)


def _protection(runs, p: Params) -> pd.Series:
    """log2[(treated/untreated with resistance) / (treated/untreated without)] at 24 h.

    A one-cell pseudocount keeps strains that were removed from dividing by zero.
    """
    fl = p.floor

    def x(key):
        return snapshot(runs[key], 24.0) + fl

    score = np.log2((x(("resistant", True)) / x(("resistant", False))) /
                    (x(("sensitive", True)) / x(("sensitive", False))))
    carrier = runs[("resistant", True)]["carrier"]
    names = np.array([f"k{i + 1:02d}" for i in range(p.n_strains)])
    return pd.Series(score[~carrier], index=names[~carrier], name="log2_protection")


@lru_cache(maxsize=512)
def _cached_sim(p: Params) -> dict:
    return simulate(p)


def run(params: Params | None = None, **changes: Any) -> Result:
    """Run one matched experiment and return every readout.

    Anything not named keeps its baseline value::

        run()                       # the baseline experiment
        run(growth_mean=0.2)        # slower community, everything else unchanged
        run(dose=70, n_resistant=4)

    The four runs are the manuscript's design: resistant and all-sensitive communities,
    each treated and untreated. The untreated runs carry no drug, so they are cached and
    reused across a sweep.
    """
    p = replace(params or BASELINE, **changes)
    if p.n_resistant <= 0:
        raise ValueError("n_resistant must be at least 1; the all-sensitive community "
                         "is built automatically as the matched control")
    sensitive = replace(p, n_resistant=0)
    runs = {
        ("resistant", True): _cached_sim(p),
        ("resistant", False): _cached_sim(replace(p, dose=0.0)),
        ("sensitive", True): _cached_sim(sensitive),
        ("sensitive", False): _cached_sim(replace(sensitive, dose=0.0)),
    }
    return Result(params=p, runs=runs, metrics=_metrics(runs, p),
                  timeseries=_timeseries(runs, p), protection=_protection(runs, p))


# --------------------------------------------------------------------------- sweeps
def _row(changes: dict) -> pd.Series:
    """One sweep row: the parameters that were changed, then every readout.

    A changed parameter that is also a readout name (`n_resistant`, say) is kept only
    once, so the row index stays unique.
    """
    r = run(**changes)
    extra = pd.Series({"protection_median": float(r.protection.median()),
                       "protection_min": float(r.protection.min()),
                       "protection_max": float(r.protection.max())})
    head = pd.Series({k: v for k, v in changes.items()
                      if k not in r.metrics.index and k not in extra.index})
    return pd.concat([head, r.metrics, extra])


def sweep(parameter: str, values, cache: str | None = None,
          recompute: bool = False, **fixed: Any) -> pd.DataFrame:
    """Vary one parameter over a list of values; everything else stays at baseline.

        sweep("n_resistant", [1, 2, 4, 8, 16], vmax=500)

    `cache` names a CSV: the sweep is computed once and reloaded afterwards. Delete the
    file, or pass `recompute=True`, to run it again.
    """
    if cache and Path(cache).exists() and not recompute:
        return pd.read_csv(cache)
    df = pd.DataFrame([_row({parameter: v, **fixed}) for v in values])
    if cache:
        df.to_csv(cache, index=False)
    return df


def grid(p1: str, v1, p2: str, v2, cache: str | None = None,
         recompute: bool = False, **fixed: Any) -> pd.DataFrame:
    """Vary two parameters over a rectangle, for a heatmap.

        grid("gamma", [1, 2, 4], "vmax", [100, 1000], cache="panelC.csv")
    """
    if cache and Path(cache).exists() and not recompute:
        return pd.read_csv(cache)
    df = pd.DataFrame([_row({p1: a, p2: b, **fixed}) for a in v1 for b in v2])
    if cache:
        df.to_csv(cache, index=False)
    return df


# --------------------------------------------------------------------------- checks
TESTED_RANGES = {
    # parameter: (low, high) - the span over which the model has been exercised and
    # behaves sensibly. Outside these the solver may struggle or the biology stops
    # making sense; see README.md.
    "growth_mean": (0.02, 1.8),
    "growth_cv": (0.0, 0.6),
    "food": (0.2, 5.0),
    "half_saturation": (0.01, 0.5),
    "death_rate": (0.0005, 0.02),
    "n_resistant": (1, 57),
    "focal_multiplier": (1.0, 30.0),
    "resistance_cost": (0.0, 0.5),
    "dose": (5.0, 500.0),
    "gamma": (0.3, 30.0),
    "hill": (1.0, 4.0),
    "susceptible_ic50": (1.0, 20.0),
    "ic50_cv": (0.1, 1.2),
    "vmax": (10.0, 1e5),
    "km": (0.5, 50.0),
    "initial_biomass": (1e-4, 0.05),
    "dilution": (10.0, 1000.0),
    "carrier_lysis": (0.0, 1.0),
    "free_enzyme_activity": (0.0, 10.0),
    "enzyme_decay": (0.0, 2.0),
    "lag_h": (0.0, 4.0),
}


def check(verbose: bool = True) -> pd.DataFrame:
    """Sanity checks on the model's behaviour. Run once per session, not on import."""
    rows = []

    def add(name, ok, detail):
        rows.append(dict(check=name, passed=bool(ok), detail=detail))

    base = run()

    # 1. mass balance
    worst = base.metrics.budget_error
    add("Mass balance closes: food supplied = biomass + dead + discarded",
        worst < 1e-6, f"max |error| = {worst:.1e}")

    # 2. no negative state
    mn = min(min(r["N"].min(), r["R"].min(), r["a"].min()) for r in base.runs.values())
    add("No negative biomass, nutrient or drug", mn > -1e-8, f"min = {mn:.1e}")

    # 3. zero dose is the untreated control
    z = run(dose=0.0)
    d = float(np.max(np.abs(snapshot(z.treated, 48.0) - snapshot(z.control, 48.0))))
    add("Dose 0 reproduces the untreated control exactly", d < 1e-12, f"max difference = {d:.1e}")

    # 4. no resistance -> resistant and sensitive communities are identical
    one = run(n_resistant=1, resistance_cost=0.0, vmax=0.0)
    a1 = snapshot(one.runs[("resistant", True)], 48.0)
    a2 = snapshot(one.runs[("sensitive", True)], 48.0)
    d = float(np.max(np.abs(a1 - a2)))
    add("A carrier that neither degrades nor pays a cost still differs only by its immunity",
        a1[one.params.focal_index] >= a2[one.params.focal_index] - 1e-12,
        f"max difference across strains = {d:.1e} (k30 immune, so it can only do better)")

    # 5. killing needs growth
    stat = run(pregrowth_h=40.0)
    h = stat.params.pregrowth_h + stat.params.treatment_h
    ratio = snapshot(stat.runs[("sensitive", True)], h).sum() / \
        snapshot(stat.runs[("sensitive", False)], h).sum()
    add("A drug added in stationary phase kills almost nobody",
        ratio > 0.95, f"treated/untreated biomass = {ratio:.4f}")

    # 6. the transfer is an exact split
    c = base.control
    pre, post = c["pre_transfer"], c["post_transfer"]
    keep = pre / base.params.dilution >= base.params.floor
    add("The transfer is an exact 1:dilution split of the survivors",
        np.allclose(post[keep], pre[keep] / base.params.dilution, rtol=1e-10) and
        np.all(post[~keep] == 0.0),
        f"1:{base.params.dilution:g}; {int((~keep).sum())} strains fell below one cell "
        f"and were removed")

    # 7. extinction is irreversible
    bad = 0
    for r in base.runs.values():
        ever_zero = (r["N"] <= 0).any(axis=0)
        for j in np.flatnonzero(ever_zero):
            first = int(np.argmax(r["N"][:, j] <= 0))
            if np.any(r["N"][first:, j] > 0):
                bad += 1
    add("A strain removed at the one-cell floor never returns", bad == 0,
        f"{bad} strains reappeared")

    # 8. carriers are immune
    hard = run(dose=500.0, gamma=20.0, vmax=0.0)
    x = snapshot(hard.treated, 48.0)
    add("Resistant strains survive a dose that removes the susceptible ones",
        x[hard.runs[("resistant", True)]["carrier"]].sum() > 0,
        f"k30 biomass {x[hard.params.focal_index]:.3g} OD, "
        f"susceptible richness {richness(detected(x, hard.params)) - 1}")

    # 9. more resistant biomass clears the drug faster
    slow = [run(n_resistant=k, vmax=20.0).metrics.drug_at_24h for k in (1, 8)]
    add("More resistant biomass leaves less drug behind", slow[1] <= slow[0] + 1e-9,
        f"drug at 24 h: 1 carrier {slow[0]:.3g}, 8 carriers {slow[1]:.3g} ug/ml")

    # 10. saturating degradation makes clearance time scale with dose
    # 10. degradation saturates: well above Km the removal rate stops tracking the dose.
    # This is what keeps the dose dependence in the model; with first-order degradation
    # the rate would be proportional to the dose and the two ratios below would match.
    p0 = BASELINE
    mu_, ic50_, carrier_, n0_ = community(p0)
    N_ = np.where(carrier_, 0.01, 0.0)                 # a fixed amount of resistant biomass

    def removal(a):
        return -_rates(N_, p0.food, a, mu_, ic50_, carrier_, p0)[3]
    hi = removal(10 * p0.km) / removal(5 * p0.km)      # both well above Km -> nearly flat
    lo = removal(0.1 * p0.km) / removal(0.05 * p0.km)  # both well below Km -> nearly linear
    add("Degradation saturates above Km and is first-order below it",
        hi < 1.2 and lo > 1.8,
        f"doubling the drug multiplies the removal rate by {hi:.2f} above Km "
        f"and {lo:.2f} below it (first-order everywhere would give 2.00)")

    # 11. resilience of an untreated run against itself is 1
    add("Untreated against untreated gives resilience 1",
        abs(1 - bray_curtis(snapshot(c, 48.0), snapshot(c, 48.0)) - 1.0) < 1e-12,
        "compositional resilience = 1.000")

    # 12. It is degradation that rescues the community, not the mere presence of a
    # surviving resistant strain. A carrier that degrades nothing still competes for
    # nutrient, which slows everyone's growth and therefore slows killing, so a small
    # positive protection score remains. That residual must be far weaker than enzymatic
    # degradation, and it must not rescue the susceptible members.
    no_enz = run(vmax=0.0)
    with_enz = base
    add("Degradation, not the survival of a resistant strain, rescues the community",
        float(with_enz.protection.median()) > 2 * float(no_enz.protection.median()) and
        with_enz.metrics.protected_richness > no_enz.metrics.protected_richness + 10,
        f"median protection {with_enz.protection.median():.1f} with the enzyme vs "
        f"{no_enz.protection.median():.1f} without; susceptible survivors "
        f"{with_enz.metrics.protected_richness:.0f} vs {no_enz.metrics.protected_richness:.0f} of 56")

    # 13. Shannon cannot exceed ln(number of strains)
    worstH = max(shannon(snapshot(r, 4.0)) for r in base.runs.values())
    cap = np.log(base.params.n_strains)
    add("Shannon diversity stays below ln(n_strains)", worstH <= cap + 1e-9,
        f"max H = {worstH:.3f}, ln({base.params.n_strains}) = {cap:.3f}")

    # 14. solver convergence
    tight = run(rtol=1e-10, atol=1e-12)
    rel = float(np.max(np.abs(snapshot(base.treated, 48.0) - snapshot(tight.treated, 48.0))) /
                max(snapshot(tight.treated, 48.0).max(), 1e-300))
    add("Tighter solver tolerances give the same answer", rel < 1e-4,
        f"max relative difference = {rel:.1e}")

    # 15. results are not an artefact of one seed
    direction = [run(seed=s).metrics.protected_shannon_ratio >=
                 run(seed=s, vmax=0.0).metrics.protected_shannon_ratio for s in (1, 7, 42, 99)]
    add("Degradation helps the susceptible members for every seed tested",
        all(direction), f"{sum(direction)}/4 seeds")

    # 16. a carrier stripped of immunity, degradation, release and cost is just another
    # susceptible strain, so the resistant community must match the all-sensitive one
    plain = run(carrier_lysis=1.0, vmax=0.0, free_enzyme_activity=0.0, resistance_cost=0.0)
    a1 = snapshot(plain.runs[("resistant", True)], 48.0)
    a2 = snapshot(plain.runs[("sensitive", True)], 48.0)
    a1_24 = snapshot(plain.runs[("resistant", True)], 24.0)
    a2_24 = snapshot(plain.runs[("sensitive", True)], 24.0)
    d = float(np.max(np.abs(a1_24 - a2_24)) / max(a2_24.max(), 1e-300))
    add("A carrier with no immunity, degradation, release or cost is an ordinary strain",
        d < 1e-6 and np.array_equal(a1 > 0, a2 > 0),
        f"max relative difference at 24 h = {d:.1e}")

    # 17. free enzyme exists only when carriers can be lysed
    e_base = max(float(r["E"].max()) for r in base.runs.values())
    lysed = run(carrier_lysis=1.0)
    e_lys = float(lysed.treated["E"].max())
    no_rel = run(carrier_lysis=1.0, free_enzyme_activity=0.0)
    add("Free enzyme appears only when carriers are lysed, and only it can help then",
        e_base == 0.0 and e_lys > 0.0 and
        lysed.metrics.protected_richness >= no_rel.metrics.protected_richness,
        f"max free enzyme: immune carriers {e_base:.1e}, lysed carriers {e_lys:.1e} OD; "
        f"susceptible survivors with/without release {lysed.metrics.protected_richness:.0f}/"
        f"{no_rel.metrics.protected_richness:.0f}")

    # 18. losing private protection costs the carrier: monotone in carrier_lysis
    cb = [run(carrier_lysis=b).metrics.carrier_biomass_24h for b in (0.0, 0.3, 0.6, 0.9)]
    add("The less private protection, the less carrier biomass at 24 h",
        all(x >= y for x, y in zip(cb, cb[1:])),
        "carrier OD at 24 h for carrier_lysis 0, 0.3, 0.6, 0.9: " +
        ", ".join(f"{x:.2g}" for x in cb))

    # 19. lag follows the Baranyi-Roberts curve: while nutrient is plentiful, the log-biomass
    # deficit at time t is -ln[exp(-mu t)(1 - exp(-mu L)) + exp(-mu L)], which tends to mu L
    lag, t4, mu_ = 1.0, 4.0, BASELINE.growth_mean
    n_nolag = snapshot(run(dose=0.0).control, t4).sum()
    n_lag = snapshot(run(dose=0.0, lag_h=lag).control, t4).sum()
    expected = -np.log(np.exp(-mu_ * t4) * -np.expm1(-mu_ * lag) + np.exp(-mu_ * lag))
    shift = np.log(n_nolag / n_lag) / expected
    add("Lag follows the Baranyi-Roberts curve",
        0.9 < shift < 1.1,
        f"log-biomass deficit at 4 h for a 1 h lag = {shift:.3f} x the Baranyi prediction "
        f"(strains' growth rates differ from the mean, so not exactly 1)")

    # 20. cells in lag are drug-tolerant, because killing needs growth
    def early_survival(lag_h):
        s = run(pregrowth_h=0.5, lag_h=lag_h, vmax=0.0).runs
        return snapshot(s[("sensitive", True)], 2.5).sum() / \
            snapshot(s[("sensitive", False)], 2.5).sum()
    s0, s3 = early_survival(0.0), early_survival(3.0)
    add("Cells still in lag when the drug arrives are killed less",
        s3 > s0, f"treated/untreated biomass 2 h into treatment: no lag {s0:.3f}, "
        f"3 h lag {s3:.3f}")

    df = pd.DataFrame(rows)
    if verbose:
        with pd.option_context("display.width", 200, "display.max_colwidth", 95):
            print(df.to_string(index=False))
        print(f"\n{int(df.passed.sum())}/{len(df)} checks passed")
    return df


if __name__ == "__main__":
    check()
