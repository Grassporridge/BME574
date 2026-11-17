"""
When does antibiotic treatment permanently remove susceptible strains?
======================================================================

Simulation code for test_extinction_conditions/. Three factors are varied:

    dose        carbenicillin added at 4 h (ug/ml)
    exposure_h  hours the drug is present before it is removed (see below)
    vmax        per-carrier drug degradation rate (ug/ml/h per OD of carrier biomass)

and the readout is the number of susceptible (non-carrier) strains that are EXTINCT at the end of
the experiment: driven below the extinction floor, after which the model removes them for good.

Exposure time
-------------
In the protocol the drug is removed by the 1:100 transfer at 24 h. Changing `treatment_h` in the
model moves that transfer (and the whole schedule) with it, so it changes two things at once:
how long the drug acts, and when the transfer bottleneck happens. To vary exposure alone, this
module removes the drug at 4 h + exposure_h (all remaining drug set to zero, cells and medium
untouched, e.g. by adding excess beta-lactamase or a wash) and keeps the schedule fixed:
transfer at 24 h, readout at 48 h. exposure_h >= 20 is the unmodified protocol.

Nothing in resilience_model.py is changed; this reuses its integrator.
"""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from resilience_model import BASELINE, Params, community, _segment, simulate  # noqa: E402

PROTOCOL_EXPOSURE = BASELINE.treatment_h            # 20 h: drug present until the transfer


def simulate_exposure(p: Params, exposure_h: float) -> dict:
    """Pre-growth -> drug for `exposure_h` -> drug removed -> transfer at the usual time -> recovery.

    With exposure_h >= p.treatment_h this is exactly resilience_model.simulate(p).
    """
    mu, ic50, carrier, n0 = community(p)
    alive = n0 > p.floor
    t_treat = p.pregrowth_h
    t_transfer = p.pregrowth_h + p.treatment_h
    t_end = t_transfer + p.recovery_h
    t_off = min(t_treat + max(exposure_h, 0.0), t_transfer)

    # phase 1: drug-free pre-growth
    *_, N, R, a, E, _ = _segment(n0, p.food, 0.0, 0.0, 0.0, t_treat, 0.0, mu, ic50, carrier, p, alive)
    sus_min = N[~carrier].copy()

    # phase 2a: drug present
    ts, Ns, _, As, *_, N, R, a, E, _ = _segment(N, R, p.dose, E, t_treat, t_off, 0.0,
                                                 mu, ic50, carrier, p, alive)
    if len(Ns):
        sus_min = np.minimum(sus_min, np.min(np.asarray(Ns)[:, ~carrier], axis=0))
    below = np.flatnonzero(np.asarray(As) < np.median(ic50)) if len(As) else np.array([])
    clearance = float(ts[below[0]] - t_treat) if below.size and p.dose > np.median(ic50) else np.nan
    drug_at_removal = float(a)

    # phase 2b: drug removed, same culture, until the transfer
    if t_off < t_transfer - 1e-12:
        _, Ns, *_, N, R, a, E, _ = _segment(N, R, 0.0, E, t_off, t_transfer, 0.0,
                                            mu, ic50, carrier, p, alive)
        sus_min = np.minimum(sus_min, np.min(np.asarray(Ns)[:, ~carrier], axis=0))

    # transfer: 1:dilution into fresh drug-free medium; strains below one floor unit are lost
    N = N / p.dilution
    E = E / p.dilution
    gone = (N < p.floor) & alive
    N[gone] = 0.0
    alive = alive & ~gone
    sus_min = np.minimum(sus_min, N[~carrier])

    # phase 3: recovery
    *_, N, R, a, E, _ = _segment(N, p.food, 0.0, E, t_transfer, t_end, t_transfer,
                                 mu, ic50, carrier, p, alive)
    return dict(N_end=N, carrier=carrier, mu=mu, ic50=ic50, sus_min=sus_min,
                clearance_h=clearance, drug_at_removal=drug_at_removal)


def one(job: dict) -> dict:
    """One treated run -> one result row. `job` holds dose, exposure_h, vmax, n_resistant, floor."""
    p = replace(BASELINE, dose=float(job["dose"]), vmax=float(job["vmax"]),
                n_resistant=int(job["n_resistant"]), cells_per_OD=1.0 / job["floor"])
    s = simulate_exposure(p, job["exposure_h"])
    car = s["carrier"]
    n_sus = int((~car).sum())
    extinct = int((s["N_end"][~car] <= 0).sum())
    tot = s["N_end"].sum()
    return {**job, "n_susceptible": n_sus, "extinct": extinct,
            "extinct_frac": extinct / n_sus if n_sus else np.nan,
            "surviving": n_sus - extinct,
            # alive but below the 1e-4 sequencing detection limit (swamped, not extinct)
            "undetected_alive": int(((s["N_end"] > 0) & (s["N_end"] < p.detection_rel * tot))[~car].sum()),
            "carrier_share": float(s["N_end"][car].sum() / tot) if tot > 0 else np.nan,
            "min_susceptible_od": float(np.min(s["sus_min"])),
            "clearance_h": s["clearance_h"], "drug_at_removal": s["drug_at_removal"]}


def run_jobs(jobs: list[dict], cache: Path | None = None, workers: int | None = None,
             recompute: bool = False) -> pd.DataFrame:
    """Run many jobs in parallel; cached to CSV."""
    if cache is not None and cache.exists() and not recompute:
        return pd.read_csv(cache)
    workers = workers or max(1, (os.cpu_count() or 2) - 1)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        rows = list(ex.map(one, jobs, chunksize=4))
    df = pd.DataFrame(rows)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(cache, index=False)
    return df


def self_test() -> None:
    """exposure >= 20 h must reproduce the model's own simulate() exactly; no drug -> no extinction."""
    for ch in [{}, dict(dose=70.0), dict(vmax=0.0), dict(n_resistant=5, vmax=300.0)]:
        p = replace(BASELINE, **ch)
        ref = simulate(p)["N"][-1]
        mine = simulate_exposure(p, 20.0)["N_end"]
        assert np.allclose(mine, ref, rtol=1e-9, atol=1e-15), ch
        assert np.array_equal(mine > 0, ref > 0), ch
    # untreated: nobody goes extinct, so any extinction below is caused by the drug
    for k in (1, 2, 5, 10, 20):
        z = simulate_exposure(replace(BASELINE, dose=0.0, n_resistant=k), 20.0)
        assert (z["N_end"] > 0).all(), k
    # removing the drug at 0 h is the same as never adding it
    z0 = simulate_exposure(replace(BASELINE, dose=0.0), 20.0)["N_end"]
    e0 = simulate_exposure(BASELINE, 0.0)["N_end"]
    assert np.allclose(z0, e0, rtol=1e-9), "exposure 0 should equal no drug"
    print("scan.self_test: exposure >= 20 h reproduces simulate(); untreated runs lose no strain; "
          "exposure 0 = no drug")


if __name__ == "__main__":
    self_test()
