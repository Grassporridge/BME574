"""Runs every scan for test_extinction_conditions/ and writes the CSVs. Safe to re-run (cached)."""
import itertools
import sys
from pathlib import Path

from scan import run_jobs, self_test

OUT = Path(__file__).resolve().parent / "extinction_results"

# ---------------------------------------------------------------- the three factors
DOSES = [5, 7, 10, 14, 20, 28, 35, 50, 70, 100, 140, 200, 280, 400, 500]       # ug/ml
EXPOSURES = [0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20]      # h (20 = protocol)
VMAXES = [0, 3, 10, 30, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000]
CARRIERS = [1, 2, 5, 10, 20]
FLOORS = {"one_cell": 1e-9, "1e-12": 1e-12}                                     # OD

# the fixed values used when a factor is NOT being varied (the protocol / model defaults)
FIXED = dict(dose=50.0, exposure_h=20.0, vmax=10000.0)


def grid_jobs(floor):
    return [dict(dose=float(d), exposure_h=float(t), vmax=float(v), n_resistant=k, floor=floor)
            for k, d, t, v in itertools.product(CARRIERS, DOSES, EXPOSURES, VMAXES)]


def sensitive_jobs(floor):
    # the all-sensitive community (no carrier at all) for reference: n_resistant = 0
    return [dict(dose=float(d), exposure_h=float(t), vmax=0.0, n_resistant=0, floor=floor)
            for d, t in itertools.product(DOSES, EXPOSURES)]


if __name__ == "__main__":
    self_test()
    recompute = "--recompute" in sys.argv
    for name, fl in FLOORS.items():
        g = run_jobs(grid_jobs(fl), OUT / f"grid_{name}.csv", recompute=recompute)
        s = run_jobs(sensitive_jobs(fl), OUT / f"all_sensitive_{name}.csv", recompute=recompute)
        print(name, len(g), "grid runs,", len(s), "all-sensitive runs")
