"""Does extinction move in the intuitive direction along every one-factor slice of the grid?

Expected:
    more degradation (vmax)        -> fewer extinctions   (non-increasing)
    higher dose                    -> more extinctions    (non-decreasing)
    longer exposure                -> more extinctions    (non-decreasing)
    more carriers (same vmax each) -> smaller fraction of susceptibles extinct (non-increasing)
Each violation is listed so it can be investigated.
"""
import sys
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent / "extinction_results"
RULES = [("vmax", "extinct", -1), ("dose", "extinct", +1), ("exposure_h", "extinct", +1),
         ("n_resistant", "extinct_frac", -1)]


def violations(df):
    rows = []
    keys = ["n_resistant", "dose", "exposure_h", "vmax"]
    for factor, col, sign in RULES:
        others = [k for k in keys if k != factor]
        for vals, g in df.groupby(others):
            g = g.sort_values(factor)
            y = g[col].to_numpy() * sign
            x = g[factor].to_numpy()
            for i in range(1, len(y)):
                if y[i] < y[i - 1] - 1e-12:
                    rows.append({"factor": factor, **dict(zip(others, vals)),
                                 f"{factor}_from": x[i - 1], f"{factor}_to": x[i],
                                 f"{col}_from": g[col].iloc[i - 1], f"{col}_to": g[col].iloc[i]})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "one_cell"
    df = pd.read_csv(OUT / f"grid_{name}.csv")
    v = violations(df)
    print(f"{name}: {len(df)} runs; {len(v)} violations of the intuitive direction")
    if len(v):
        with pd.option_context("display.width", 250, "display.max_rows", 200):
            print(v.groupby("factor").size())
            print(v.head(60).to_string(index=False))
    v.to_csv(OUT / f"intuition_violations_{name}.csv", index=False)
