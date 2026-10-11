"""Antibiotic persistence/decay analysis at fixed initial dose."""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model import Config
from antibiotic_analysis_common import (
    add_common_cli, assert_fixed_community, community_metrics, community_timeseries,
    ensure_dir, run_condition, save_figure,
)

DEFAULT_DECAYS = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10]


def main():
    parser = add_common_cli(argparse.ArgumentParser(description=__doc__), 'analysis_decay')
    parser.add_argument('--decays', type=float, nargs='+', default=DEFAULT_DECAYS)
    parser.add_argument('--dose', type=float, default=1.0)
    args = parser.parse_args()
    out = ensure_dir(args.output)

    base = Config(seed=args.seed, initial_antibiotic=args.dose)
    rows = []; community_rows = []; runs = []; ref = None
    for kdeg in args.decays:
        p, r = run_condition(base, antibiotic_decay=kdeg)
        if ref is None: ref = p
        else: assert_fixed_community(ref, p)
        row = community_metrics(p, r)
        row['decay_rate'] = kdeg
        row['half_life'] = np.inf if kdeg == 0 else np.log(2) / kdeg
        rows.append(row)
        community_rows.append(community_timeseries(f'kdeg={kdeg:g}', p, r).assign(decay_rate=kdeg))
        runs.append((kdeg, p, r))

    df = pd.DataFrame(rows).sort_values('decay_rate')
    df.to_csv(out / 'decay_summary.csv', index=False)
    pd.concat(community_rows, ignore_index=True).to_csv(out / 'decay_community_timeseries.csv', index=False)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for kdeg, p, r in runs:
        axs[0].plot(r['t'], r['A'], label=f'kdeg={kdeg:g}')
        axs[1].plot(r['t'], r['N'].sum(axis=1), label=f'kdeg={kdeg:g}')
    axs[0].set(xlabel='Time', ylabel='Antibiotic concentration', title=f'Antibiotic decay (A0={args.dose:g})')
    axs[1].set(xlabel='Time', ylabel='Total living biomass', title='Population response to exposure duration')
    axs[0].legend(fontsize=8); axs[1].legend(fontsize=8)
    save_figure(fig, out / '01_decay_timecourses.png')

    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    x = df['decay_rate']
    axs[0].plot(x, df['final_biomass'], marker='o'); axs[0].set(ylabel='Final biomass')
    axs[1].plot(x, df['survivors'], marker='o'); axs[1].set(ylabel='Surviving species')
    axs[2].plot(x, df['final_weighted_IC50'], marker='o'); axs[2].set(ylabel='Final weighted IC50')
    for ax in axs: ax.set(xlabel='Antibiotic decay rate kdeg')
    save_figure(fig, out / '02_decay_endpoints.png')

    print(df[['decay_rate', 'half_life', 'final_biomass', 'survivors', 'final_weighted_IC50']].to_string(index=False))
    print(f'Wrote decay analysis to: {out.resolve()}')


if __name__ == '__main__':
    main()
