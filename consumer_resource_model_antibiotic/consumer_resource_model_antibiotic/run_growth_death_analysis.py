"""Mechanistic analysis of growth, antibiotic death, and net growth for representative bacteria."""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model import Config
from antibiotic_analysis_common import (
    add_common_cli, ensure_dir, representative_species, run_condition, save_figure,
)


def main():
    parser = add_common_cli(argparse.ArgumentParser(description=__doc__), 'analysis_growth_death')
    parser.add_argument('--dose', type=float, default=1.0)
    parser.add_argument('--decay', type=float, default=0.02)
    args = parser.parse_args()
    out = ensure_dir(args.output)

    base = Config(seed=args.seed, initial_antibiotic=args.dose, antibiotic_decay=args.decay)
    p, r = run_condition(base)
    reps = representative_species(p)

    rows = []
    fig, axs = plt.subplots(len(reps), 2, figsize=(12, 3.5 * len(reps)), squeeze=False)
    for row_idx, (label, i) in enumerate(reps.items()):
        species = p['table']['species'].iloc[i]
        ic50 = p['IC50'][i]
        ax = axs[row_idx, 0]
        ax.plot(r['t'], r['gross'][:, i], label='Resource-driven growth g_i')
        ax.plot(r['t'], r['death_rates'][:, i], label='Antibiotic death d_i')
        ax.plot(r['t'], r['net'][:, i], label='Net g_i - d_i')
        ax.axhline(0, ls='--', lw=1)
        ax.set(xlabel='Time', ylabel='Per-capita rate', title=f'{label}: {species}, IC50={ic50:.3g}')
        ax.legend(fontsize=8)

        ax2 = axs[row_idx, 1]
        ax2.plot(r['t'], np.maximum(r['N'][:, i], r['cutoff'] / 10), label='Species biomass')
        ax2.set_yscale('log')
        ax2.set(xlabel='Time', ylabel='Biomass', title=f'{species} biomass under A0={args.dose:g}')
        axA = ax2.twinx()
        axA.plot(r['t'], r['A'], ls='--', label='Antibiotic A(t)')
        axA.set_ylabel('Antibiotic concentration')

        for k, time in enumerate(r['t']):
            rows.append({
                'resistance_group': label, 'species': species, 'IC50': ic50, 'time': time,
                'biomass': r['N'][k, i], 'gross_growth_rate': r['gross'][k, i],
                'antibiotic_death_rate': r['death_rates'][k, i], 'net_growth_rate': r['net'][k, i],
                'antibiotic': r['A'][k],
            })

    save_figure(fig, out / '01_representative_growth_death_balance.png')
    pd.DataFrame(rows).to_csv(out / 'representative_species_timeseries.csv', index=False)

    rep_table = p['table'].iloc[list(reps.values())].copy()
    rep_table.insert(0, 'resistance_group', list(reps.keys()))
    rep_table.to_csv(out / 'representative_species_traits.csv', index=False)
    print(rep_table[['resistance_group', 'species', 'IC50', 'mu_max', 'can_consume_R0', 'can_consume_M1', 'can_consume_M2']].to_string(index=False))
    print(f'Wrote mechanistic growth/death analysis to: {out.resolve()}')


if __name__ == '__main__':
    main()
