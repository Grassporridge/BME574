"""No-antibiotic control: resource competition and cross-feeding only."""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model import Config, export
from antibiotic_analysis_common import (
    add_common_cli, community_metrics, community_timeseries, ensure_dir,
    run_condition, save_figure, species_outcomes, weighted_mean_trait,
)


def main():
    parser = add_common_cli(argparse.ArgumentParser(description=__doc__), 'analysis_control')
    args = parser.parse_args()
    out = ensure_dir(args.output)

    base = Config(seed=args.seed)
    p, r = run_condition(base, initial_antibiotic=0.0)

    # Full export is useful for the reference control.
    export(p, r, out / 'control_run')
    pd.DataFrame([community_metrics(p, r)]).to_csv(out / 'control_summary.csv', index=False)
    species_outcomes('control', p, r).to_csv(out / 'control_species_outcomes.csv', index=False)
    community_timeseries('control', p, r).to_csv(out / 'control_community_timeseries.csv', index=False)

    design_cols = [
        'species', 'food_category', 'food_category_code', 'resistance', 'IC50',
        'growth_score', 'intrinsic_mu_max', 'resistance_fitness_cost',
        'resistance_growth_multiplier', 'mu_max', 'initial_biomass',
    ]
    p['table'][design_cols].to_csv(out / 'control_species_design.csv', index=False)

    # Design check: balanced feeding categories and explicit resistance-growth cost.
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    category_order = ['R0', 'M1', 'M2', 'R0+M1', 'R0+M2', 'M1+M2', 'R0+M1+M2']
    counts = p['table']['food_category'].value_counts().reindex(category_order, fill_value=0)
    axs[0].bar(np.arange(len(category_order)), counts.values)
    axs[0].set(xticks=np.arange(len(category_order)), xticklabels=category_order,
               ylabel='Number of species', title='Balanced resource-use categories')
    axs[0].tick_params(axis='x', rotation=35)

    axs[1].scatter(p['table']['resistance'], p['table']['resistance_growth_multiplier'])
    axs[1].set(xlabel='Resistance score (rho)', ylabel='Growth multiplier',
               title=f'Resistance fitness cost (max={p["config"].resistance_cost_max:.2f})')
    save_figure(fig, out / '00_species_design.png')

    t = r['t']; N = r['N']; total = N.sum(axis=1)
    mean_ic50 = weighted_mean_trait(N, p['IC50'])
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    for i in range(N.shape[1]):
        axs[0, 0].plot(t, np.maximum(N[:, i], r['cutoff'] / 10), lw=0.8)
    axs[0, 0].axhline(r['cutoff'], ls='--', lw=1)
    axs[0, 0].set(yscale='log', xlabel='Time', ylabel='Species biomass', title='No-antibiotic species dynamics')

    axs[0, 1].plot(t, total, label='Total biomass')
    for j, name in enumerate(['R0', 'M1', 'M2']):
        axs[0, 1].plot(t, r['R'][:, j], label=name)
    axs[0, 1].set(xlabel='Time', ylabel='Energy-equivalent concentration', title='Biomass and resources')
    axs[0, 1].legend()

    shares = np.divide(100 * N, total[:, None], out=np.zeros_like(N), where=total[:, None] > 0)
    axs[1, 0].stackplot(t, shares.T)
    axs[1, 0].set(xlabel='Time', ylabel='Community composition (%)', ylim=(0, 100), title='Control composition')

    axs[1, 1].plot(t, mean_ic50)
    axs[1, 1].set(xlabel='Time', ylabel='Biomass-weighted mean IC50', title='Resistance trait drift without antibiotic')
    save_figure(fig, out / '01_control_dynamics.png')

    print(pd.DataFrame([community_metrics(p, r)]).to_string(index=False))
    print(f'Wrote control analysis to: {out.resolve()}')


if __name__ == '__main__':
    main()
