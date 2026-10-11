"""Reproducible global sensitivity screen for antibiotic/ecological parameters.

This is an exploratory parameter-importance screen, not a formal Sobol analysis.
It samples parameter sets, runs the model, then reports Spearman rank association
between each input and each model outcome.
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from model import Config
from antibiotic_analysis_common import add_common_cli, community_metrics, ensure_dir, run_condition, save_figure

INPUTS = ['initial_antibiotic', 'antibiotic_decay', 'antibiotic_kmax', 'antibiotic_hill', 'initial_resource', 'leakage', 'resistance_cost_max']
OUTPUTS = ['living_energy_fraction', 'survivors', 'final_shannon', 'final_weighted_IC50', 'cumulative_antibiotic_death_loss']


def log_uniform(rng, low, high, n):
    return np.exp(rng.uniform(np.log(low), np.log(high), n))


def main():
    parser = add_common_cli(argparse.ArgumentParser(description=__doc__), 'analysis_sensitivity')
    parser.add_argument('--samples', type=int, default=40, help='Number of random parameter sets (default: 40)')
    args = parser.parse_args()
    out = ensure_dir(args.output)

    rng = np.random.default_rng(args.seed + 900000)
    n = args.samples
    design = pd.DataFrame({
        'initial_antibiotic': rng.uniform(0.0, 2.0, n),
        'antibiotic_decay': log_uniform(rng, 0.005, 0.10, n),
        'antibiotic_kmax': rng.uniform(0.05, 0.40, n),
        'antibiotic_hill': rng.uniform(1.0, 4.0, n),
        'initial_resource': rng.uniform(0.5, 2.0, n),
        'leakage': rng.uniform(0.05, 0.40, n),
        'resistance_cost_max': rng.uniform(0.0, 0.40, n),
    })

    base = Config(seed=args.seed)
    rows = []
    for run_id, row in design.iterrows():
        changes = {k: float(row[k]) for k in INPUTS}
        p, r = run_condition(base, **changes)
        metrics = community_metrics(p, r)
        rows.append({'run_id': run_id, **changes, **metrics})
        print(f'Completed sensitivity run {run_id + 1}/{n}', end='\r', flush=True)
    print()
    results = pd.DataFrame(rows)
    results.to_csv(out / 'sensitivity_runs.csv', index=False)

    corr = pd.DataFrame(index=INPUTS, columns=OUTPUTS, dtype=float)
    pvals = pd.DataFrame(index=INPUTS, columns=OUTPUTS, dtype=float)
    for x in INPUTS:
        for y in OUTPUTS:
            rho, pval = spearmanr(results[x], results[y])
            corr.loc[x, y] = rho
            pvals.loc[x, y] = pval
    corr.to_csv(out / 'sensitivity_spearman_rho.csv')
    pvals.to_csv(out / 'sensitivity_spearman_pvalues.csv')

    fig, ax = plt.subplots(figsize=(10, 5.5))
    im = ax.imshow(corr.to_numpy(dtype=float), vmin=-1, vmax=1, aspect='auto', interpolation='nearest')
    ax.set(xticks=np.arange(len(OUTPUTS)), xticklabels=[x.replace('_', ' ') for x in OUTPUTS],
           yticks=np.arange(len(INPUTS)), yticklabels=[x.replace('_', ' ') for x in INPUTS],
           title='Global sensitivity screen: Spearman rank correlations')
    ax.tick_params(axis='x', rotation=30)
    for i in range(len(INPUTS)):
        for j in range(len(OUTPUTS)):
            ax.text(j, i, f'{corr.iloc[i, j]:.2f}', ha='center', va='center', fontsize=8)
    fig.colorbar(im, ax=ax, label='Spearman rho')
    save_figure(fig, out / '01_sensitivity_correlation_heatmap.png')

    # Show strongest drivers of two key outcomes.
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, outcome in zip(axs, ['living_energy_fraction', 'final_weighted_IC50']):
        s = corr[outcome].sort_values()
        ax.barh(s.index, s.values)
        ax.axvline(0, lw=1)
        ax.set(xlabel='Spearman rho', title=outcome.replace('_', ' '), xlim=(-1, 1))
    save_figure(fig, out / '02_sensitivity_key_outcomes.png')

    print('\nSpearman rank correlations:')
    print(corr.round(3).to_string())
    print(f'Wrote sensitivity analysis to: {out.resolve()}')


if __name__ == '__main__':
    main()
