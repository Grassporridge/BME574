"""Test whether antibiotic exposure enriches high-IC50 bacteria in a fixed community."""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from model import Config
from antibiotic_analysis_common import (
    add_common_cli, assert_fixed_community, ensure_dir, run_condition,
    save_figure, weighted_mean_trait,
)

DEFAULT_DOSES = [0.0, 0.25, 0.5, 1.0, 2.0]


def main():
    parser = add_common_cli(argparse.ArgumentParser(description=__doc__), 'analysis_resistance_selection')
    parser.add_argument('--doses', type=float, nargs='+', default=DEFAULT_DOSES)
    args = parser.parse_args()
    out = ensure_dir(args.output)

    base = Config(seed=args.seed)
    runs = {}
    ref = None
    for dose in args.doses:
        p, r = run_condition(base, initial_antibiotic=dose)
        if ref is None: ref = p
        else: assert_fixed_community(ref, p)
        runs[dose] = (p, r)

    if 0.0 not in runs:
        p0, r0 = run_condition(base, initial_antibiotic=0.0)
        runs[0.0] = (p0, r0)
        if ref is not None: assert_fixed_community(ref, p0)

    p0, r0 = runs[0.0]
    control_final = np.maximum(r0['N'][-1], 1e-300)
    rows = []
    corr_rows = []
    for dose in sorted(runs):
        p, r = runs[dose]
        ratio = r['N'][-1] / control_final
        log_ratio = np.log10(np.maximum(ratio, 1e-12))
        if np.isclose(dose, 0.0) or np.allclose(log_ratio, log_ratio[0]):
            rho, pval = np.nan, np.nan
        else:
            rho, pval = spearmanr(p['IC50'], log_ratio)
        corr_rows.append({'dose': dose, 'spearman_rho_IC50_vs_log10_relative_final_biomass': rho, 'p_value': pval})
        for i, species in enumerate(p['table']['species']):
            rows.append({
                'dose': dose, 'species': species, 'resistance': p['resistance'][i], 'IC50': p['IC50'][i],
                'final_biomass': r['N'][-1, i], 'control_final_biomass': control_final[i],
                'relative_final_biomass': ratio[i], 'log10_relative_final_biomass': log_ratio[i],
                'survived': np.isnan(r['extinction'][i]),
            })

    species_df = pd.DataFrame(rows)
    corr_df = pd.DataFrame(corr_rows)
    species_df.to_csv(out / 'resistance_species_outcomes.csv', index=False)
    corr_df.to_csv(out / 'resistance_selection_correlations.csv', index=False)

    nonzero = [d for d in sorted(runs) if d > 0]
    fig, axs = plt.subplots(1, len(nonzero), figsize=(4.2 * len(nonzero), 4), squeeze=False)
    for ax, dose in zip(axs.flat, nonzero):
        g = species_df[np.isclose(species_df['dose'], dose)]
        ax.scatter(g['IC50'], g['log10_relative_final_biomass'])
        ax.axhline(0, ls='--', lw=1)
        rho = corr_df.loc[np.isclose(corr_df['dose'], dose), 'spearman_rho_IC50_vs_log10_relative_final_biomass'].iloc[0]
        ax.set(xlabel='Species IC50', ylabel='log10(final biomass / control)', title=f'A0={dose:g}, Spearman rho={rho:.2f}')
    save_figure(fig, out / '01_IC50_vs_antibiotic_effect.png')

    fig, ax = plt.subplots(figsize=(8, 5))
    for dose in sorted(runs):
        p, r = runs[dose]
        ax.plot(r['t'], weighted_mean_trait(r['N'], p['IC50']), label=f'A0={dose:g}')
    ax.set(xlabel='Time', ylabel='Biomass-weighted mean IC50', title='Selection for resistance through time')
    ax.legend()
    save_figure(fig, out / '02_weighted_IC50_over_time.png')

    fig, ax = plt.subplots(figsize=(8, 5))
    final_means = []
    for dose in sorted(runs):
        p, r = runs[dose]
        final_means.append(weighted_mean_trait(r['N'], p['IC50'])[-1])
    ax.plot(sorted(runs), final_means, marker='o')
    ax.set(xlabel='Initial antibiotic A0', ylabel='Final biomass-weighted mean IC50', title='Resistance enrichment with dose')
    save_figure(fig, out / '03_final_resistance_vs_dose.png')

    print(corr_df.to_string(index=False))
    print(f'Wrote resistance-selection analysis to: {out.resolve()}')


if __name__ == '__main__':
    main()
