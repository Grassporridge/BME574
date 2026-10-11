"""Matched antibiotic dose-response analysis with one fixed bacterial community."""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model import Config
from antibiotic_analysis_common import (
    add_common_cli, assert_fixed_community, community_metrics, community_timeseries,
    ensure_dir, run_condition, save_figure, species_outcomes,
)

DEFAULT_DOSES = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0]


def main():
    parser = add_common_cli(argparse.ArgumentParser(description=__doc__), 'analysis_dose_response')
    parser.add_argument('--doses', type=float, nargs='+', default=DEFAULT_DOSES)
    parser.add_argument('--decay', type=float, default=0.02, help='Antibiotic decay rate')
    args = parser.parse_args()
    out = ensure_dir(args.output)

    base = Config(seed=args.seed, antibiotic_decay=args.decay)
    runs = []
    metrics = []
    species_rows = []
    community_rows = []
    reference_p = None

    for dose in args.doses:
        label = f'A0={dose:g}'
        p, r = run_condition(base, initial_antibiotic=dose)
        if reference_p is None:
            reference_p = p
        else:
            assert_fixed_community(reference_p, p)
        runs.append((dose, p, r))
        row = community_metrics(p, r); row['dose'] = dose; metrics.append(row)
        species_rows.append(species_outcomes(label, p, r).assign(dose=dose))
        community_rows.append(community_timeseries(label, p, r).assign(dose=dose))

    metrics_df = pd.DataFrame(metrics).sort_values('dose')
    species_df = pd.concat(species_rows, ignore_index=True)
    community_df = pd.concat(community_rows, ignore_index=True)
    metrics_df.to_csv(out / 'dose_summary.csv', index=False)
    species_df.to_csv(out / 'dose_species_outcomes.csv', index=False)
    community_df.to_csv(out / 'dose_community_timeseries.csv', index=False)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for dose, p, r in runs:
        axs[0].plot(r['t'], r['N'].sum(axis=1), label=f'A0={dose:g}')
        axs[1].plot(r['t'], r['A'], label=f'A0={dose:g}')
    axs[0].set(xlabel='Time', ylabel='Total living biomass', title='Population response to antibiotic dose')
    axs[1].set(xlabel='Time', ylabel='Antibiotic concentration', title='Antibiotic exposure histories')
    axs[0].legend(fontsize=8); axs[1].legend(fontsize=8)
    save_figure(fig, out / '01_dose_timecourses.png')

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    x = metrics_df['dose']
    axs[0, 0].plot(x, metrics_df['final_biomass'], marker='o')
    axs[0, 0].set(ylabel='Final biomass')
    axs[0, 1].plot(x, metrics_df['survivors'], marker='o')
    axs[0, 1].set(ylabel='Surviving species')
    axs[1, 0].plot(x, metrics_df['final_shannon'], marker='o')
    axs[1, 0].set(ylabel='Final Shannon diversity')
    axs[1, 1].plot(x, metrics_df['final_weighted_IC50'], marker='o')
    axs[1, 1].set(ylabel='Final biomass-weighted IC50')
    for ax in axs.flat:
        ax.set(xlabel='Initial antibiotic A0')
    save_figure(fig, out / '02_dose_endpoints.png')

    # Heatmap: final biomass relative to each species' no-antibiotic control outcome.
    control = species_df[np.isclose(species_df['dose'], 0.0)].set_index('species')['final_biomass']
    ordered_species = reference_p['table'].sort_values('IC50')['species'].tolist()
    heat = []
    nonzero_doses = [d for d in args.doses if d > 0]
    for d in nonzero_doses:
        g = species_df[np.isclose(species_df['dose'], d)].set_index('species')['final_biomass']
        ratio = (g / control).reindex(ordered_species)
        heat.append(np.log10(np.maximum(ratio.to_numpy(), 1e-12)))
    if heat:
        fig, ax = plt.subplots(figsize=(12, max(3.5, 0.55 * len(nonzero_doses) + 2)))
        im = ax.imshow(np.asarray(heat), aspect='auto', interpolation='nearest')
        ax.set(yticks=np.arange(len(nonzero_doses)), yticklabels=[f'A0={d:g}' for d in nonzero_doses],
               xticks=np.arange(len(ordered_species)), xticklabels=ordered_species,
               xlabel='Species, ordered from low to high IC50', ylabel='Dose',
               title='Species final biomass relative to no-antibiotic control (log10 ratio)')
        ax.tick_params(axis='x', rotation=90, labelsize=7)
        fig.colorbar(im, ax=ax, label='log10(treatment/control final biomass)')
        save_figure(fig, out / '03_species_dose_heatmap.png')

    print(metrics_df[['dose', 'final_biomass', 'survivors', 'final_shannon', 'final_weighted_IC50']].to_string(index=False))
    print(f'Wrote dose-response analysis to: {out.resolve()}')


if __name__ == '__main__':
    main()
