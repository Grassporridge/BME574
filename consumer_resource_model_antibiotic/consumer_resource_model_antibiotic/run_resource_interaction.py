"""Factorial test of antibiotic dose x starting resource availability."""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model import Config
from antibiotic_analysis_common import add_common_cli, community_metrics, ensure_dir, run_condition, save_figure

DEFAULT_DOSES = [0.0, 0.25, 0.5, 1.0, 2.0]
DEFAULT_RESOURCES = [0.25, 0.5, 1.0, 1.5, 2.0]


def heatmap(ax, pivot, title, cbar_label):
    arr = pivot.to_numpy(dtype=float)
    im = ax.imshow(arr, aspect='auto', origin='lower', interpolation='nearest')
    ax.set(xticks=np.arange(len(pivot.columns)), xticklabels=[f'{x:g}' for x in pivot.columns],
           yticks=np.arange(len(pivot.index)), yticklabels=[f'{y:g}' for y in pivot.index],
           xlabel='Initial antibiotic A0', ylabel='Initial resource R0', title=title)
    return im, cbar_label


def main():
    parser = add_common_cli(argparse.ArgumentParser(description=__doc__), 'analysis_resource_interaction')
    parser.add_argument('--doses', type=float, nargs='+', default=DEFAULT_DOSES)
    parser.add_argument('--resources', type=float, nargs='+', default=DEFAULT_RESOURCES)
    args = parser.parse_args()
    out = ensure_dir(args.output)

    base = Config(seed=args.seed)
    rows = []
    for resource in args.resources:
        for dose in args.doses:
            p, r = run_condition(base, initial_resource=resource, initial_antibiotic=dose)
            row = community_metrics(p, r)
            row['dose'] = dose; row['initial_resource'] = resource
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(out / 'resource_antibiotic_factorial_summary.csv', index=False)

    metrics = [
        ('living_energy_fraction', 'Final living fraction', 'Living biomass / initial biomass+resource'),
        ('survivors', 'Species survival', 'Surviving species'),
        ('IC50_enrichment_ratio', 'Resistance enrichment', 'Final / initial weighted IC50'),
    ]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (metric, title, cbar_label) in zip(axs, metrics):
        pivot = df.pivot(index='initial_resource', columns='dose', values=metric).sort_index().sort_index(axis=1)
        im, label = heatmap(ax, pivot, title, cbar_label)
        fig.colorbar(im, ax=ax, label=label)
    save_figure(fig, out / '01_resource_by_antibiotic_heatmaps.png')

    fig, ax = plt.subplots(figsize=(8, 5))
    for resource, g in df.groupby('initial_resource'):
        g = g.sort_values('dose')
        ax.plot(g['dose'], g['living_energy_fraction'], marker='o', label=f'R0={resource:g}')
    ax.set(xlabel='Initial antibiotic A0', ylabel='Final living energy fraction', title='Does resource availability modify antibiotic impact?')
    ax.legend()
    save_figure(fig, out / '02_resource_dose_interaction.png')

    print(df[['initial_resource', 'dose', 'living_energy_fraction', 'survivors', 'final_weighted_IC50']].to_string(index=False))
    print(f'Wrote resource-interaction analysis to: {out.resolve()}')


if __name__ == '__main__':
    main()
