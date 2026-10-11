"""Analyze the tradeoff between antibiotic resistance and its growth fitness cost.

The intrinsic bacterial community is held fixed.  Only the maximum resistance
fitness cost is varied.  Each cost is evaluated both without antibiotic and at
one treatment dose so the ecological cost/benefit tradeoff can be seen directly.
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from model import Config
from antibiotic_analysis_common import (
    add_common_cli, assert_fixed_intrinsic_community, community_metrics,
    ensure_dir, run_condition, save_figure, species_outcomes,
)

DEFAULT_COSTS = [0.0, 0.10, 0.20, 0.30, 0.40]


def main():
    parser = add_common_cli(
        argparse.ArgumentParser(description=__doc__),
        'analysis_resistance_cost',
    )
    parser.add_argument('--costs', type=float, nargs='+', default=DEFAULT_COSTS,
                        help='Maximum fractional growth cost at resistance rho=1')
    parser.add_argument('--dose', type=float, default=1.0,
                        help='Antibiotic dose used for the treatment arm')
    args = parser.parse_args()
    out = ensure_dir(args.output)

    if any(c < 0 or c >= 1 for c in args.costs):
        raise ValueError('All resistance costs must satisfy 0 <= cost < 1.')

    base = Config(seed=args.seed)
    runs = {}
    rows = []
    species_rows = []
    intrinsic_ref = None

    for dose_label, dose in [('control', 0.0), ('antibiotic', args.dose)]:
        for cost in args.costs:
            p, r = run_condition(
                base,
                initial_antibiotic=dose,
                resistance_cost_max=cost,
            )
            if intrinsic_ref is None:
                intrinsic_ref = p
            else:
                assert_fixed_intrinsic_community(intrinsic_ref, p)

            runs[(dose_label, cost)] = (p, r)
            m = community_metrics(p, r)
            m.update({'arm': dose_label, 'dose': dose, 'resistance_cost_max': cost})
            rows.append(m)
            species_rows.append(
                species_outcomes(f'{dose_label}, cost={cost:g}', p, r)
                .assign(arm=dose_label, dose=dose, resistance_cost_max=cost)
            )

    summary = pd.DataFrame(rows).sort_values(['arm', 'resistance_cost_max'])
    species = pd.concat(species_rows, ignore_index=True)
    summary.to_csv(out / 'resistance_cost_summary.csv', index=False)
    species.to_csv(out / 'resistance_cost_species_outcomes.csv', index=False)

    # Community endpoints versus resistance cost.
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    outcomes = [
        ('final_biomass', 'Final biomass'),
        ('survivors', 'Surviving species'),
        ('final_shannon', 'Final Shannon diversity'),
        ('final_weighted_IC50', 'Final biomass-weighted IC50'),
    ]
    for ax, (column, ylabel) in zip(axs.flat, outcomes):
        for arm in ['control', 'antibiotic']:
            g = summary[summary['arm'] == arm].sort_values('resistance_cost_max')
            ax.plot(g['resistance_cost_max'], g[column], marker='o', label=arm)
        ax.set(xlabel='Maximum resistance fitness cost', ylabel=ylabel)
    axs[0, 0].legend()
    fig.suptitle(f'Resistance fitness-cost tradeoff (antibiotic arm A0={args.dose:g})')
    save_figure(fig, out / '01_resistance_cost_endpoints.png')

    # Show how the growth multiplier depends on resistance at each tested cost.
    rho = intrinsic_ref['resistance']
    order = np.argsort(rho)
    fig, ax = plt.subplots(figsize=(8, 5))
    for cost in args.costs:
        ax.plot(rho[order], 1.0 - cost * rho[order], marker='o', ms=3,
                label=f'cost={cost:g}')
    ax.set(xlabel='Resistance score (rho)', ylabel='Growth multiplier',
           title='Implemented resistance-growth tradeoff')
    ax.legend(fontsize=8)
    save_figure(fig, out / '02_growth_penalty_curves.png')

    # For each cost, compare antibiotic-treated final biomass with its matched control.
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for cost in args.costs:
        control = species[(species['arm'] == 'control') &
                          np.isclose(species['resistance_cost_max'], cost)]\
                          .set_index('species')['final_biomass']
        treated = species[(species['arm'] == 'antibiotic') &
                          np.isclose(species['resistance_cost_max'], cost)]\
                          .set_index('species')
        ratio = treated['final_biomass'] / np.maximum(control, 1e-300)
        order = np.argsort(treated['resistance'].to_numpy())
        ax.plot(treated['resistance'].to_numpy()[order],
                np.log10(np.maximum(ratio.to_numpy()[order], 1e-12)),
                marker='o', ms=3, label=f'cost={cost:g}')
    ax.axhline(0, ls='--', lw=1)
    ax.set(xlabel='Resistance score (rho)',
           ylabel='log10(treated final biomass / control)',
           title=f'Antibiotic benefit versus resistance under different fitness costs (A0={args.dose:g})')
    ax.legend(fontsize=8)
    save_figure(fig, out / '03_cost_changes_antibiotic_selection.png')

    print(summary[['arm', 'resistance_cost_max', 'final_biomass', 'survivors',
                   'final_shannon', 'final_weighted_IC50']].to_string(index=False))
    print(f'Wrote resistance-cost analysis to: {out.resolve()}')


if __name__ == '__main__':
    main()
