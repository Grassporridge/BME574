"""Shared utilities for matched antibiotic analyses.

All analyses keep the random community fixed with the same seed unless a script
explicitly varies a biological parameter.  This lets treatment differences be
attributed to the parameter being tested rather than to a newly drawn community.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from model import Config, generate, simulate, summarize


def ensure_dir(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_figure(fig, path, dpi=150):
    """Save a figure using an absolute Windows-safe path."""
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(str(path), dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    return path


def shannon_from_biomass(N):
    N = np.asarray(N, dtype=float)
    total = N.sum(axis=-1)
    shares = np.divide(N, total[..., None], out=np.zeros_like(N), where=total[..., None] > 0)
    H = -np.sum(shares * np.log(np.maximum(shares, 1e-300)), axis=-1)
    return np.where(total > 0, H, np.nan)


def weighted_mean_trait(N, trait):
    N = np.asarray(N, dtype=float)
    trait = np.asarray(trait, dtype=float)
    total = N.sum(axis=1)
    return np.divide(N @ trait, total, out=np.full(len(total), np.nan), where=total > 0)


def community_metrics(p, r):
    total = r['N'].sum(axis=1)
    H = shannon_from_biomass(r['N'])
    ic50_mean = weighted_mean_trait(r['N'], p['IC50'])
    initial_total_energy = p['config'].initial_biomass + p['config'].initial_resource
    return {
        'initial_antibiotic': p['config'].initial_antibiotic,
        'antibiotic_decay': p['config'].antibiotic_decay,
        'antibiotic_kmax': p['config'].antibiotic_kmax,
        'antibiotic_hill': p['config'].antibiotic_hill,
        'resistance_cost_max': p['config'].resistance_cost_max,
        'initial_resource': p['config'].initial_resource,
        'final_antibiotic': float(r['A'][-1]),
        'antibiotic_auc': float(np.trapezoid(r['A'], r['t'])),
        'final_biomass': float(total[-1]),
        'living_energy_fraction': float(total[-1] / initial_total_energy),
        'minimum_biomass': float(total.min()),
        'time_of_minimum_biomass': float(r['t'][np.argmin(total)]),
        'biomass_auc': float(np.trapezoid(total, r['t'])),
        'survivors': int(np.sum(np.isnan(r['extinction']))),
        'final_shannon': float(H[-1]),
        'initial_weighted_IC50': float(ic50_mean[0]),
        'final_weighted_IC50': float(ic50_mean[-1]),
        'IC50_enrichment_ratio': float(ic50_mean[-1] / ic50_mean[0]) if ic50_mean[0] > 0 else np.nan,
        'cumulative_antibiotic_death_loss': float(r['death_loss'][-1]),
        'max_budget_error': float(np.max(np.abs(r['budget_error']))),
    }


def run_condition(base_config: Config, **changes):
    config = replace(base_config, **changes)
    p = generate(config)
    r = simulate(p)
    return p, r


def community_timeseries(label, p, r):
    total = r['N'].sum(axis=1)
    H = shannon_from_biomass(r['N'])
    mean_ic50 = weighted_mean_trait(r['N'], p['IC50'])
    return pd.DataFrame({
        'condition': label,
        'time': r['t'],
        'total_biomass': total,
        'R0': r['R'][:, 0],
        'M1': r['R'][:, 1],
        'M2': r['R'][:, 2],
        'antibiotic': r['A'],
        'shannon': H,
        'weighted_mean_IC50': mean_ic50,
        'surviving_species': (r['N'] > 0).sum(axis=1),
        'cumulative_antibiotic_death_loss': r['death_loss'],
    })


def species_outcomes(label, p, r):
    df = summarize(p, r).copy()
    df.insert(0, 'condition', label)
    return df


def representative_species(p):
    """Return low-, median-, and high-resistance species indices."""
    ic50 = np.asarray(p['IC50'])
    order = np.argsort(ic50)
    idx = {
        'low resistance': int(order[0]),
        'median resistance': int(order[len(order) // 2]),
        'high resistance': int(order[-1]),
    }
    return idx


def assert_fixed_community(reference_p, other_p, include_ic50=True):
    """Raise if matched analyses accidentally redraw species traits."""
    cols = [
        'growth_score', 'intrinsic_mu_max', 'mu_max',
        'resistance_growth_multiplier', 'resistance_fitness_cost',
        'food_category_code', 'total_uptake_capacity',
        'M1_secretion_share', 'M2_secretion_share', 'initial_biomass',
        'leakage_fraction', 'K_R0', 'K_M1', 'K_M2',
        'uptake_max_R0', 'uptake_max_M1', 'uptake_max_M2',
        'can_consume_R0', 'can_consume_M1', 'can_consume_M2',
        'resistance',
    ]
    if include_ic50:
        cols.append('IC50')
    a = reference_p['table'][cols].to_numpy(dtype=float)
    b = other_p['table'][cols].to_numpy(dtype=float)
    if not np.allclose(a, b, equal_nan=True):
        raise AssertionError('Community traits changed between matched conditions.')


def assert_fixed_intrinsic_community(reference_p, other_p, include_ic50=True):
    """Check traits that should stay fixed when resistance cost itself is varied."""
    cols = [
        'growth_score', 'intrinsic_mu_max',
        'M1_secretion_share', 'M2_secretion_share', 'initial_biomass',
        'leakage_fraction', 'food_category_code',
        'K_R0', 'K_M1', 'K_M2',
        'preference_R0', 'preference_M1', 'preference_M2',
        'can_consume_R0', 'can_consume_M1', 'can_consume_M2',
        'resistance',
    ]
    if include_ic50:
        cols.append('IC50')
    a = reference_p['table'][cols].to_numpy(dtype=float)
    b = other_p['table'][cols].to_numpy(dtype=float)
    if not np.allclose(a, b, equal_nan=True):
        raise AssertionError('Intrinsic community traits changed between matched conditions.')


def write_config(config, path):
    from dataclasses import asdict
    Path(path).write_text(json.dumps(asdict(config), indent=2))


def add_common_cli(parser, default_output):
    parser.add_argument('--seed', type=int, default=42, help='Random community seed (default: 42)')
    parser.add_argument('--output', default=default_output, help='Output directory')
    return parser
