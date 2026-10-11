"""Sanity checks for the balanced antibiotic consumer-resource model."""
from dataclasses import replace
import numpy as np

from model import Config, generate, simulate, death_rates


def same_traits(a, b):
    cols = [
        'growth_score', 'intrinsic_mu_max', 'mu_max',
        'resistance_growth_multiplier', 'resistance_fitness_cost',
        'food_category_code', 'total_uptake_capacity',
        'M1_secretion_share', 'M2_secretion_share', 'initial_biomass',
        'leakage_fraction', 'resistance', 'IC50',
        'K_R0', 'K_M1', 'K_M2',
        'uptake_max_R0', 'uptake_max_M1', 'uptake_max_M2',
        'can_consume_R0', 'can_consume_M1', 'can_consume_M2',
    ]
    return np.allclose(a['table'][cols].to_numpy(float), b['table'][cols].to_numpy(float))


c = Config(seed=42)
p = generate(c)

# Balanced food categories: all seven are represented for n>=7 and counts differ by at most one.
counts = p['table']['food_category_code'].value_counts().reindex(range(1, 8), fill_value=0)
if c.n_species >= 7:
    assert np.all(counts.to_numpy() > 0)
assert counts.max() - counts.min() <= 1

# Resistance fitness-cost mapping.
rho = p['resistance']
expected_multiplier = 1.0 - c.resistance_cost_max * rho
assert np.allclose(p['growth_multiplier'], expected_multiplier)
assert np.allclose(p['table'].resistance_growth_multiplier.to_numpy(), expected_multiplier)
assert np.allclose(
    p['table'].mu_max.to_numpy(),
    p['table'].intrinsic_mu_max.to_numpy() * expected_multiplier,
)

# Capacity normalization: at saturating allowed resources, gross growth = effective mu_max.
assert np.allclose(np.sum(p['uptake'] * p['retain'], axis=1), p['table'].mu_max)

# Turning off the cost restores intrinsic mu_max without changing resistance/category identity.
p_no_cost = generate(replace(c, resistance_cost_max=0.0))
assert np.allclose(p_no_cost['table'].mu_max, p_no_cost['table'].intrinsic_mu_max)
assert np.allclose(p_no_cost['resistance'], p['resistance'])
assert np.array_equal(p_no_cost['food_category_code'], p['food_category_code'])

# Antibiotic is the only mortality source.
assert np.allclose(death_rates(0.0, p), 0.0)
d = death_rates(c.initial_antibiotic, p)
assert np.all(d >= 0)
assert np.all(d <= c.antibiotic_kmax + 1e-12)
# Higher IC50 should never increase killing at the same A.
order = np.argsort(p['IC50'])
assert np.all(np.diff(d[order]) <= 1e-12)

# Matched dose/decay conditions must keep the bacterial community fixed.
p_control = generate(replace(c, initial_antibiotic=0.0))
p_high = generate(replace(c, initial_antibiotic=2.0))
p_fast_decay = generate(replace(c, antibiotic_decay=0.1))
assert same_traits(p_control, p_high)
assert same_traits(p_control, p_fast_decay)

# Antibiotic follows the specified first-order decay exactly (within solver tolerance).
r = simulate(p)
expected_A = c.initial_antibiotic * np.exp(-c.antibiotic_decay * r['t'])
assert np.allclose(r['A'], expected_A, rtol=2e-7, atol=1e-10)

# Energy accounting and extinction handling.
assert np.max(np.abs(r['budget_error'])) < 1e-6
for i, te in enumerate(r['extinction']):
    if np.isfinite(te):
        assert np.all(r['N'][r['t'] > te, i] == 0)

# No-antibiotic control has zero antibiotic death and zero death-loss sink.
r0 = simulate(p_control)
assert np.allclose(r0['death_rates'], 0.0)
assert np.allclose(r0['death_loss'], 0.0, atol=1e-12)
assert np.max(np.abs(r0['budget_error'])) < 1e-6

print('All balanced antibiotic-model checks passed.')
print('Food-category counts:')
print(counts.to_string())
print(f'Default maximum resistance fitness cost: {c.resistance_cost_max:.2f}')
