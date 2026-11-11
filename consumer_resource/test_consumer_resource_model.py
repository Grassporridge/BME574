"""Tests for the baseline consumer-resource model. Run with:  pytest -q"""

import numpy as np

from consumer_resource_model import (SimulationConfig, community_from_dict, community_to_dict,
                                     crossfeeding_potential, rhs, shannon_diversity, richness,
                                     simulate, validate_params)
from generate_community import N_RESOURCES, N_SPECIES, generate_community

SEED = 42


def make_community():
    return generate_community(seed=SEED)


def initial_state(c):
    p = c.ecology
    return np.concatenate([p.initial_abundance, p.initial_resources])


def test_rhs_output_dimensions():
    c = make_community()
    dy = rhs(0.0, initial_state(c), c.ecology)
    assert dy.shape == (N_SPECIES + N_RESOURCES,)
    assert np.all(np.isfinite(dy))


def test_initial_states_nonnegative():
    c = make_community()
    assert np.all(initial_state(c) >= 0)
    assert np.all(c.ecology.initial_abundance > 0)


def test_matrix_shapes():
    p = make_community().ecology
    S, M = N_SPECIES, N_RESOURCES
    for name in ("preference", "half_saturation", "efficiency", "production"):
        assert getattr(p, name).shape == (S, M), name
    for name in ("mu_max", "mortality", "self_limitation", "leakage", "initial_abundance"):
        assert getattr(p, name).shape == (S,), name
    for name in ("supply", "decay", "initial_resources"):
        assert getattr(p, name).shape == (M,), name
    assert crossfeeding_potential(p).shape == (S, S)
    assert np.allclose(p.preference.sum(axis=1), 1.0)
    assert np.allclose(p.production.sum(axis=1), 1.0)


def test_susceptibility_in_unit_interval():
    c = make_community()
    s = c.future.susceptibility
    assert s.shape == (N_SPECIES,)
    assert np.all(s >= 0.0) and np.all(s <= 1.0)


def test_generation_is_reproducible_and_seed_dependent():
    a, b, d = generate_community(seed=1), generate_community(seed=1), generate_community(seed=2)
    assert np.array_equal(a.ecology.mu_max, b.ecology.mu_max)
    assert np.array_equal(a.future.susceptibility, b.future.susceptibility)
    assert not np.array_equal(a.ecology.mu_max, d.ecology.mu_max)


def test_json_roundtrip():
    c = make_community()
    c2 = community_from_dict(community_to_dict(c))
    assert np.allclose(c.ecology.preference, c2.ecology.preference)
    assert np.allclose(c.future.susceptibility, c2.future.susceptibility)
    validate_params(c2)


def test_susceptibility_is_not_used_by_dynamics():
    """Changing susceptibility must not change the baseline simulation."""
    c = make_community()
    cfg = SimulationConfig(t_end=20.0, n_points=21)
    r1 = simulate(c, cfg)
    c.future.susceptibility = np.ones(N_SPECIES)
    r2 = simulate(c, cfg)
    assert np.array_equal(r1.N, r2.N)


def test_short_simulation_runs():
    c = make_community()
    res = simulate(c, SimulationConfig(t_end=50.0, n_points=51))
    assert res.success
    assert res.N.shape == (51, N_SPECIES)
    assert res.R.shape == (51, N_RESOURCES)
    assert np.all(np.isfinite(res.N)) and np.all(np.isfinite(res.R))


def test_no_state_substantially_negative():
    c = make_community()
    res = simulate(c, SimulationConfig(t_end=300.0, n_points=301))
    assert res.N.min() > -1e-8
    assert res.R.min() > -1e-8


def test_diversity_helpers():
    N = np.array([[1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    assert list(richness(N)) == [2, 1]
    assert np.allclose(shannon_diversity(N), [np.log(2), 0.0])
