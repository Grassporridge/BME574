from dataclasses import dataclass, asdict
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.stats import truncnorm


@dataclass
class Config:
    # Legacy scenario flags kept so the existing analysis_tools.py still imports/runs.
    # starvation_mortality is no longer used in the ODE: death is antibiotic-only.
    variable_affinity: bool = True
    starvation_mortality: bool = False
    structured_community: bool = False

    affinity_log_sd: float = .7
    starvation_mean: float = .06
    starvation_sd: float = .015
    seed: int = 42
    n_species: int = 25
    t_end: float = 300.
    n_times: int = 1201
    growth_mean: float = .5
    growth_sd: float = .15
    growth_reference: float = 1.
    death_mean: float = .01
    death_sd: float = .003
    initial_biomass: float = .01
    initial_resource: float = 1.
    half_saturation: float = .1
    leakage: float = .2
    extinction_fraction: float = 1e-6

    # Antibiotic parameters
    initial_antibiotic: float = 1.0
    antibiotic_decay: float = 0.02
    antibiotic_kmax: float = 0.20
    antibiotic_hill: float = 2.0

    # Continuous resistance range, represented through species-specific IC50 values
    ic50_min: float = 0.05
    ic50_max: float = 2.0

    # Fitness cost of resistance. At resistance=1, maximum growth capacity is
    # multiplied by (1 - resistance_cost_max). Default 0.20 = 20% maximum cost.
    resistance_cost_max: float = 0.20


FOODS = ['R0', 'M1', 'M2']


def truncated(rng, mean, sd, low, high, size):
    return truncnorm.rvs(
        (low - mean) / sd,
        (high - mean) / sd,
        loc=mean,
        scale=sd,
        size=size,
        random_state=rng,
    )


def generate(c=Config()):
    if c.n_species < 3 or c.n_times < 2 or c.t_end <= 0:
        raise ValueError('Require >=3 species, >=2 output times, and positive duration')
    if not (0 <= c.leakage < 1 and 0 < c.extinction_fraction < 1):
        raise ValueError('Require 0<=leakage<1 and 0<extinction_fraction<1')
    if min(c.growth_sd, c.death_sd, c.growth_reference,
           c.initial_biomass, c.half_saturation) <= 0 or c.initial_resource < 0:
        raise ValueError('Invalid positive scale or initial resource')
    if c.initial_antibiotic < 0 or c.antibiotic_decay < 0 or c.antibiotic_kmax < 0:
        raise ValueError('Antibiotic concentration, decay, and kmax must be nonnegative')
    if c.antibiotic_hill <= 0:
        raise ValueError('antibiotic_hill must be positive')
    if c.ic50_min <= 0 or c.ic50_max < c.ic50_min:
        raise ValueError('Require 0 < ic50_min <= ic50_max')
    if not (0 <= c.resistance_cost_max < 1):
        raise ValueError('Require 0 <= resistance_cost_max < 1')

    rng = np.random.default_rng(c.seed)
    score = truncated(rng, c.growth_mean, c.growth_sd, 0, 1, c.n_species)

    # Preserve the original RNG sequence for community construction, but baseline
    # and starvation mortality are intentionally not used in the new model.
    _legacy_death_draw = truncated(rng, c.death_mean, c.death_sd, 0, np.inf, c.n_species)

    # Keep the old RNG consumption so the other random traits remain reproducible,
    # but replace the random food-category frequencies with a balanced design.
    _legacy_codes_draw = rng.integers(1, 8, c.n_species)

    # Seven nonempty food subsets encoded as bit masks:
    # 1=R0, 2=M1, 3=R0+M1, 4=M2, 5=R0+M2,
    # 6=M1+M2, 7=R0+M1+M2. Repeating 1..7 makes category
    # counts differ by at most one; a separate RNG shuffles which species
    # receives which category without changing the counts.
    codes = np.resize(np.arange(1, 8, dtype=int), c.n_species)
    food_rng = np.random.default_rng(c.seed + 300000)
    food_rng.shuffle(codes)
    mask = ((codes[:, None] >> np.arange(3)) & 1).astype(bool)
    weights = rng.dirichlet(np.ones(3), c.n_species) * mask
    weights /= weights.sum(axis=1, keepdims=True)

    # Retained energy fractions: R0 leaks, metabolites do not.
    retain = np.array([1 - c.leakage, 1., 1.])

    # Preliminary capacity; final capacity is recomputed below after applying
    # the resistance fitness cost.
    capacity = score * c.growth_reference / (weights @ retain)
    uptake = capacity[:, None] * weights
    split = rng.dirichlet([1., 1.], c.n_species)

    # Extra draws use an independent generator so original community draws stay fixed.
    extra = np.random.default_rng(c.seed + 100000)
    Kdraw = extra.lognormal(
        np.log(c.half_saturation),
        c.affinity_log_sd,
        (c.n_species, 3),
    )
    Kdraw = np.clip(Kdraw, .01, 1.)

    # Keep these legacy draws only to preserve the old RNG sequence for biomass/structure.
    _legacy_starve_draw = truncated(
        extra, c.starvation_mean, c.starvation_sd, 0, .15, c.n_species
    )
    biomass_draw = extra.lognormal(0, 1., c.n_species)

    # Use a separate RNG for resistance so adding antibiotics does not alter K/biomass draws.
    antibiotic_rng = np.random.default_rng(c.seed + 200000)
    resistance = antibiotic_rng.uniform(0, 1, c.n_species)
    IC50 = c.ic50_min * (c.ic50_max / c.ic50_min) ** resistance

    # Linear fitness cost of resistance. This affects growth even when A=0.
    # rho=0 -> no cost; rho=1 -> resistance_cost_max fractional reduction
    # in maximum resource-supported growth.
    resistance_fitness_cost = c.resistance_cost_max * resistance
    growth_multiplier = 1.0 - resistance_fitness_cost
    intrinsic_mu_max = score * c.growth_reference
    effective_mu_max = intrinsic_mu_max * growth_multiplier

    # Optional legacy structure for leakage/secretion and initial biomass.
    # Feeding categories remain balanced even when this flag is enabled.
    if c.structured_community:
        types = extra.choice([0, 1, 2, 3], c.n_species, p=[.25, .25, .25, .25])
        # Choose an R0 consumer to guarantee at least one R0 user secretes both metabolites.
        r0_users = np.flatnonzero(mask[:, 0])
        if len(r0_users):
            types[r0_users[0]] = 3
        leak = np.where(types == 0, 0., c.leakage)
        split = extra.dirichlet([1., 1.], c.n_species)
        split[types == 1] = [1., 0.]
        split[types == 2] = [0., 1.]
    else:
        leak = np.full(c.n_species, c.leakage)

    # Category codes remain the balanced codes generated above.
    codes = (mask * np.array([1, 2, 4])).sum(axis=1).astype(int)

    retain = np.column_stack([1 - leak, np.ones((c.n_species, 2))])

    # Normalize uptake capacity so growth at saturating permitted resources is
    # the resistance-cost-adjusted mu_max. Therefore resistance carries a
    # growth penalty even in the no-antibiotic control.
    capacity = effective_mu_max / np.sum(weights * retain, axis=1)
    uptake = capacity[:, None] * weights

    initial = (
        c.initial_biomass * biomass_draw / biomass_draw.sum()
        if c.structured_community
        else np.full(c.n_species, c.initial_biomass / c.n_species)
    )
    K = Kdraw if c.variable_affinity else np.full_like(Kdraw, c.half_saturation)

    # These legacy columns are retained only so the old analysis_tools.py does not break.
    # They are zero because mortality now comes only from antibiotic exposure.
    baseline_death = np.zeros(c.n_species)
    starvation_increment = np.zeros(c.n_species)

    table = pd.DataFrame({
        'species': [f'S{i+1:02d}' for i in range(c.n_species)],
        'growth_score': score,
        'intrinsic_mu_max': intrinsic_mu_max,
        'mu_max': effective_mu_max,
        'resistance_growth_multiplier': growth_multiplier,
        'resistance_fitness_cost': resistance_fitness_cost,
        'death_rate': baseline_death,
        'total_uptake_capacity': capacity,
        'M1_secretion_share': split[:, 0],
        'M2_secretion_share': split[:, 1],
        'initial_biomass': initial,
        'starvation_death_increment': starvation_increment,
        'leakage_fraction': leak,
        'food_category_code': codes,
        'resistance': resistance,
        'IC50': IC50,
    })


    category_names = {
        1: 'R0',
        2: 'M1',
        3: 'R0+M1',
        4: 'M2',
        5: 'R0+M2',
        6: 'M1+M2',
        7: 'R0+M1+M2',
    }
    table['food_category'] = [category_names[int(code)] for code in codes]

    for j, f in enumerate(FOODS):
        table[f'can_consume_{f}'] = mask[:, j]
        table[f'K_{f}'] = K[:, j]
        table[f'preference_{f}'] = weights[:, j]
        table[f'uptake_max_{f}'] = uptake[:, j]

    return {
        'config': c,
        'table': table,
        'uptake': uptake,
        'split': split,
        'retain': retain,
        'K': K,
        'leak': leak,
        'initial': initial,
        'starve': starvation_increment,  # legacy compatibility only
        'resistance': resistance,
        'IC50': IC50,
        'resistance_fitness_cost': resistance_fitness_cost,
        'growth_multiplier': growth_multiplier,
        'food_category_code': codes,
    }


def fluxes(N, R, p):
    u = p['uptake'] * np.maximum(R, 0) / (p['K'] + np.maximum(R, 0))
    gross = np.sum(u * p['retain'], axis=1)
    consumed = np.maximum(N, 0)[:, None] * u
    secretion = np.sum(
        consumed[:, 0, None] * p['leak'][:, None] * p['split'],
        axis=0,
    )
    return gross, consumed, secretion


def death_rates(A, p):
    """Species-specific antibiotic killing rate at antibiotic concentration A."""
    c = p['config']
    A = max(float(A), 0.0)
    IC50 = p['IC50']
    h = c.antibiotic_hill
    return c.antibiotic_kmax * A**h / (IC50**h + A**h)


def simulate(p):
    c = p['config']
    n = c.n_species
    N0 = p['initial']
    cutoff = c.initial_biomass / n * c.extinction_fraction

    # State layout:
    # [N_1 ... N_n, R0, M1, M2, A, cumulative_death_loss, cutoff_loss]
    y = np.r_[
        N0,
        c.initial_resource,
        0.,
        0.,
        c.initial_antibiotic,
        0.,
        0.,
    ]

    active = np.ones(n, dtype=bool)
    extinction = np.full(n, np.nan)
    times = np.linspace(0, c.t_end, c.n_times)
    states = np.empty((len(times), n + 6))
    t = 0.
    filled = np.zeros(len(times), dtype=bool)

    while t < c.t_end:
        def rhs(t, y):
            N = np.maximum(y[:n], 0) * active
            R = y[n:n+3]
            A = max(y[n+3], 0.0)

            gross, consumed, secreted = fluxes(N, R, p)
            death = death_rates(A, p)
            dA = -c.antibiotic_decay * A

            return np.r_[
                N * (gross - death),
                -consumed.sum(axis=0) + np.r_[0., secreted],
                dA,
                N @ death,
                0.,
            ]

        events = []
        ids = np.flatnonzero(active)
        for i in ids:
            def event(t, y, i=i):
                return y[i] - cutoff
            event.terminal = True
            event.direction = -1
            events.append(event)

        sol = solve_ivp(
            rhs,
            (t, c.t_end),
            y,
            method='DOP853',
            rtol=1e-8,
            atol=1e-13,
            dense_output=True,
            events=events or None,
            max_step=2.,
        )
        if not sol.success:
            raise RuntimeError(sol.message)

        end = sol.t[-1]
        sel = (times >= t) & (times <= end) & (~filled)
        if sel.any():
            states[sel] = sol.sol(times[sel]).T
            filled[sel] = True

        y = sol.y[:, -1].copy()
        t = end
        if t >= c.t_end:
            break

        hits = [ids[k] for k, e in enumerate(sol.t_events) if len(e)]
        hits = list(set(
            hits + list(np.flatnonzero(active & (y[:n] <= cutoff * (1 + 1e-8))))
        ))
        if not hits:
            raise RuntimeError('Extinction event without a crossing')

        y[-1] += y[hits].sum()
        y[hits] = 0
        active[hits] = False
        extinction[hits] = t

    assert filled.all()

    N = np.maximum(states[:, :n], 0)
    R = np.maximum(states[:, n:n+3], 0)
    A = np.maximum(states[:, n+3], 0)

    gross = np.array([
        fluxes(N[k], R[k], p)[0]
        for k in range(len(times))
    ])
    deaths = np.array([
        death_rates(A[k], p)
        for k in range(len(times))
    ])

    result = {
        't': times,
        'N': N,
        'R': R,
        'A': A,
        'gross': gross,
        'net': gross - deaths,
        'death_rates': deaths,
        'mortality': N * deaths,
        'extinction': extinction,
        'death_loss': states[:, -2],
        'cutoff_loss': states[:, -1],
        'cutoff': cutoff,
    }

    # Antibiotic is not included in this biomass/resource energy budget.
    budget = N.sum(axis=1) + R.sum(axis=1) + states[:, -2] + states[:, -1]
    result['budget_error'] = budget - (c.initial_biomass + c.initial_resource)

    assert np.max(np.abs(result['budget_error'])) < 1e-6
    assert states[:, :n+4].min() > -1e-9
    return result


def summarize(p, r):
    df = p['table'].copy()
    df['final_biomass'] = r['N'][-1]
    df['peak_biomass'] = r['N'].max(axis=0)
    df['biomass_time_integral'] = np.trapezoid(r['N'], r['t'], axis=0)
    df['extinction_time'] = r['extinction']
    df['survived'] = np.isnan(r['extinction'])
    df['final_percent'] = 100 * r['N'][-1] / max(r['N'][-1].sum(), 1e-300)
    return df


def export(p, r, out):
    from pathlib import Path
    import json

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    summarize(p, r).to_csv(out / 'species_summary.csv', index=False)

    frames = []
    for i, name in enumerate(p['table'].species):
        frames.append(pd.DataFrame({
            'time': r['t'],
            'species': name,
            'biomass': r['N'][:, i],
            'percent': 100 * r['N'][:, i] / np.maximum(r['N'].sum(axis=1), 1e-300),
            'gross_growth_rate': r['gross'][:, i],
            'death_rate': r['death_rates'][:, i],
            'net_growth_rate': r['net'][:, i],
            'biomass_death_flux': r['mortality'][:, i],
        }))
    pd.concat(frames).to_csv(out / 'species_timeseries.csv', index=False)

    pd.DataFrame({
        'time': r['t'],
        **{f: r['R'][:, j] for j, f in enumerate(FOODS)},
        'antibiotic': r['A'],
        'total_biomass': r['N'].sum(axis=1),
        'cumulative_death_loss': r['death_loss'],
        'cumulative_cutoff_loss': r['cutoff_loss'],
        'budget_error': r['budget_error'],
    }).to_csv(out / 'community_timeseries.csv', index=False)

    (out / 'config.json').write_text(json.dumps(asdict(p['config']), indent=2))
