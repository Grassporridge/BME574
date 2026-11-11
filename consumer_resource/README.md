# Baseline multispecies consumer-resource model (no antibiotics)

A modular, matrix-based ODE model of 25 microbial species competing for and cross-feeding on 3 shared resources. This is the **no-antibiotic baseline**. Susceptibility traits are generated and saved but are *not used* by the dynamics.

## Files

| File | Purpose |
|---|---|
| `consumer_resource_model.py` | Dataclasses, Monod growth, uptake, cross-feeding, ODE right-hand side, `simulate()`, diagnostics, JSON I/O |
| `generate_community.py` | Seeded generator for the 25-species, 3-resource community; writes JSON + CSV tables |
| `run_baseline.py` | Loads or generates a community, runs the baseline, saves the time course, prints diagnostics |
| `01_baseline_analysis.ipynb` | Plots and checks (abundances, resources, biomass, richness, Shannon, matrices, extinctions, nonnegativity, steady state) |
| `test_consumer_resource_model.py` | pytest tests |

## Quick start

```bash
pip install numpy scipy pandas matplotlib pytest jupyter
python generate_community.py --seed 42      # writes community.json and community_*.csv
python run_baseline.py                       # writes baseline_timecourse.csv, prints diagnostics
pytest -q
jupyter notebook 01_baseline_analysis.ipynb
```

`run_baseline.py` options: `--community PATH`, `--seed N`, `--regenerate`, `--t-end T`, `--n-points K`, `--out PATH`.

## Model

Indices: species *i* = 1..S (S = 25), resources *j*, *k* = 1..M (M = 3).

```
Monod_ij       = R_j / (K_ij + R_j)
growth_i       = sum_j mu_max_i * U_ij * Monod_ij
uptake_ij      = (mu_max_i * U_ij / e_ij) * Monod_ij          (per unit biomass)

dN_i/dt = N_i * ( growth_i - mortality_i - self_limitation_i * N_i )

dR_j/dt = supply_j
          - sum_i N_i * uptake_ij                              (consumption)
          + sum_i N_i * leakage_i * P_ij * sum_k uptake_ik     (cross-feeding)
          - decay_j * R_j
```

Because each row of `U` sums to 1, `mu_max_i` is the species' true maximum growth rate when all resources are saturating. The yield `e_ij` converts uptake into biomass: `e_ij * uptake_ij = mu_max_i * U_ij * Monod_ij`.

**Mass bookkeeping.** Of each unit of resource taken up by species *i*: `e_ij` becomes biomass, `leakage_i` is released as metabolites (split across resources by `P_i·`), and the remainder is lost (respiration/waste). `validate_params` requires `max_j e_ij + leakage_i <= 1` so the model never creates mass.

**Cross-feeding** happens because metabolites released by one species enter the shared resource pools that other species consume. There are no explicit species-to-species links. `crossfeeding_potential(p)` returns a derived S×S matrix `A[i,k] = leakage_i * sum_m P_im U_km` for inspection.

**Nonnegativity.** Growth and loss terms are proportional to `N_i`, so `N_i = 0` is invariant. Monod terms clip `R` at 0 to protect against solver round-off. Stiff-capable `LSODA` with tight tolerances is the default.

## Parameters

### Species vectors (shape S)

| Parameter | Meaning | Default generation (placeholder) |
|---|---|---|
| `mu_max` | Max growth rate, 1/time | Uniform(0.4, 1.0) |
| `mortality` | Density-independent loss, 1/time | Uniform(0.05, 0.15) |
| `self_limitation` | Density-dependent loss coefficient (`s_i` in `-s_i N_i^2`) | 0.02 for all |
| `leakage` | Fraction of uptake released as byproducts | Uniform(0.1, 0.3) |
| `initial_abundance` | Starting biomass | 10^Uniform(-2, -1) |
| `susceptibility` | **Future trait, unused**; 0 = strongly resistant, 1 = highly susceptible | Beta(2, 2) |

### Species × resource matrices (shape S×M)

| Parameter | Meaning | Default generation (placeholder) |
|---|---|---|
| `preference` (U) | Resource-use preference; rows sum to 1 | Dirichlet(0.6, 0.6, 0.6): mostly specialists |
| `half_saturation` (K) | Monod half-saturation constant, resource units | 0.5 × LogNormal(0, 0.5) |
| `efficiency` (e) | Biomass yield per unit resource | Uniform(0.2, 0.5) |
| `production` (P) | Which metabolites species *i* releases; rows sum to 1 | Dirichlet(0.5, 0.5, 0.5) |

### Resource vectors (shape M)

| Parameter | Meaning | Default |
|---|---|---|
| `supply` | External input rate | (2.0, 0.2, 0.2): R1 is the primary nutrient, R2/R3 are mainly cross-fed |
| `decay` | First-order removal rate | 0.05 each |
| `initial_resources` | Starting concentration | (1.0, 0.1, 0.1) |

### Simulation (`SimulationConfig`)

`t_end` (2000), `n_points` (2001), `method` (LSODA), `rtol` (1e-8), `atol` (1e-12), `extinction_threshold` (1e-6; used only by diagnostics such as richness and Shannon diversity, never by the solver).

To change the resource-use or cross-feeding structure, edit `preference` and `production` (or pass different arguments to `generate_community`); the equations never need to change.

## Reference baseline result (seed 42)

With the defaults, 8 of 25 species persist at `t = 2000`, with a Shannon index of about 1.77 and total biomass near 9.5. Residual derivatives are around 1e-6, so the system is effectively at a steady state. Some marginal species decline very slowly toward extinction, which is why the default horizon is long. The loss of most species is expected: with only 3 resources, competitive exclusion limits stable coexistence unless cross-feeding and self-limitation open extra niches. If you want higher richness, try larger `self_limitation`, more resources, or more distinct supply and cross-feeding structure. These are not parameter sweeps built into this version, only knobs.

## Planned antibiotic extension

Susceptibility is stored in `FutureTraits`, deliberately separate from `EcologicalParams`. A later version can add an `AntibioticParams` dataclass and an antibiotic concentration state `A`, then add a kill term such as `- kill_max * susceptibility_i * A / (A50 + A)` to `dN_i/dt`, with an antibiotic degradation/decay equation for `A`. Resistant traits (e.g. a degradation contribution per species) can reuse the same per-species vector structure. Because the dynamics read only `EcologicalParams`, existing code and the saved `community.json` remain valid. Note the sign convention: **0 = resistant, 1 = susceptible**.

## Assumptions that are placeholders

Everything numeric above is a toy default and should be replaced with biological data when available:

- All distributions for growth rates, mortality, K, yields, leakage, preferences and production profiles
- The value of `self_limitation` (a stand-in for unmodeled crowding, spatial structure, phage, or other limits; set to 0 to remove)
- The number of resources and which are externally supplied versus cross-fed
- Linear first-order resource decay and constant supply (chemostat-like, no biomass dilution)
- Leakage as a fixed fraction of uptake, with production independent of which resource was consumed
- Additive, non-substitutable growth on multiple resources (no diauxie or regulation)
- The Beta(2, 2) susceptibility distribution and the independence of susceptibility from all other traits (real traits may be correlated, e.g. growth-rate trade-offs)
- Arbitrary time and biomass units
- The 1e-6 extinction threshold (no demographic stochasticity or immigration)
