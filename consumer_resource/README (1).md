# Batch bacterial consumer–resource simulation

This project simulates 25 species in a well-mixed, fixed-volume batch with one initial food R0 and two secreted metabolites M1 and M2. It is an exploratory Monod consumer–resource model inspired by MiCRM, not a calibrated bacterial culture prediction. Time, biomass, and chemical energy concentrations are normalized, not hours or measured concentrations.

## Run

```bash
python -m pip install -r requirements.txt
python run_simulation.py --seed 42 --replicates 20 --output outputs
jupyter lab analysis.ipynb
```

Run from this project folder. The notebook imports the accompanying scripts; keep them together. Run all notebook cells to reproduce baseline figures, repeat-community results, and cutoff sensitivity. CLI output and notebook output use the same seeds and equations. Set `--replicates 0` to skip the ensemble. Baseline outputs are included. Run `python check_model.py` for capacity, energy-accounting, analytic no-food mortality and irreversible-extinction checks. The delivered notebook cells were executed directly in a shared Python namespace and rich outputs captured; a Jupyter kernel was not available in the generation environment. Software versions used for the delivered run are in `outputs/environment.json`.

Files: `model.py` contains parameters, sampling, equations, event-driven extinction, and CSV exports; `plots.py` produces figures and repeated communities; `run_simulation.py` is the command-line entry point; `analysis.ipynb` walks through analysis; `outputs/` contains figures, data and configuration.

## Equations

Indices: i = 1,…,25 species; alpha = 0,1,2 for R0,M1,M2. N_i is living biomass concentration; R_alpha is extracellular chemical concentration in common energy-equivalent units.

Let P_i,alpha be food allocation weights, summing to 1 for every species, with zero on forbidden foods. Let s_i be the relative growth-capacity score and mu_i = s_i mu_ref the maximum gross growth rate when all permitted foods saturate uptake. Let l = 0.2 be the R0 leakage fraction, and f = (1-l,1,1) the retained fractions.

To ensure that maximum growth means the same thing for all food profiles:

\[
V_i=\frac{\mu_i}{\sum_\alpha P_{i\alpha}f_\alpha},\qquad
v_{i\alpha}=V_iP_{i\alpha}.
\]

V_i is total maximum uptake capacity; v_i,alpha is the maximum uptake rate per unit biomass for one food. Species with more leakage exposure need more uptake to attain the same maximum gross growth.

Monod uptake and gross growth:

\[
u_{i\alpha}(R_\alpha)=v_{i\alpha}\frac{R_\alpha}{K+R_\alpha},\qquad
\mu_i^{\mathrm{gross}}(t)=\sum_\alpha f_\alpha u_{i\alpha}(R_\alpha).
\]

K is the half-saturation concentration: uptake is half its maximum when R_alpha = K. Foods are substitutable; contributions add, without a requirement to obtain multiple foods together. There is no preferential switching, lag, intracellular storage, inhibition, or oxygen/pH limitation.

Living biomass:

\[
\frac{dN_i}{dt}=N_i[\mu_i^{\mathrm{gross}}(t)-d_i].
\]

d_i is constant per-capita mortality. Net growth is gross growth minus mortality, not a new random draw at each time. Biomass production flux is N_i mu_i^gross; mortality flux is N_i d_i. Death rates remain constant even when foods change. Nutrient limitation acts through growth.

Chemicals:

\[
\frac{dR_0}{dt}=-\sum_iN_iu_{i0},
\]
\[
\frac{dR_1}{dt}=-\sum_iN_iu_{i1}+l\sum_iN_iu_{i0}q_{i1},
\qquad
\frac{dR_2}{dt}=-\sum_iN_iu_{i2}+l\sum_iN_iu_{i0}q_{i2}.
\]

q_i1 + q_i2 = 1 divides leaked R0 energy between the two metabolites. All R0 consumers secrete both metabolites almost surely; secretion amounts differ. Species may consume their own products. Consuming metabolites produces biomass but no further modeled metabolite. There is no chemical input, renewal, abiotic decay, or dilution. Metabolites can still disappear through uptake. This feed-forward secretion structure is a simplifying choice, not a universal MiCRM requirement.

Energy accounting uses unit biomass yield per retained energy-equivalent unit. It does not imply that actual chemicals have equal mass yields or that all real uptake energy becomes biomass. A physical fit would require conversion factors, respiration/maintenance losses, and measured yields. We use explicit mortality and do not add a separate maintenance penalty.

For cumulative mortality L, dL/dt = sum_i d_i N_i. C is biomass removed by numerical extinction events. The checked accounting identity is:

\[
\sum_iN_i+R_0+R_1+R_2+L+C=N_{\mathrm{total}}(0)+R_0(0).
\]

Dead biomass is not returned to usable chemicals; L is an accounting sink. C records the tiny biomass discarded at cutoffs. This identity prevents artificial energy creation in cross-feeding.

## Parameters and distributions

Values below are exploratory choices agreed for this project, not experimentally estimated or universal defaults.

| Parameter | Meaning | Default and distribution | Allowed/model range |
|---|---|---|---|
| n_species | Number of species | 25 | Positive integer; >=3 recommended for guaranteed food coverage |
| seed | Reproducible random state | 42 | Nonnegative integer |
| t_end | End time | 300 | Positive normalized time |
| n_times | Output times | 1201, equally spaced | Integer >=2 |
| s_i | Relative maximum-growth score | Normal(mean 0.5, SD 0.15), truncated to [0,1] | 0–1 |
| mu_ref | Converts score to maximum growth | 1 | Positive inverse time |
| mu_i | Maximum gross growth | s_i mu_ref | 0–mu_ref |
| d_i | Per-capita death rate | Normal(mean 0.01, SD 0.003), truncated above 0 | Positive inverse time; no imposed upper bound |
| P_i,alpha | Food allocation weights | Dirichlet(1,1,1), masked then normalized | 0–1; row sum 1 |
| food permissions | Which foods can be consumed | Uniform among 7 nonempty subsets; first 3 species forced to R0-only, M1-only, M2-only | At least one food per species |
| V_i | Maximum total uptake | Derived from mu_i and food allocation | 0–mu_ref/(1-l), here up to 1.25 |
| v_i,alpha | Food-specific maximum uptake | V_i P_i,alpha | Nonnegative; zero if forbidden |
| K | Half-saturation concentration | Shared 0.1 | Positive |
| l | Fraction of R0 energy secreted | Shared 0.2 | 0 <= l < 1 |
| q_i1,q_i2 | Partition of leaked energy | Dirichlet(1,1); each marginal uniform on [0,1] | Nonnegative, sum 1 |
| retained energy yield | Biomass conversion in normalized units | 1 | Fixed model convention |
| total initial biomass | Starting living biomass | 0.01 | Positive |
| N_i(0) | Equal initial species biomass | 0.01/25 = 0.0004 | Positive |
| R0(0) | Initial starting food | 1 | Nonnegative |
| M1(0),M2(0) | Initial metabolites | 0,0 | Nonnegative |
| extinction_fraction | Cutoff relative to initial species biomass | 1e-6 | Positive, <1 recommended |
| N_cut | Absolute extinction cutoff | 4e-10 at baseline | Derived |
| replicates | Random communities | 20; seeds 42–61 | Nonnegative integer |

Truncation slightly changes the realized mean/SD from the underlying normal parameters. With only 25 draws, empirical means and spread vary. `species_summary.csv`, `realized_parameter_spread.csv`, and notebook histograms show the actual community, including permission groups. Assignment order has no biological meaning; the three guaranteed specialists are a stated exception to fully random permission sampling.

Growth and death are sampled independently. No growth-survival tradeoff is imposed. Death is deliberately much smaller than maximum growth. For a zero-growth population with d=0.01, survival fraction after 300 units is exp(-3), about 5%. R0-only, metabolite-only and mixed consumers can have different timing of growth; equal maximum capacities do not ensure equal actual growth.

## Extinction and numerical checks

At a downward crossing of N_cut, an event stops integration, removes that species permanently, records its extinction time and discarded biomass, and resumes integration. Abundances are not renormalized. Extinction means crossing this numerical convention, not a measured detection limit. A metabolite-dependent species removed early cannot recover later. Absolute cutoff avoids calling a growing species extinct merely because its percentage falls.

DOP853 integration uses relative tolerance 1e-8, absolute tolerance 1e-13, maximum step 2. Output spacing is 0.25. The code checks near-nonnegative states and accounting error below 1e-6. Sensitivity runs use cutoffs 1e-4, 1e-6, 1e-8 times starting species biomass, holding all other draws fixed.

## Outputs and interpretation

1. Species characterization: growth/death histograms, uptake matrix, secretion shares.
2. Batch dynamics: absolute species biomass on log scale, stacked composition percentages, chemicals, total biomass and richness. Extinct species are displayed at a plotting floor in log plots, not treated as positive simulation biomass.
3. Growth/death: gross and net per-capita growth (masked after extinction), community growth and mortality fluxes, cumulative mortality and cutoff losses. Constant death rates appear in the characterization and exported time series.
4. Trait–outcome plots: maximum growth, death, R0 preference and secretion versus final biomass. Extinct species are shown at a floor.
5. Effective diversity exp(Shannon) and energy-budget error. Diversity/percentages are only meaningful while total biomass is positive; if every species is extinct, do not interpret those zero-total displays.
6. Twenty repeated communities: descriptive associations with biomass-time integral. Pooled Spearman correlations are exploratory, not causal effects or independent-sample hypothesis tests. Species interact and rows within a community are dependent.

Final survival is a threshold/time-dependent binary measure, so also examine final biomass, peak biomass, extinction time, and biomass-time integral. A rising percentage does not necessarily mean rising biomass. All species can persist above the cutoff at t=300; this is a valid result, not evidence of perpetual coexistence. More species than food types can transiently coexist in a finite batch. With finite usable resources, positive mortality and no necromass recycling, this model does not sustain nonzero biomass forever.

To study causal trait effects, perturb one parameter while holding the same community fixed, then repeat over community seeds. Correlations alone cannot establish that a property causes survival.

## Sources

Marsland et al. (2020), *The Community Simulator: A Python package for microbial ecology*, PLOS ONE. https://doi.org/10.1371/journal.pone.0230430 — MiCRM, Monod response, closed-resource dynamics, consumption/secretion matrices and energy accounting.

Liao et al. (2020), *Modeling microbial cross-feeding at intermediate scale portrays community dynamics and species coexistence*, PLOS Computational Biology. https://doi.org/10.1371/journal.pcbi.1008135 — explicit resource/metabolite dynamics, Monod uptake and mortality.
