# BME574 antibiotic–resistance community project

This repository contains the shared model, sweep infrastructure, notebook, and result-storage
conventions for the BME574 group project. The project studies how beta-lactam antibiotic treatment
and beta-lactamase (Bla) inhibition reshape a community containing susceptible and resistant
populations.

The model is based on:

- Ma HR, Xu HZ, Kim K, Anderson DJ, You L. “Private benefit of beta-lactamase dictates selection
  dynamics of combination antibiotic treatment.” *Nature Communications* 15, 8337 (2024).
- The accompanying supplementary information, especially Sections 1–3.
- The project proposal and midterm-presentation work division supplied by the group.

The local source PDFs used during model curation are outside this cloned repository in the parent
project directory. The repository itself contains the implementation and generated/shared tables.

## Project question and general idea

The main question is:

> How do antibiotic exposure, inhibitor action, resistance burden/private benefit, Bla-mediated
> antibiotic degradation, nutrient recycling, and initial susceptible/resistant composition combine
> to determine final community composition, diversity, recovery, and selection for resistance?

The workflow is:

~~~text
antibiotic + inhibitor treatment
        ↓
drug and free-Bla dynamics
        ↓
population-specific growth and lysis
        ↓
final composition, diversity, recovery, and resistance selection
~~~

The base model has two live populations:

- 'n_s': susceptible population density;
- 'n_r': resistant Bla-producing population density.

It also tracks:

- 's': shared nutrient;
- 'a': antibiotic concentration;
- 'b': free extracellular Bla.

The model is dimensionless, so numerical values should be interpreted according to the scaling in
Ma et al., not as direct clinical concentrations without an additional dimensional mapping.

## Base model equations

The five-state ODE system is:

\[
\frac{dn_s}{d\tau}=(g-l)n_s
\]

\[
\frac{dn_r}{d\tau}=(\alpha g-\beta l)n_r
\]

\[
\frac{ds}{d\tau}=(\xi l-g)n_s+(\xi\beta l-\alpha g)n_r
\]

\[
\frac{da}{d\tau}=-\kappa_bba-\phi n_ra-d_aa
\]

\[
\frac{db}{d\tau}=\beta l n_r-d_b\iota b.
\]

The closure relations are:

\[
g=\frac{s}{1+s}
\]

\[
l=\gamma\frac{a^{h_a}}{1+a^{h_a}}g
\]

\[
\iota=\frac{i^{h_i}}{1+i^{h_i}}
\]

\[
\beta=\beta_{\min}+c(1-\beta_{\min})\iota
\]

\[
\phi=\phi_{\max}(1-c\iota).
\]

Here:

- 'alpha' is the resistant growth-rate factor; '1-alpha' is the burden of resistance;
- 'beta' is the resistant lysis-rate factor; '1-beta' is the private Bla benefit;
- 'gamma' and 'h_a' determine antibiotic-driven lysis;
- 'i' and 'h_i' determine inhibitor activity 'iota';
- 'c' controls inhibitor effectiveness against intracellular/private and extracellular/public Bla;
- 'kappa_b' controls antibiotic degradation by free Bla;
- 'phi_max' controls antibiotic degradation by living resistant cells;
- 'd_b' controls inhibitor-dependent free-Bla inactivation;
- 'xi' controls nutrient recycling;
- 'd_a' is basal antibiotic decay.

The paper’s general selection criterion is:

\[
1-\beta>\frac{1-\alpha}{l/g}.
\]

Under saturating antibiotic and inhibitor concentrations, the simplified criterion is:

\[
(1-c)(1-\beta_{\min})>\frac{1-\alpha}{\gamma}.
\]

These criteria are interpretation tools. Dynamic simulations are still required because the
populations continuously change the antibiotic, Bla, and nutrient environment.

## Shared baseline values and parameter ranges

The base values and randomized ranges below come from the supplementary information. A parameter
should not be silently changed in one person’s analysis: record deviations in the relevant report
and preserve the shared baseline for comparison.

| Parameter | Base value | Range or role |
|---|---:|---|
| 'alpha' | 0.95 | 0.75–1.00 |
| 'beta_min' | 0.90 | 0–1 |
| 'c' | 0.70 | 0–1 |
| 'gamma' | 1.38 | 1.1–1.4 |
| 'h_a' | 3 | 1–5 |
| 'h_i' | 2 | 1–5 |
| 'kappa_b' | 0.35 | 0–1 |
| 'd_a' | 0.02 | fixed in the base model |
| 'd_b' | 1 | 1–10 |
| 'xi' | 0.8 | 0–1 |
| 'phi_max' | 1 | 0–5 |
| 'a0' | treatment-dependent | paper’s base dose example: 2.15; dose range 1–100 |
| 'i' | treatment-dependent | inhibitor dose range 0.1–10 |

Default initial conditions are:

\[
n_s(0)=0.2,\quad n_r(0)=0.2,\quad s(0)=4,\quad b(0)=0,
\]

with 'a(0)=a0'.

## Person 1–Person 4 task schema

The four workstreams share the same base model and endpoints. Each person owns the primary
parameters below while initially holding the other parameter blocks at their shared baseline.

### Person 1 — resistance and selection parameters

**Parameters:** 'alpha', 'beta_min', 'c'

**Purpose:** test the Ma selection criterion and determine how resistance burden, private benefit,
and inhibitor efficacy change selection for susceptible versus resistant cells.

**Core tasks:**

- sweep the three criterion parameters over the supplementary-information ranges;
- compare dynamic outcomes with the general and simplified selection criteria;
- quantify final resistant fraction, diversity, population survival, extinction, recovery, and
  bottleneck depth;
- identify parameter regions that favor susceptible cells, resistant cells, or coexistence;
- document criterion agreement and mismatches.

### Person 2 — antibiotic pharmacodynamics

**Parameters:** 'gamma', 'h_a', 'a0'

**Purpose:** determine how strongly and how steeply the antibiotic produces lysis.

**Core tasks:**

- sweep maximum lysis coefficient, antibiotic Hill coefficient, and antibiotic dose;
- distinguish dose effects from lysis-shape effects;
- generate antibiotic-response maps and representative time courses;
- provide common antibiotic settings for the other workstreams;
- report total biomass loss, final composition, peak bottleneck, recovery, and time to steady state.

### Person 3 — inhibitor pharmacodynamics

**Parameters:** 'h_i', 'i', 'd_b'

**Purpose:** determine how inhibitor dose, inhibitor response steepness, and free-Bla inactivation
alter private protection and public antibiotic degradation.

**Core tasks:**

- sweep '(i, h_i)', '(i, d_b)', and justified '(h_i, d_b)' response planes;
- track 'iota', 'beta', 'phi', free Bla, antibiotic concentration, and population trajectories;
- calculate final resistant fraction, Shannon diversity, total biomass, steady-state time,
  susceptible recovery, resistant extinction, and peak population loss;
- group representative trajectories by observed endpoint behavior;
- test selected inhibitor regimes across initial susceptible/resistant fractions;
- explicitly account for the fact that increasing 'h_i' has opposite effects below and above the
  scaled transition dose 'i=1'.

### Person 4 — degradation, nutrient dynamics, and observation/inference

**Parameters:** 'kappa_b', 'd_a', 'phi_max', 'xi'

**Purpose:** determine how quickly antibiotic clears and how nutrient recycling affects population
recovery and ecological outcomes.

**Core tasks:**

- sweep free-Bla antibiotic degradation, basal antibiotic decay, resistant-cell degradation, and
  nutrient recycling;
- quantify crash-and-recovery versus crash-and-die behavior;
- evaluate robustness of the other workstreams’ conclusions to degradation and nutrient settings;
- add synthetic observation noise and sampling schedules where appropriate;
- assess parameter sensitivity, identifiability, and model distinguishability;
- store inference reports under 'results/reports/'.

## Shared endpoints

All workstreams should use the same definitions unless a report explicitly justifies an alternative:

- final resistant fraction:
  \[
  f_R=\frac{n_r}{n_s+n_r};
  \]
- two-population Shannon diversity:
  \[
  H=-f_S\ln(f_S)-f_R\ln(f_R);
  \]
- final total live biomass 'n_s+n_r';
- time to numerical steady state;
- resistant extinction and susceptible survival flags;
- susceptible recovery time relative to a no-treatment reference;
- minimum biomass and peak fractional population loss;
- initial composition and final composition when the analysis concerns the project’s central
  susceptible/resistant-ratio question.

Thresholds, simulation horizons, solver tolerances, and endpoint conventions belong in the relevant
analysis module and report. Do not infer extinction merely from a population being outnumbered.

## Repository map

~~~text
BME574/
├── README.md
├── LICENSE
├── notebooks/
│   └── 01_sweep_execution.ipynb
├── src/
│   ├── model/
│   │   └── ma2024_model.py
│   ├── analysis/
│   │   ├── sweep_engine.py
│   │   └── grid_runner.py
│   └── visualization/
│       └── plot_style.py
└── results/
    ├── tables/
    │   ├── treatment_grid.csv
    │   └── criterion_grid.csv
    ├── figures/
    └── reports/
~~~

### Where to find shared project components

- **Base equations and parameter definitions:** [src/model/ma2024_model.py](src/model/ma2024_model.py)
- **Base parameter values:** 'DEFAULT_PARAMS' in [src/model/ma2024_model.py](src/model/ma2024_model.py)
- **Supplementary-information randomized ranges:** 'SI_RANDOMIZED_RANGES' in
  [src/model/ma2024_model.py](src/model/ma2024_model.py)
- **Dose ranges and initial conditions:** 'SI_DOSE_RANGES' and 'SI_INITIAL_CONDITIONS' in
  [src/model/ma2024_model.py](src/model/ma2024_model.py)
- **Log-space integration and endpoint definitions:** [src/analysis/sweep_engine.py](src/analysis/sweep_engine.py)
- **Structured parameter-grid execution:** [src/analysis/grid_runner.py](src/analysis/grid_runner.py)
- **Shared color and plot conventions:** [src/visualization/plot_style.py](src/visualization/plot_style.py)
- **Reproducible walkthrough:** [notebooks/01_sweep_execution.ipynb](notebooks/01_sweep_execution.ipynb)
- **Existing tabular outputs:** [results/tables/](results/tables/)
- **Future generated plots:** [results/figures/](results/figures/)
- **Future written reports:** [results/reports/](results/reports/)

The 'results/reports/' directory is intentionally empty in the repository. Put member reports,
combined reports, and report-specific supplementary documents there using descriptive names such
as 'person3_inhibitor_pharmacodynamics.pdf' or 'final_group_report.pdf'. Do not overwrite another
member’s report.

## Running the project

The project is intended to run in the 'biobase' environment:

~~~bash
conda activate biobase
~~~

From the repository root, import the shared package with:

~~~python
from src.model.ma2024_model import DEFAULT_PARAMS, simulate
from src.analysis.sweep_engine import compute_endpoints
~~~

To reproduce the model’s standalone validation plots:

~~~bash
conda run -n biobase python -m src.model.ma2024_model
~~~

Those plots are written to 'results/figures/'.

To use the notebook, open 'notebooks/01_sweep_execution.ipynb' from either the repository root
or the 'notebooks/' directory. The notebook locates the repository root, imports the package under
'src/', and writes its CSV outputs to 'results/tables/'.

When adding a new analysis:

1. reuse 'src.model.ma2024_model' rather than copying the ODEs;
2. reuse or extend 'src.analysis.sweep_engine' for endpoint definitions;
3. record parameter values, ranges, solver settings, and thresholds;
4. save tables under 'results/tables/';
5. save plots under 'results/figures/' in both PNG and PDF format;
6. save reports under 'results/reports/';
7. update this README only when a shared convention or repository location changes.

## References

The primary model reference is Ma et al. (2024), DOI:
[10.1038/s41467-024-52711-w](https://doi.org/10.1038/s41467-024-52711-w).

The code comments identify the corresponding main-text and supplementary-information equations.
The full proposal, Ma paper, and supplementary PDF remain in the parent project workspace used to
curate this repository.
