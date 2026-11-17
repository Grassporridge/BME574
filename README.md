# BME574 antibiotic–resistance community project

This repository holds the models, analyses and reports for the BME574 group project. It asks how
β-lactam treatment reshapes a community of susceptible and resistant bacteria, and when a
β-lactamase (Bla) producer protects the rest of the community.

The project has had two phases:

| Phase | Model | Status | Where |
|---|---|---|---|
| 1 | Two-population model of Ma et al. (2024): susceptible and resistant cells, nutrient, drug, free Bla and an inhibitor | **complete** (Person 1–4 sweeps) | [Two-species-model-summary/](Two-species-model-summary/), plus `src/`, `notebooks/`, `results/` |
| 2 | 57-strain consumer-resource model of the barcoded *E. coli* cross-protection experiment (one *bla*TEM-1 carrier, k30, plus 56 susceptible strains) | **current** | [updated_consumer_resource_model/](updated_consumer_resource_model/) |

A generalized Lotka–Volterra model (`glv_antibiotic/`) was tried and dropped. Several earlier
25-species consumer-resource prototypes are also kept for reference; see the
[repository map](#repository-map).

**Start here:** [updated_consumer_resource_model/FINDINGS_SUMMARY.md](updated_consumer_resource_model/FINDINGS_SUMMARY.md)
(also as PDF) summarizes every phase-2 analysis so far, what each covers and what is still open.

## Project question

The original question (phase 1) was:

> How do antibiotic exposure, inhibitor action, resistance burden/private benefit, Bla-mediated
> antibiotic degradation, nutrient recycling, and initial susceptible/resistant composition combine
> to determine final community composition, diversity, recovery, and selection for resistance?

Phase 2 asks the same question of a specific experiment: the lab's barcoded community, grown for
4 h, treated with carbenicillin (50 or 70 µg/ml) until 24 h, diluted 1:100 into drug-free medium
and sampled at 4, 6, 24 and 48 h, in M9CA and in LB. In particular:

- How does a single resistant strain protect the 56 susceptible ones, and at what cost in takeover?
- How do the number, identity and growth of resistant strains change recovery?
- Why does resistance protect diversity in M9CA but not in LB?

The phase-2 model has **no β-lactamase inhibitor**, so the inhibitor part of the original question
is covered only by phase 1.

## Project status

**Phase 1.** The Person 1–Person 4 sweeps of the two-population model are complete; the
notebooks, tables and reports are in [Two-species-model-summary/](Two-species-model-summary/)
(notebooks 02–05 cover the four parameter blocks, and notebook 06 checks the paper's selection
criterion).

**Phase 2.** The model is
[updated_consumer_resource_model/resilience_model.py](updated_consumer_resource_model/resilience_model.py);
its equations, assumptions, parameters and readouts are documented in
[updated_consumer_resource_model/README (3).md](updated_consumer_resource_model/README%20(3).md).
(The `ReadME.md` file beside it is empty.) Completed analyses:

| Analysis | Folder (in `updated_consumer_resource_model/`) | Author |
|---|---|---|
| Number and identity of resistant strains, including a stochastic version | `test_resistance_distribution/` | Kylie |
| Private vs. public benefit of resistance | `test_protection_mechanism/` | Sarvesh |
| Medium ablation, M9CA vs. LB | `test_medium_ablation/` | Sarvesh |

The full write-up of the last two is `REPORT_protection_and_medium.md` (also as PDF). The
checklist below tracks what has been run.

**Not yet possible:** no experimental data (barcode counts, OD curves, measured strain traits) are
in the repository, only the lab-meeting slides. Comparing the model with the data is the main open
task.

## Running the 57-strain model

The phase-2 model needs only NumPy, SciPy, pandas and Matplotlib, so the 'ar_cp' environment
below works. From `updated_consumer_resource_model/`:

~~~python
from resilience_model import run, sweep, grid, check

check()                                  # once per session; 20/20 checks should pass
r = run()                                # the baseline experiment (M9CA, 50 µg/ml, k30 resistant)
r.metrics                                # one row of readouts at 48 h
r.timeseries                             # per-strain abundance at 4, 6, 24, 48 h
r = run(dose=70, growth_mean=1.4)        # anything not named keeps its baseline value
df = sweep("n_resistant", [1, 2, 4, 8])
df = grid("carrier_lysis", [0, 0.5, 1], "vmax", [1e3, 1e4])
~~~

Units are hours, OD and µg/ml of carbenicillin. Each analysis folder sits one level below the model
and imports it with `sys.path.insert(0, "..")`; open a notebook from its own folder and Run All.
Kylie's `test_resistance_distribution/` notebooks embed their own copy of the model and check its
hash, so they run on their own.

When adding a phase-2 analysis:

1. import `resilience_model` rather than copying its equations; add new terms to the model with a
   default that reproduces the existing behavior, plus a check in `check()`;
2. make a `test_<topic>/` folder with the notebook and a `<topic>_results/` folder for CSVs and
   `figures/`;
3. add a row to the "Characterizations" table in the model README, and tick or annotate the
   checklist below;
4. state the parameter values used, and flag any value that is assumed rather than measured.

## Possible tests: consumer-resource model

**Ready** means the test runs through 'sweep()' or 'grid()' with existing parameters; **small
change** or **larger change** means code is needed first. The example notebook already contains a
rough 'n_resistant' sweep (at 'gamma=4', 'vmax=3000'), single sweeps of 'dose', 'growth_mean' and
'gamma', and a 'gamma' × 'vmax' grid; the tests below extend those rather than repeat them.

### Community makeup

- [x] **Number of resistant strains:** 'n_resistant' 1 → 56 at several doses. How many carriers are
  needed to protect the rest, and is there a threshold? *Done at 50 µg/ml (Kylie): 0, 1, 5, 10, 20,
  30, 40, 57 carriers × 10 carrier orders, deterministic and stochastic; see
  'updated_consumer_resource_model/test_resistance_distribution/'. Other doses not yet run.*
- [ ] **Number vs. biomass of resistance:** the same total carrier biomass as one strain at 8×
  ('focal_multiplier') or eight strains at 1×. Does protection depend only on resistant biomass?
  *Ready.*
- [ ] **Starting resistant fraction:** 'focal_multiplier' 1 → 30 (k30 at 1.75% → 35%). Earlier
  detoxification versus k30 takeover. *Ready.*
- [ ] **Which strain carries resistance:** make the carrier the fastest, median or slowest grower.
  Does a fast carrier protect better but take over more? *Small change. Partly covered: random
  carrier identities in 'test_resistance_distribution/' (Kylie); the non-immune carrier's rescue
  tracks k30's growth rate across trait draws in 'test_protection_mechanism/' §7.*
- [ ] **Community size:** 'n_strains' 5, 10, 57, 200 at a fixed resistant fraction. *Ready (keep
  'focal_index' < 'n_strains').*
- [ ] **Trait heterogeneity:** 'growth_cv', 'ic50_cv'. Does more spread give graded strain loss
  instead of all-or-nothing outcomes? *Ready.*
- [ ] **Growth–sensitivity correlation:** draw IC50 correlated or anti-correlated with growth rate.
  *Small change.*
- [ ] **Robustness across trait draws:** 20+ values of 'seed'. *Ready. Partly covered: 12 seeds
  for the protection and medium-ablation scenarios; 'test_resistance_distribution/' uses seed 42
  only.*

### Resistance traits

- [x] **Cost of resistance:** 'resistance_cost' 0 → 0.5. *Done: × dose 10–500 µg/ml, rescue limits,
  four carrier types, 12 seeds; see 'updated_consumer_resource_model/test_cost_and_duration/'.*
- [ ] **Detoxification strength:** 'vmax', 'km'. How much drug-removal capacity does recovery need?
  *Ready. Partly covered: 'vmax' × 'carrier_lysis' in 'test_protection_mechanism/' §3, with a
  kinetic ceiling of ≈ 3e4 µg/ml/h/OD (the default 1e4 is at the ceiling). 'km' not swept.*
- [x] **Private vs. public benefit:** carrier immune only ('vmax=0'), degrades drug only (not
  immune), or both. *Done: new 'carrier_lysis' and released-enzyme terms (Ma et al.'s β and free
  Bla); see 'updated_consumer_resource_model/test_protection_mechanism/'.*
- [ ] **Partial resistance:** carriers with a finite IC50 instead of full immunity. *Small change.
  Partly covered: 'carrier_lysis' between 0 and 1 gives partial protection as a constant fraction
  (Ma et al.'s β), not as a shifted IC50.*

### Drug and treatment protocol

- [ ] **Dose response:** 'dose' 5 → 500, including 50 vs. 70 µg/ml. Locate the transitions between
  full recovery, partial strain loss and collapse. *Ready. Partly covered: 10–500 µg/ml against
  'resistance_cost' and 'treatment_h' in 'test_cost_and_duration/'.*
- [ ] **Growth phase at drug addition:** 'pregrowth_h' 0 → 8. Exponential vs. near-stationary
  treatment. *Ready.*
- [x] **Treatment length and transfer timing:** 'treatment_h'. Is the drug cleared before the 1:100
  transfer? *Done: 0.5–44 h × dose, weaker-protection scenarios, × cost; see
  'updated_consumer_resource_model/test_cost_and_duration/'. Note that the model's built-in
  'metrics' read at clock hours 24 and 48, so non-default schedules need schedule-relative readouts.*
- [ ] **Killing shape:** 'gamma', 'hill'. *Ready.*
- [ ] **Repeated treatment cycles:** several treat → transfer cycles. Does strain loss accumulate?
  *Small change.*
- [ ] **β-lactamase inhibitor (clavulanate).** *Larger change.*

### Environment (M9CA vs. LB)

- [x] **Growth rate × food:** 'grid("growth_mean", …, "food", …)', separating the two differences
  the LB preset changes together. *Done in 'test_medium_ablation/' (sections 1 and 4).*
- [ ] **Inoculum density:** 'initial_biomass', which sets the density at drug addition. *Ready.*
- [x] **Lag phase:** a new 'lag_h' parameter, the realistic version of the inoculum test.
  *Done: Baranyi–Roberts lag after inoculation and transfer.*
- [x] **Medium ablation:** switch growth, food, 'gamma', 'vmax' and lag one at a time from M9CA
  toward LB. Which difference causes LB's strain loss? *Done: one at a time, pairwise, and a full
  2⁵ factorial with Shapley attribution; see 'updated_consumer_resource_model/test_medium_ablation/'.*
- [ ] **Transfer bottleneck:** 'dilution' 10 → 1000. How much loss happens at the transfer?
  *Ready.*
- [ ] **Starvation death:** 'death_rate'. *Ready.*

### Measurement and readouts

- [ ] **Detection limit:** 'detection_rel' 1e-5 → 1e-3. *Ready.*
- [ ] **Readout time:** 24 h vs. 48 h, plus longer 'recovery_h' or extra transfers. Is loss
  permanent? *Ready (extra transfers is a small change).*
- [ ] **Protection-score distribution:** compare with the lab-meeting results at 50 and 70 µg/ml.
  *Ready.*

### Analytical checks and structural changes

- [ ] **Small cases (1S, 1R + 1S, 1R + nS):** 'n_strains' = 1, 2, 5, 10, to test derived survival
  conditions against simulation in the style of notebook 06. *Ready.*
- [ ] **Time rescaling:** scale all rates by one factor to test whether growth rate alone can change
  the outcome. *Ready.*
- [ ] **Stochastic transfer and persisters:** binomial sampling at the split and a persister
  subpopulation. *Small change. Partly covered: demographic noise (tau-leaping) in
  'test_resistance_distribution/stochastic_community_resilience.ipynb' (Kylie), negligible at
  10⁹ cells/OD; transfer sampling and persisters not yet added.*
- [ ] **Cross-feeding:** compare with the earlier consumer-resource model's metabolite structure.
  *Larger change.*
- [ ] **Pathogen invasion (PA14) after recovery.** *Larger change.*

## Shared endpoints

### Phase 2 (57-strain model)

These come from `run().metrics`. All are at 48 h (after the 1:100 transfer and 24 h drug-free
regrowth), and "treated" is compared with a matched untreated run of the same community unless
stated. "Protected" strains are those without the resistance gene.

| Readout | Definition | Use it for |
|---|---|---|
| `protected_richness` | protected strains above the detection limit (10⁻⁴ relative abundance) | survival, comparable with the lab's richness |
| `protected_extinct` | protected strains truly gone (below one cell) | telling extinction from non-detection |
| `protected_biomass_ratio` | protected biomass, treated ÷ untreated | how much of the susceptible community regrew |
| `compositional_resilience` | 1 − Bray–Curtis distance between treated and untreated communities | recovery of the composition |
| `resistant_share` | carrier fraction of the treated community | takeover by the resistant strain |
| `protected_change_4_6h` | protected biomass at 6 h ÷ at 4 h | the lab's early-decline panel |
| `drug_clearance_h` | hours from drug addition until the drug falls below the median susceptible IC50 | mechanism: how fast the community detoxifies |

Conventions:

- **Always report takeover with recovery.** Survival, composition and carrier share can move in
  opposite directions; a community can keep every strain while the carrier takes over.
- **Report detected and extinct strains separately.** A strain below the detection limit is not
  necessarily extinct.
- **Say which reference you use.** The metrics above compare treated with untreated at 48 h.
  `test_resistance_distribution/` instead compares each community with itself at 4 h (Shannon
  retention, Bray–Curtis to baseline), adjusted by the untreated control. Both are valid, but
  they are not interchangeable.

### Phase 1 (two-population model)

All workstreams should use the same definitions unless a report explicitly justifies an alternative:

- final resistant fraction: fᵣ = nᵣ/(nₛ + nᵣ);
- two-population Shannon diversity: H = −fₛ ln(fₛ) − fᵣ ln(fᵣ);
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
├── README.md                          this file
├── environment.yml                    the 'ar_cp' conda environment
│
├── updated_consumer_resource_model/   PHASE 2: current work
│   ├── resilience_model.py            the 57-strain model
│   ├── README (3).md                  model documentation (ReadME.md is empty)
│   ├── FINDINGS_SUMMARY.md / .pdf     summary of all phase-2 analyses
│   ├── REPORT_protection_and_medium.md / .pdf
│   ├── example_usage_executed.ipynb
│   ├── test_resistance_distribution/  number and identity of resistant strains (Kylie)
│   ├── test_protection_mechanism/     private vs. public benefit
│   └── test_medium_ablation/          M9CA vs. LB
│
├── Two-species-model-summary/         PHASE 1: final notebooks, sweep package, results, reports
├── src/                               PHASE 1: shared package (model, sweep engine, plotting)
├── notebooks/                         PHASE 1: sweep notebooks
├── results/                           PHASE 1: tables, figures, Person 3 and sweep reports
├── docs/                              PHASE 1: Person 3 plan
│
├── glv_antibiotic/                    dropped gLV model (kept for reference)
├── consumer_resource_model/           earlier 25-species prototype (Kylie, Oct)
├── consumer_resource_model_antibiotic/  25-species prototype with antibiotic and resistance cost
├── consumer_resource/                 25-species prototypes: baseline cross-feeding and batch models
├── consumer_resource_aligned/         25-species no-antibiotic control matched to the prototype above
│
├── param_sweep/                       copies of phase-1 sweep scripts and notebooks
├── trash/                             superseded notebooks
├── presentation_midterm/              placeholder
└── ma2024_model.py, *.zip             loose copies; the maintained versions are in src/ and
                                       Two-species-model-summary/
~~~

The 25-species prototypes, `param_sweep/`, `trash/`, the zip files and the root-level
`ma2024_model.py` (which differs from `src/model/ma2024_model.py`) are superseded. They are left
in place until the group agrees what to delete.

## Phase 1 reference: the two-population model

This section documents the Ma et al. model used in phase 1. It is not used by the phase-2 model.

The model has two live populations, 'n_s' (susceptible) and 'n_r' (resistant, Bla-producing), and
tracks a shared nutrient 's', antibiotic 'a' and free extracellular Bla 'b'. It is dimensionless,
so numerical values should be interpreted according to the scaling in Ma et al., not as direct
clinical concentrations without an additional dimensional mapping. The source PDFs (paper,
supplement, proposal) are in the parent project directory, outside this repository.

### Equations

The five-state ODE system is:

~~~text
dnₛ/dτ = (g − l)nₛ
dnᵣ/dτ = (αg − βl)nᵣ
ds/dτ  = (ξl − g)nₛ + (ξβl − αg)nᵣ
da/dτ  = −κ_bba − φnᵣa − dₐa
db/dτ  = βlnᵣ − dᵦιb
~~~

The closure relations are:

~~~text
g = s/(1 + s)
l = γ[aʰᵃ/(1 + aʰᵃ)]g
ι = iʰⁱ/(1 + iʰⁱ)
β = βₘᵢₙ + c(1 − βₘᵢₙ)ι
φ = φₘₐₓ(1 − cι)
~~~

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

~~~text
1 − β > (1 − α)/(l/g)
~~~

Under saturating antibiotic and inhibitor concentrations, the simplified criterion is:

~~~text
(1 − c)(1 − βₘᵢₙ) > (1 − α)/γ
~~~

These criteria are interpretation tools. Dynamic simulations are still required because the
populations continuously change the antibiotic, Bla, and nutrient environment.

### Shared baseline values and parameter ranges

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

~~~text
nₛ(0) = 0.2,  nᵣ(0) = 0.2,  s(0) = 4,  b(0) = 0,  a(0) = a₀
~~~

### Where to find the two-species components

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
- **Tabular outputs:** [results/tables/](results/tables/)
- **Generated plots (PNG and PDF):** [results/figures/](results/figures/)
- **Written reports:** [results/reports/](results/reports/)
- **Person 3 work plan and coordination items:** [docs/Person_3_Plan.md](docs/Person_3_Plan.md)

Member reports, combined reports, and report-specific supplementary documents belong in
'results/reports/' with descriptive names. Do not overwrite another member’s report.

### Running the two-species code

The project is intended to run in the minimal 'ar_cp' environment defined in
[environment.yml](environment.yml). It contains Python 3.11, NumPy, SciPy, pandas, Matplotlib,
reportlab for PDF report generation, and ipykernel for notebook-kernel support.

~~~bash
conda env create -f environment.yml
conda activate ar_cp
~~~

From the repository root, import the shared package with:

~~~python
from src.model.ma2024_model import DEFAULT_PARAMS, simulate
from src.analysis.sweep_engine import compute_endpoints
~~~

To reproduce the model’s standalone validation plots:

~~~bash
conda run -n ar_cp python -m src.model.ma2024_model
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
7. update this README when a shared convention or repository location changes.

## References

- Ma HR, Xu HZ, Kim K, Anderson DJ, You L. “Private benefit of beta-lactamase dictates selection
  dynamics of combination antibiotic treatment.” *Nature Communications* 15, 8337 (2024).
  DOI: [10.1038/s41467-024-52711-w](https://doi.org/10.1038/s41467-024-52711-w). The code comments
  identify the corresponding main-text and supplementary-information equations.
- The phase-2 model's references (TEM-1 kinetics, killing ∝ growth, lag phase) are listed in
  [updated_consumer_resource_model/README (3).md](updated_consumer_resource_model/README%20(3).md).
- The lab-meeting slides describing the experiment (18 Sep 2026) are in the parent project
  directory, outside this repository.
