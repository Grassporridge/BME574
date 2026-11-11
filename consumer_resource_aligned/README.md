# Antibiotic-aligned consumer-resource control

This is a new copy based on the uploaded BME574 repository. It matches the antibiotic
model with `initial_antibiotic=0`, with no antibiotic variable or killing equation.
It retains balanced feeding categories, resistance/IC50 traits, resistance growth
cost (default 0.20), Monod uptake, and R0-to-M1/M2 cross-feeding. Basal and starvation
death are zero. Legacy mortality settings do not change the ODE. IC50 has no direct
effect here. Resistance is fixed, not evolving; no beta-lactamase protection is modeled.

## Run

Extract this folder, open a terminal in it, then:

```bash
python -m pip install -r requirements.txt
python check_model.py
python run_simulation.py --seed 42 --replicates 10 --composition --output outputs
```

For a shorter run, use `--replicates 0` and omit `--composition`.
Alternatively open `revised_analysis.ipynb` and Run All. It embeds the updated
model and analysis code, and writes to `outputs_notebook`.

## What is included

- model.py: community generation, equations, integration, species summaries and CSV exports.
- analysis_tools.py: five scenarios, diversity, trait plots, composition/resource dynamics,
  growth/death fluxes, accounting/extinction checks, and across-seed comparisons.
- run_simulation.py: CLI workflow; optional starting-composition experiment.
- run_initial_composition.py: same species with 13 starting mixtures and five figures.
- check_model.py: matched traits, no-food behavior, energy budgets and scenario checks.
- check_alignment.py: numerical comparison to your antibiotic model (provide its path).
- revised_analysis.ipynb: self-contained notebook, cutoff sensitivity and composition experiment.

## Scenario meanings

The former starvation experiment is replaced with an explicit resistance-cost experiment:
baseline = fixed K, no cost; affinity = variable K, no cost; cost = fixed K, cost 0.20;
both = variable K, cost 0.20; structured = both plus heterogeneous initial abundances
and secretion/leakage types. All use balanced food categories. `both` matches the
antibiotic model's default ecological settings with A=0. Structured also matches its
corresponding flag, including its RNG sequence. Changing cost intentionally changes
uptake coefficients while intrinsic community traits remain fixed.

## Important interpretation

Biomass cannot decline from mortality in this control. Species ordinarily persist;
resource exhaustion stops growth instead of killing cells. The old starvation-driven
extinction results are intentionally not reproduced. Richness can stay fixed while
Shannon diversity and composition change. Time/biomass are normalized, not calibrated.
Total living biomass + all three resources + loss counters stays at the initial budget.
Previously saved outputs are not included; regenerate figures and tables with this code.

## Optional direct comparison

```bash
python check_alignment.py --antibiotic-model /path/to/consumer_resource_model_antibiotic/model.py
```

This checks default, fixed-affinity/no-cost, and structured controls with the same
seed and duration. The comparison source is your uploaded antibiotic model snapshot;
later changes to that source may require updating the alignment.
