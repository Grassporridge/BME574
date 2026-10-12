# gLV antibiotic analysis

This is a standalone extension of the antibiotic analyses in
`consumer_resource_model_antibiotic/`. It uses a randomly sampled, diagonally
stable generalized Lotka–Volterra (gLV) community instead of an explicit
consumer–resource system. The antibiotic acts through direct species-specific
Hill killing and first-order decay. There are no β-lactamase, enzyme-mediated
drug removal, or inhibitor terms.

## Equations and interpretation

For species `i = 1,…,n`:

```text
dXᵢ/dt = Xᵢ [ rᵢ(1 − cρᵢ) + Σⱼ AᵢⱼXⱼ − kmax A(t)^h/(IC50ᵢ^h + A(t)^h) ]
dA/dt  = −kdeg A,              A(0) = A₀
IC50ᵢ = IC50min (IC50max/IC50min)^ρᵢ
```

`Xᵢ` is species abundance; `rᵢ` intrinsic growth; `Aᵢⱼ` the effect of species
`j` on species `i` (`Aᵢᵢ < 0` is self-limitation); `ρᵢ ∈ [0,1]` is a sampled
resistance score; `c` is the maximum fractional growth cost at `ρ=1`; `A(t)` is
antibiotic concentration; `IC50ᵢ` is the species-specific half-maximal killing
concentration; `kmax` is the maximum killing rate; `h` is the Hill coefficient;
`kdeg` is drug decay; and `A₀` is initial drug concentration.

The interaction matrix is sampled constructively as
`A = W⁻¹(K − P)`, with positive diagonal `W`, skew-symmetric `K`, and symmetric
positive-definite `P`. This gives `WA + AᵀW = −2P ≺ 0` and therefore a direct
diagonal-stability certificate. Rows are scaled so `Aᵢᵢ = −1`. A positive
target equilibrium `x*` is sampled, and `r = −Ax*`; consequently the
zero-antibiotic, zero-resistance-cost gLV system has the feasible equilibrium
`x*`. The run saves `A`, `W`, `x*`, species traits, and numerical certificate
checks in `outputs/`.

The stability certificate applies to the autonomous gLV interaction matrix. It
does not, by itself, establish global stability of the full system under
time-varying antibiotic exposure, nor does it imply that the randomly generated
parameters are empirically calibrated.

## Analyses and outputs

Run `python run_analysis.py` from this folder. Results are written to
`outputs/`, plots to `outputs/figures/` (each plot is saved as both PNG and PDF),
and the report to `outputs/reports/gLV_Antibiotic_Analysis_Report.pdf`.

The suite includes:

1. No-antibiotic reference and sampled species traits.
2. Matched antibiotic dose response (`A₀ = 0, 0.1, 0.25, 0.5, 1, 2`), endpoint
   summaries, species outcomes, and a species-by-dose heatmap against matched
   untreated abundance.
3. Antibiotic-decay sweep (`kdeg = 0, 0.005, 0.01, 0.02, 0.05, 0.1`).
4. Resistance-selection comparison against the matched no-drug control.
5. Resistance-cost sweep (`c = 0, 0.1, 0.2, 0.3, 0.4`) with treated and untreated
   cases.
6. Representative-species growth-versus-killing mechanism plots.
7. Initial-biomass-by-dose factorial. This replaces the CRM resource-by-dose
   analysis: gLV has no explicit resource variable, while initial biomass
   changes the strength of density-dependent interaction feedback.
8. Interaction-strength-by-dose factorial using separate stable matrices.
9. A 40-community Spearman rank-correlation sensitivity screen over dose,
   decay, maximum killing, Hill coefficient, resistance cost, initial biomass,
   and interaction strength.
10. A ten-community seed robustness check over a dose series.
11. A controlled resistance-assignment experiment: keep the same interaction
    matrix, growth vector, equilibrium, and exact multiset of 25 resistance
    scores, then reassign those scores to species. Cases place high resistance
    on high/low intrinsic-growth species and high/low baseline-equilibrium
    abundance species, alongside 200 random permutations. Each assignment is
    compared to its own untreated control, accounting for the resistance cost.

Numerical ranges are exploratory normalized scenario choices, not fitted
biological values. The sensitivity screen is not a formal Sobol analysis, and
the resistant community members are sampled initially rather than evolving
resistance during the simulation.

## Files

- `model.py`: parameterization, stable community sampling, antibiotic dynamics,
  simulation, and summary metrics.
- `run_analysis.py`: complete analysis/figure/report workflow.
- `outputs/scenario_summary.csv`: pooled primary sweep summaries.
- `outputs/sensitivity_runs.csv` and `sensitivity_spearman_rho.csv`: sensitivity
  ensemble and rank correlations.
- `outputs/ensemble_robustness.csv`: across-seed robustness results.
- `outputs/resistance_assignment_permutations.csv`: structured and 200 random
  fixed-community resistance allocations and matched treatment/control metrics.
- `outputs/resistance_assignment_species_outcomes.csv`: species-level resistance
  assignments and outcomes for the named allocations and one example random draw.
- `outputs/reports/gLV_Final_Observations.md`: data-backed summary across all
  analyses, including interpretation and limitations.
- `outputs/figures/`: all analysis plots, in PNG and PDF.
- `outputs/reports/`: generated PDF report.

The repository's root `environment.yml` defines the `ar_cp` environment used by
the analyses. You can also run locally with an existing environment, e.g.
`conda run -n biobase python run_analysis.py`.
