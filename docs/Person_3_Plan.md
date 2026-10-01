# Person 3 work plan: inhibitor pharmacodynamics

## 1. Scope and objective

This plan translates the **Mid term presentation** section of 'Group_Proposal.docx' into the work owned by Person 3.

The central question for this workstream is:

> How do inhibitor dose, inhibitor response steepness, and Bla-inactivation rate change antibiotic selection between susceptible and resistant populations?

The three primary parameters are:

- \(h_i\): inhibitor Hill coefficient;
- \(i\): dimensionless inhibitor dose/concentration;
- \(d_b\): inhibitor-dependent inactivation/turnover rate of free Bla.

They are grouped because they act through the same pharmacodynamic pathway. The first two determine the inhibitor activity \(\iota\), and all three affect the resistant-cell private benefit and the public antibiotic-degradation pathway. The work should report not only the final resistant fraction, but also whether the population reaches a resistant-dominant, susceptible-dominant, or mixed endpoint and how quickly it gets there.

This workstream is not responsible for redefining the whole Ma model or independently sweeping the other people’s parameters. Those parameters should initially be held at the shared baseline values, then varied only in explicitly labelled coordination/robustness runs.

## 2. Source material and model to implement

Use the following local files as the authoritative project inputs:

- [Group_Proposal.docx](./Group_Proposal.docx), especially the **Mid term presentation** section;
- [Ma_paper.pdf](./Ma_paper.pdf), especially the two-population model and selection-dynamics figures;
- [ma_supp.pdf](./ma_supp.pdf), especially Supplementary Information Sections 1–3 and Supplementary Figures S2–S5.

The baseline dimensionless ODE system is:

\[
\begin{aligned}
\frac{dn_s}{d\tau} &= (g-l)n_s,\\
\frac{dn_r}{d\tau} &= (\alpha g-\beta l)n_r,\\
\frac{ds}{d\tau} &= (\xi l-g)n_s+(\xi\beta l-\alpha g)n_r,\\
\frac{da}{d\tau} &= -\kappa_bba-\phi n_ra-d_aa,\\
\frac{db}{d\tau} &= \beta l n_r-d_b\iota b,
\end{aligned}
\]

with

\[
g=\frac{s}{1+s},\qquad
l=\gamma\frac{a^{h_a}}{1+a^{h_a}}g,
\]

\[
\iota=\frac{i^{h_i}}{1+i^{h_i}},\qquad
\beta=\beta_{\min}+c(1-\beta_{\min})\iota,
\]

\[
\phi=\phi_{\max}(1-c\iota).
\]

Interpretation:

- larger \(\iota\) means stronger inhibitor action;
- larger \(\beta\) means resistant cells lose more of their private lysis protection;
- smaller \(\phi\) means resistant cells provide less antibiotic degradation;
- larger \(d_b\) removes free Bla faster, shortening the duration of extracellular protection;
- \(h_i\) changes the steepness of the inhibitor response, not simply its strength.

Important: because \(i\) is scaled by the inhibitor half-maximal concentration, increasing \(h_i\) increases \(\iota\) when \(i>1\), but decreases \(\iota\) when \(0<i<1\). Every \(h_i\) result must therefore report the dose range and include the corresponding \(\iota(i,h_i)\) curve.

### Shared baseline values

Unless the group agrees on updated values, initialize the model with the supplementary-information values:

| Parameter | Baseline | Supplementary range/notes |
|---|---:|---|
| \(\alpha\) | 0.95 | 0.75–1 |
| \(\beta_{\min}\) | 0.90 | 0–1 |
| \(\xi\) | 0.80 | 0–1 |
| \(d_a\) | 0.02 | fixed in base model |
| \(\kappa_b\) | 0.35 | 0–1 |
| \(d_b\) | 1 | 1–10 in the paper’s randomized simulations |
| \(\gamma\) | 1.38 | 1.1–1.4 |
| \(h_a\) | 3 | 1–5 |
| \(i\) | dose sweep | paper uses 0.1–10 |
| \(h_i\) | 2 | 1–5 |
| \(\phi_{\max}\) | 1 | 0–5 |
| \(c\) | 0.70 | 0–1 |

Initial conditions for the base model are \(n_s(0)=n_r(0)=0.2\), \(s(0)=4\), \(b(0)=0\), and an antibiotic initial condition \(1<a(0)<100\). The exact \(a(0)\) used in each figure must be recorded because Person 2 owns the antibiotic-dose/lysis parameters.

## 3. Concrete responsibilities

### A. Implement and verify the Person 3 model pathway

1. Implement the five-state Ma model in a modular solver, with \((h_i,i,d_b)\) exposed as named inputs.
2. Implement helper functions for \(g,l,\iota,\beta,\phi\), so the inhibitor pathway can be unit-tested independently of the ODE solver.
3. Confirm that the code reproduces the expected qualitative behavior of the Ma model:
   - antibiotic exposure initially suppresses both populations;
   - resistant cells and Bla can reduce antibiotic exposure;
   - stronger inhibitor action reduces private protection and Bla-mediated protection;
   - treatment can end with either resistant or susceptible enrichment depending on parameters.
4. Reproduce at least one representative Ma-style time course and one two-dimensional treatment-response map before running new sweeps. This is a validation checkpoint, not a requirement to reproduce every paper figure.

### B. Characterize the inhibitor response itself

Before interpreting population dynamics, calculate and plot \(\iota(i,h_i)\), \(\beta(i,h_i)\), and \(\phi(i,h_i)\). This separates direct pharmacodynamic effects from ecological feedback.

Required checks:

- low-dose regime: \(i<1\);
- transition regime: \(i\approx1\);
- high-dose regime: \(i>1\);
- low, baseline, and high \(h_i\);
- low, baseline, and high \(d_b\).

This prevents a misleading conclusion that a larger \(h_i\) always means stronger inhibition.

### C. Run the grouped parameter study

Use one common simulation grid and one common plotting function for the three pairwise parameter planes:

1. **Dose–steepness plane:** \((i,h_i)\) at baseline \(d_b\).
2. **Dose–Bla-inactivation plane:** \((i,d_b)\) at baseline \(h_i\).
3. **Steepness–Bla-inactivation plane:** \((h_i,d_b)\) at a small, baseline, and/or large dose.

The third plane should be included only where it adds information; if it is nearly redundant with the first two, report that redundancy instead of expanding the number of figures.

Each plane should be evaluated with the same endpoint panel:

- final resistant fraction \(f_R\);
- final Shannon diversity of the two live populations;
- final total live biomass;
- time to numerical steady state;
- susceptible recovery status and time to recovery;
- resistant extinction status;
- peak population loss.

For each plane, save the raw parameter grid and all endpoint arrays so that color scales, thresholds, and representative points can be changed without rerunning the ODEs.

### D. Group trajectories by dynamical outcome

Do not make a separate trajectory figure for every parameter combination. From each parameter plane, select representative points based on the simulated outcomes:

- **susceptible-dominant endpoint:** low \(f_R\), with \(n_s\) remaining viable;
- **mixed/coexistence endpoint:** both populations remain above the extinction threshold and \(f_R\) is intermediate;
- **resistant-dominant endpoint:** high \(f_R\), whether or not total biomass recovers;
- **near-extinction or failed-recovery endpoint:** both populations remain below the viability criterion or total biomass fails to recover.

For each selected point, show \(n_s,n_r,s,a,b\), \(f_R(\tau)\), and \(\iota,\beta,\phi\) where useful. If one or more outcome classes do not occur in a parameter plane, state that explicitly rather than forcing a representative example.

### E. Connect inhibitor behavior to the initial susceptible/resistant ratio

The group’s central question asks how initial susceptible/resistant composition affects final diversity and time to the endpoint. Person 3 should contribute a focused composition check using a small set of inhibitor regimes selected from the response atlas:

- one regime predicted to favor susceptible cells;
- one mixed or near-boundary regime;
- one regime predicted to favor resistant cells, if present.

At each regime, vary the initial resistant fraction while keeping the initial total biomass fixed. For example, use \(f_{R,0}=0.05,0.25,0.50,0.75,0.95\), unless the group chooses a different shared grid. Report whether changing the initial ratio changes only the final magnitude of \(f_R\), or changes the direction of selection and endpoint class. This directly addresses Supplementary Figure S4A and the project’s stated question.

## 4. Endpoint definitions to agree on before sweeping

Write these definitions into the code and methods notes before producing final figures.

### Resistant fraction and diversity

\[
f_R(\tau)=\frac{n_r(\tau)}{n_s(\tau)+n_r(\tau)}.
\]

For two live populations, calculate Shannon diversity as

\[
H(\tau)=-f_S\ln f_S-f_R\ln f_R,
\qquad f_S=1-f_R,
\]

with the zero-abundance terms defined by continuity. Report both \(H\) and \(f_R\): maximum diversity does not necessarily mean high total biomass or successful treatment.

### Steady state

Define the numerical steady state as the first time after which the state changes remain below a documented tolerance for a documented time window. Verify that the reported endpoint is not simply the end of a simulation that stopped too early. If no convergence occurs, label the run as nonconverged rather than silently using its final time point.

### Extinction and viability

Use a scale-aware threshold, for example a fixed fraction of the initial total biomass or a threshold agreed upon by the group. Report the threshold in every figure caption. Distinguish:

- resistant extinction: \(n_r\) below threshold;
- susceptible extinction: \(n_s\) below threshold;
- whole-community collapse: \(n_s+n_r\) below threshold.

### Recovery

Define the no-treatment trajectory as the reference for viable recovery. Susceptible recovery time is the first time after the treatment-induced minimum at which \(n_s\) returns to a chosen fraction of the corresponding no-treatment susceptible abundance and remains there for the required window. Do not call recovery “carrying capacity” without specifying the reference, because the nutrient-based model does not use a separately imposed logistic carrying-capacity parameter.

### Peak population loss

Use a clearly stated normalization, preferably

\[
1-\frac{\min_\tau(n_s+n_r)}{n_s(0)+n_r(0)},
\]

and also retain the minimum absolute biomass. This distinguishes a deep bottleneck from a merely high final resistant fraction.

## 5. Figure and analysis package

The minimum Person 3 figure package should be:

1. **Model/pathway validation:** one inhibitor activity plot and one representative Ma-style trajectory.
2. **Pairwise inhibitor response atlas:** up to three identically formatted phase diagrams for \((i,h_i)\), \((i,d_b)\), and \((h_i,d_b)\), using final \(f_R\) as the primary map and companion maps for diversity, recovery, or biomass where they reveal a different conclusion.
3. **Outcome-class trajectories:** one grouped panel each for susceptible-dominant, mixed, resistant-dominant, and failed-recovery behavior when those classes are present.
4. **Initial-composition analysis:** \(f_{R,0}\) versus final \(f_R\), final diversity, and time to steady state for selected inhibitor regimes.
5. **Mechanistic bridge:** a compact plot linking \(i,h_i,d_b\) to \(\iota,\beta,\phi\), \(b(\tau)\), and \(a(\tau)\), showing why the endpoint changes.

If the phase diagrams are strongly redundant, combine them into a single response-atlas figure with shared axes and move redundant endpoint maps to supplementary material. The goal is to expose distinct behaviors, not to maximize the number of heatmaps.

## 6. Analysis questions and expected interpretations

Answer these questions in the results notes:

1. Does increasing inhibitor activity consistently lower resistant private protection through \(\beta\)?
2. Does it also reduce resistant-mediated antibiotic degradation through \(\phi\), and does that create a delayed loss of both populations?
3. Can two inhibitor settings produce similar final total biomass but different resistant fractions, as in the Ma paper’s dose-response results?
4. Does \(d_b\) mainly alter the duration of free-Bla protection, the final composition, or both?
5. Does \(h_i\) matter independently of \(\iota\), or is its apparent effect mostly a consequence of where the dose lies relative to \(i=1\)?
6. Are there parameter regions in which the resistant population is outnumbered but not extinguished? Report this separately from resistant clearance.
7. Are there regions in which total biomass recovers while diversity or composition becomes unfavorable?
8. Does the initial susceptible/resistant ratio change the direction of selection, or primarily shift the endpoint fraction and time scale?

Use the Ma selection criterion as an interpretation aid, not as a replacement for the dynamic simulation. The general criterion is

\[
1-\beta > \frac{1-\alpha}{l/g},
\]

and the saturating-dose simplification is

\[
(1-c)(1-\beta_{\min})>\frac{1-\alpha}{\gamma}.
\]

Because Person 3 changes the time-dependent inhibitor activity and Bla turnover, the full ODE trajectories are needed to capture the changing antibiotic and Bla environment.

## 7A. Expanded local analysis requested for this project iteration

The following additions define the exhaustive Person 3 analysis package. They were first developed locally and are now included in this repository for review; results remain conditional on the shared assumptions and dependency choices.

### A. Selection-relative endpoint metrics

In addition to the endpoint resistant fraction, report metrics relative to the initial composition:

- Δf_R = f_R(final) − f_R(0);
- change in log resistant-to-susceptible ratio, Δlog(n_r/n_s);
- selection direction classified as resistant-enriched, susceptible-enriched, or neutral using a documented tolerance;
- final Shannon diversity and its change from the initial value;
- final total biomass, minimum total biomass, peak population loss, free-Bla area under the curve, and recovery time.

This prevents a high final f_R from being interpreted as resistant selection when it is only the consequence of starting with a high resistant fraction or when both populations have collapsed.

### B. Shared antibiotic anchors and effective inhibitor variables

Use the Person 2 antibiotic-dose anchors as a comparison set, with the Ma-style baseline a₀ = 2.15, an intermediate anchor, and a high anchor near 100 when the group confirms them. Person 3 owns the inhibitor response at each anchor; Person 3 does not unilaterally redefine the antibiotic-dose sweep.

For every run retain the effective variables:

- ι = iʰⁱ/(1 + iʰⁱ);
- β, the resistant-cell private-lysis factor;
- φ, the resistant-cell public antibiotic-degradation factor;
- d_bι, the effective free-Bla inactivation rate.

This makes clear that i and h_i are not independently identifiable from a single static ι value, while d_b changes the time scale of free-Bla loss even when ι is held fixed.

### C. Mechanistic ablation panel

For selected low-selection, near-neutral, and high-selection regimes, compare:

1. the full inhibitor pathway;
2. private-protection-only variation through β;
3. public-degradation-only variation through φ;
4. free-Bla-turnover-only variation through d_bι;
5. an inhibitor-independent control.

Report the change in resistant fraction, log-ratio selection, free-Bla exposure, antibiotic exposure, bottleneck depth, and recovery. These ablations are the bridge between the pharmacodynamic functions and the ecological endpoint.

### D. Expanded response atlas

For each confirmed antibiotic anchor, generate identically formatted maps for:

- (i, h_i), at baseline d_b;
- (i, d_b), at baseline h_i;
- (h_i, d_b), at low, baseline, and high i.

The primary map is Δf_R or Δlog(n_r/n_s), with companion maps for final f_R, diversity, free-Bla exposure, recovery, total biomass, and peak loss. Retain every grid row in machine-readable tables. If the same structure appears across anchors, summarize that redundancy quantitatively rather than treating repeated visual patterns as independent evidence.

### E. Initial-composition and robustness analyses

At fixed initial total biomass, vary f_R(0) over a broad grid in selected susceptible-enriched, near-neutral, and resistant-enriched regimes. Determine whether composition changes the direction of selection, only shifts the endpoint fraction, or changes time scales and recovery.

Run a reduced Latin-hypercube/random robustness sample over the non-Person 3 parameter ranges from the supplement. Record the high-dose minus low-dose selection effect, rank correlations, and the fraction of samples retaining the baseline interpretation. Keep this separate from the primary Person 3 conclusion because it depends on the other people’s parameters.

### F. Identifiability and confounding check

Construct matched-ι settings with different (i, h_i) pairs. Compare trajectories and endpoints within matched-ι groups. Then vary d_b at fixed ι to quantify which observations distinguish inhibitor activity from free-Bla turnover. The expected interpretation is that endpoint composition may be similar for matched ι, while the b(t), a(t), and recovery trajectories retain information about d_b.

### G. Explicit dependencies on the other project roles

- Person 1 supplies α, β_min, c, the shared selection criterion, and the robustness ranges or alternative growth assumptions.
- Person 2 supplies the final a₀/γ/h_a anchors, antibiotic time-course convention, and the cross-comparison dose grid.
- Person 4 supplies κ_b, d_a, φ_max, ξ uncertainty, direct a/b observations, and the sampling schedule or experimental observables.
- The whole group must agree on state ordering, solver tolerances, simulation horizon, extinction/recovery thresholds, and the shared output schema.

Runs requiring these inputs are labelled as coordination/robustness analyses rather than presented as independent Person 3-owned conclusions. The report preserves these caveats pending group coordination.

## 7. Implementation and data products

Create or contribute the following reusable components in the shared repository:

- 'model.py' or equivalent: Ma base-model RHS and named pharmacodynamic functions;
- 'simulate.py': one-run solver with initial-condition and parameter overrides;
- 'sweeps.py': pairwise grids and initial-composition sweeps;
- 'metrics.py': \(f_R\), Shannon diversity, steady-state, extinction, recovery, and bottleneck metrics;
- 'figures_person3.py': standardized phase diagrams and grouped trajectory panels;
- a parameter/configuration file containing every baseline value, grid, solver tolerance, stopping rule, and threshold;
- machine-readable output files containing parameter values, convergence status, and all endpoint metrics;
- a short methods note explaining which equations are reproduced from Ma/its supplement and which plotting or endpoint definitions are project additions.

Every simulation record should retain:

'h_i', 'i', 'd_b', all frozen parameters, initial conditions, solver settings, convergence status, final state, final \(f_R\), final \(H\), time to steady state, extinction flags, recovery time, and peak loss.

## 8. Suggested execution order

### Phase 1 — Shared model contract

- Confirm the exact equations, units/scaling, \(a(0)\), parameter defaults, solver, and endpoint thresholds with the group.
- Make sure all people use the same state ordering and parameter names.
- Decide whether the shared branch uses the base model only or also the paper’s revised single-population fitting model. Person 3’s main analysis should use the base two-population model.

### Phase 2 — Validation

- Unit-test \(\iota,\beta,\phi\) against hand-calculated values.
- Test no-inhibitor/low-inhibitor and high-inhibitor limiting behavior.
- Reproduce a representative Ma trajectory and treatment-response map.
- Check solver convergence and nonnegative state variables.

### Phase 3 — Inhibitor response atlas

- Generate the common pairwise grids.
- Calculate the full endpoint panel for every grid point.
- Identify outcome classes from the simulations.
- Select representative trajectories algorithmically from those classes.

### Phase 4 — Initial-composition check and robustness

- Run the selected inhibitor regimes over the initial resistant fraction grid.
- Repeat a reduced version of the atlas at one alternative antibiotic dose or lysis setting supplied by Person 2.
- Test whether conclusions survive modest changes in simulation horizon, steady-state tolerance, and extinction threshold.

### Phase 5 — Integration

- Provide the group with the final figures, raw data, parameter table, and concise interpretation.
- Explain which conclusions are specific to inhibitor pharmacodynamics and which depend on Person 1, 2, or 4 parameters.
- Add a limitations paragraph covering dimensionless scaling, no spatial structure, no horizontal gene transfer in the base model, and the ambiguity between endpoint dominance and true extinction.

## 9. Completion checklist

- [x] Base Ma ODE model runs with the shared initial conditions and defaults.
- [x] \(\iota,\beta,\phi\) unit tests pass, including the \(i<1\) versus \(i>1\) behavior of \(h_i\).
- [x] One Ma-style trajectory and one response map are qualitatively reproduced.
- [x] Pairwise \((i,h_i)\), \((i,d_b)\), and justified \((h_i,d_b)\) sweeps are complete or documented as redundant.
- [x] Every sweep reports \(f_R\), diversity, biomass, steady-state status, extinction, recovery, and peak loss.
- [x] Representative trajectories are grouped by observed outcome rather than by arbitrary parameter values.
- [x] Initial resistant-fraction analysis is complete for selected inhibitor regimes.
- [x] Thresholds, convergence rules, grids, and frozen parameters are documented.
- [x] Raw simulation outputs and figure-generation code are saved.
- [ ] Results are handed to the group with a short interpretation of the inhibitor mechanism and its limitations.

## 10. Local expanded-analysis status

The expanded analysis was run with `biobase`; its code and results are now included in this repository. The environment file documents the minimal `ar_cp` dependencies for other group members.

### Completed locally

- [x] Added selection-relative metrics: Δf_R, Δlog(n_r/n_s), selection direction, diversity change, biomass AUC, free-Bla AUC, recovery, peak loss, and pathway variables.
- [x] Used a log-population solver formulation for n_s and n_r so high-dose bottlenecks do not erase the resistant/susceptible composition signal numerically.
- [x] Ran the (i, h_i), (i, d_b), and h_i–d_b planes at low, baseline, and high inhibitor doses for antibiotic anchors a₀ = 2.15, 10, and 100.
- [x] Ran algorithmic representative selection for susceptible-enriched, near-neutral, resistant-enriched, and deepest-bottleneck regimes.
- [x] Ran fixed-total-biomass initial-composition analysis over f_R(0) = 0.05–0.95 for selected regimes.
- [x] Ran full/private-only/public-only/free-Bla-turnover-only/inhibitor-independent mechanistic ablations.
- [x] Ran a 48-sample Latin-hypercube robustness analysis over the supplementary non-Person 3 parameter ranges.
- [x] Ran matched-ι identifiability checks and a d_b sweep at fixed ι.
- [x] Generated the expanded figure package, each plot in both PNG and PDF, plus machine-readable tables and `results/reports/Person_3_Expanded_Report.pdf`.
- [x] Verified 3,159 atlas simulations completed with no solver failures.

### Requires input from other persons before shared conclusions

- **Person 1:** confirm α, β_min, c, the selection criterion convention, and whether the randomized ranges or growth assumptions should be changed.
- **Person 2:** confirm the final a₀, γ, and h_a anchors, the treatment/time-course convention, and whether the high-dose horizon is appropriate.
- **Person 4:** confirm κ_b, d_a, φ_max, ξ uncertainty, direct a/b observables, and the measurement schedule to which recovery/Bla metrics should be compared.
- **Whole group:** approve the state ordering, solver tolerances, simulation horizons, extinction and recovery thresholds, and the shared output schema.

These coordination items are deliberately listed separately because changing them can alter the interpretation of Person 3’s results. The expanded script is `src/analysis/person3_expanded_analysis.py`; its outputs are under `results/figures/`, `results/tables/`, and `results/reports/`.

## 11. Matched breadth to the other analysis notebooks

To make Person 3’s package comparable in breadth—not parameter values—to the completed analyses from the other members, the following supplementary analysis was added. The three inhibitor parameters remain the focus; the other model parameters are held at the shared baseline except in the previously documented robustness and coordination checks.

### Completed matched analyses

- [x] Ran a full 21 × 21 × 21 structured joint grid at a₀ = 2.15: 9,261 parameter combinations over log-spaced i ∈ [0.1, 10] and linearly spaced hᵢ ∈ [1, 5], dᵦ ∈ [1, 10]. The existing three-anchor response atlas and prior Person 3 analyses are retained as well.
- [x] Saved per-run endpoints and diagnostics: Δf_R and Δlog(nᵣ/nₛ), diversity, final and minimum populations, peak loss, biomass and pathway AUCs, antibiotic clearance, susceptible viability/recovery (both treatment-relative and no-treatment-relative), resistant extinction, and shared steady-state status.
- [x] Added all three pairwise inhibitor planes for eight endpoints: final diversity, time to steady state, final f_R, susceptible recovery time, peak population loss, final total population, final nₛ, and final nᵣ. These are projections of the complete 3D grid, with the third parameter averaged over and explicitly identified in the data/report.
- [x] Added the team-notebook-style analysis families: marginal binned mean ± SD; Spearman rank associations; standardized main-effect and pairwise-interaction regressions; endpoint tradeoffs; absolute population/composition summaries; regime counts; convergence diagnostics; and selection-criterion-versus-dynamic-outcome validation.
- [x] Added grouped representative trajectories selected algorithmically from full-grid outcomes, plus low/middle/high one-factor inhibitor trajectories and clearance summaries. Full endpoint simulations use τ = 300; representative transient plots show the first 60 τ for legibility.
- [x] Checked the log-population implementation against direct linear-state integration on 12 randomized parameter combinations; all 12 agree within the documented numerical tolerance.
- [x] Rebuilt the expanded report with the original Person 3 analyses and the matched comprehensive appendix. Every generated plot has PNG and PDF versions under `results/figures/`.

### What the completed inhibitor grid says (conditional on this baseline)

Across the 9,261 runs, final resistant fraction ranges from about 0.407 to 0.509, while all runs remain mixed under the predeclared dominance threshold (f_R < 0.25 or > 0.75); this grid contains neither resistant dominance nor extinction. The direction and size of selection still vary: Δf_R ranges from approximately −0.093 to +0.009. Inhibitor dose i is the strongest monotonic correlate across the grid (Spearman ρ ≈ −0.91 for Δf_R); hᵢ and dᵦ have weaker marginal rank associations over these chosen ranges. Static selection-criterion direction agrees with the dynamic Δf_R direction in approximately 78.7% of runs, so it is informative but not a substitute for the trajectories. These are conditional model results, not experimental predictions; the agreed parameter ranges, endpoint thresholds, and shared baselines should be confirmed before group-level claims.

The detailed results and interpretation are in [Person_3_Expanded_Report.pdf](../results/reports/Person_3_Expanded_Report.pdf). The reproducible matched-analysis entry point is `src/analysis/person3_comprehensive_analysis.py`; full-grid and summary tables are in `results/tables/` with the `comprehensive_` prefix. Both scripts and all generated artifacts are included in the repository.
