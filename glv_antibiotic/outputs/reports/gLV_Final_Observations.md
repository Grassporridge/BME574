# Final observations: antibiotic-treated gLV system

## Scope and model interpretation

This report summarizes a 25-species gLV community with a randomly sampled, diagonally stable interaction matrix and direct species-specific Hill killing. Antibiotic concentration decays by first-order kinetics. There are no explicit resource, beta-lactamase, inhibitor, or evolving-resistance states. Resistance score rho is assigned at initialization; its fitness cost is included in intrinsic growth.

The sampled matrix passed the diagonal-stability certificate: maximum eigenvalue of WA + A^T W = -2 (must be negative), minimum eigenvalue of P = 1. The target equilibrium residual was 0. This certifies the autonomous gLV interaction matrix, not the full antibiotic-forced dynamics.

## Main observations

- **Dose response (one fixed community and resistance assignment):** final biomass was 4.665 without antibiotic, 4.207 at A0=1, and 3.037 at A0=2. At dose 2, 24 of 25 species remained above the numerical extinction threshold versus 25 in control; Shannon diversity shifted from 3.214 to 2.828.
- **Resistance enrichment in the dose series:** biomass-weighted rho was 0.520 in control and 0.628 at dose 2. This indicates ecological enrichment of initially more resistant members in this sampled community, not evolution of resistance.
- **Drug persistence:** with A0=1, final biomass was 0.077 at kdeg=0 and 4.665 at kdeg=0.1. Faster drug decay shortens the killing period; this is a matched comparison with other settings fixed.
- **Resistance cost:** in the untreated control, raising the maximum cost from 0 to 0.4 changed final biomass 5.203 -> 4.127. Under A0=1, the corresponding values were 4.841 -> 3.502. Interpret treated values alongside their same-cost controls because cost acts even without antibiotic.
- **Initial biomass x dose:** at A0=1, final biomass ranged from 4.137 to 4.349 across the four starting biomass levels. In gLV this is a density-feedback context test, not a resource-availability test.
- **Interaction-strength x dose:** at A0=1, the four separately sampled stable communities produced final biomass from 2.426 to 3.248 over interaction scales 0.1–0.55. Since each scale uses a different sampled matrix and associated feasible growth vector, this is an ensemble comparison, not a one-parameter perturbation of a single A.
- **Sensitivity screen:** 40 joint parameter/community draws gave the largest absolute Spearman associations with final biomass for final_shannon (rho=0.86), survivors (rho=0.73), final_weighted_resistance (rho=-0.65). These are exploratory screening signals, not causal effects or formal variance-based sensitivity indices.
- **Across-community robustness:** across ten independently sampled communities, final-biomass median (range) by dose was A0=0: 4.602 (4.381–4.772); A0=0.5: 4.213 (3.913–4.700); A0=1: 3.470 (2.714–3.986); A0=2: 2.555 (1.310–3.042). These communities also resample resistance assignment, so this spread includes both ecological and trait-assignment differences.

## Which species receive higher initial resistance?

To isolate placement, this experiment held A, r, x*, the multiset of 25 rho values, dose, decay, and cost fixed while changing only which species received each rho. Each assignment was run both untreated and treated, so antibiotic response is measured relative to that assignment's own baseline and the resistance cost is not mistaken for drug benefit.
- Pairing the highest rho values with the highest intrinsic-growth species gave a treated/control final-biomass ratio of 0.854; pairing high rho with the lowest-growth species gave 0.999.
- Pairing high rho with the highest versus lowest baseline-equilibrium-abundance species gave treated/control ratios 0.824 and 0.967, respectively.
- Across 200 random permutations, the treated/control biomass ratio ranged from 0.793 to 0.983, with median 0.883 and 5th–95th percentiles 0.827–0.937. The lowest-ratio placement was random_198 (0.793); the highest was random_101 (0.983).
- The assignment-specific response can differ because species vary in intrinsic growth, equilibrium abundance, and interaction-network position. The structured high/low-growth and high/low-equilibrium-abundance placements ask whether resistance is allocated to different ecological roles; the random-permutation distribution estimates variation across arbitrary allocations.

## Limitations and interpretation

All parameter values are exploratory normalized choices, not fitted biological estimates. Numerical values of kmax and time scales are not directly comparable between CRM and gLV without scale calibration. Results depend on one primary sampled community except for the explicitly labeled ensemble screens. Resistance is a pre-existing species trait, not within-species adaptation; no mutation, horizontal transfer, or changing resistance is modeled. The random-permutation experiment holds the rho multiset fixed but includes the associated resistance-growth cost, and its treated/control matching is therefore essential.
