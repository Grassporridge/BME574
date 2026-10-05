
# BME574 parameter sweeps — Ma et al. 2024 antibiotic + inhibitor model

Setup: `pip install -r requirements.txt`, then open the notebooks from the `notebooks/` folder.

| Notebook | What it does | Run time |
|---|---|---|
| 00_model_verification | 34 checks of the model and engine against the paper | ~1 min |
| 01_run_sweeps | Runs every sweep and writes `results/tables/` (only needed to regenerate the tables) | ~25 min on 2 cores |
| 02_alpha_betamin_c | Resistant-strain traits | seconds |
| 03_gamma_ha_a0 | Antibiotic lysis and dose | seconds |
| 04_inhibitor_i_hi_db | Inhibitor dose and action | ~1 min |
| 05_kappab_da_phimax_xi | Degradation and recycling | ~1 min |
| 06_criterion_validation | Paper's selection criterion vs simulation | seconds |
| 07_summary | One-notebook summary of everything | seconds |
| 08_dose_by_trait | Dose plane × α, β_min, c: where a dose removes R and keeps S (recomputes its two tables if missing, ~25 min) | ~1 min |
| 09_fine_grid_low_degradation | Fine grid for low φ_max and κ_b, with d_a (recomputes its table if missing, ~1.5 min) | ~1 min |
| 10_killing_vs_protection | Is the relation between killing rate (γ) and R's protection linear? (recomputes its two tables if missing, ~1 min) | ~1 min |
| 11_gamma_vs_phimax | Is the relation between killing rate (γ) and degradation by living R (φ_max) linear? (recomputes its two tables if missing, ~30 s) | seconds |

Notebooks 02–11 read the saved tables, so they run without re-running 01.

Code: `src/model/ma2024_model.py` (equations, SI base parameters), `src/analysis/sweep_engine.py` (integration, endpoints, outcome labels), `src/analysis/grid_runner.py` (grids), `src/analysis/dose_window.py` (selective dose window), `src/analysis/kill_protection.py` (killing rate vs. protection), `src/analysis/kill_degradation.py` (killing rate vs. φ_max), `src/analysis/plotting.py` (figures).
