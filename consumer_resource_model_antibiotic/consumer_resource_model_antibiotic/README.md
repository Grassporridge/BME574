# Antibiotic analysis workflow

This version of the model adds two controlled design features:

1. **Balanced resource-use categories.** The seven possible non-empty feeding categories are represented as evenly as possible. With 25 species, every category contains either 3 or 4 species.
2. **A fitness cost of resistance.** Resistance still raises IC50, but it now also reduces maximum resource-supported growth. The default maximum cost is 20% at resistance `rho=1`.

## Core resistance tradeoff

Each species has `rho_i in [0,1]` and

`IC50_i = IC50_min * (IC50_max / IC50_min)^rho_i`.

The growth multiplier is

`growth_multiplier_i = 1 - resistance_cost_max * rho_i`.

Therefore

`effective_mu_max_i = intrinsic_mu_max_i * growth_multiplier_i`.

With the default `resistance_cost_max = 0.20`:

- `rho=0` -> no growth penalty
- `rho=0.5` -> 10% maximum-growth penalty
- `rho=1` -> 20% maximum-growth penalty

The numerical value 0.20 is an exploratory model assumption, not a calibrated biological constant. It can be changed in `Config` or analyzed with `run_resistance_cost_analysis.py`.

## Balanced feeding categories

The seven categories are:

- R0
- M1
- R0+M1
- M2
- R0+M2
- M1+M2
- R0+M1+M2

The category sequence is repeated as evenly as possible and then shuffled with a fixed seed. Thus category counts differ by at most one while species identities remain reproducible for a given seed.

The generated `species_summary.csv` and control design CSV include both `food_category` and `food_category_code`.

## Install and run

From the project directory:

```powershell
python -m pip install -r requirements.txt
python run_all_antibiotic_analyses.py
```

To verify the implementation first:

```powershell
python check_model.py
```

To run one experiment at a time:

```powershell
python run_control.py
python run_dose_response.py
python run_decay_analysis.py
python run_resistance_selection.py
python run_resistance_cost_analysis.py
python run_growth_death_analysis.py
python run_resource_interaction.py
python run_sensitivity_analysis.py
```

## 1. `run_control.py`

Establishes the no-antibiotic reference.

Additional design outputs in this version:

- `00_species_design.png` shows the balanced food-category counts and the resistance-growth penalty.
- `control_species_design.csv` lists each species' food category, resistance, IC50, intrinsic mu_max, fitness cost, growth multiplier, and effective mu_max.

The control is important because resistance now has a fitness cost even when antibiotic is absent.

## 2. `run_dose_response.py`

Varies only initial antibiotic dose with a fixed bacterial community.

Default doses: `0, 0.1, 0.25, 0.5, 1, 2`.

Use this to measure how dose changes total biomass, survivors, diversity, and resistance enrichment while the resource-use categories and species traits remain fixed.

## 3. `run_decay_analysis.py`

Varies antibiotic persistence while holding initial dose fixed.

`kdeg = 0` means no degradation. For positive decay, half-life is `ln(2)/kdeg`.

## 4. `run_resistance_selection.py`

Tests whether higher-IC50 species are preferentially retained under antibiotic exposure relative to their matched no-antibiotic outcomes.

Because resistance now carries a growth cost, comparing treatment with the matched control is especially important.

## 5. `run_resistance_cost_analysis.py`

New analysis specifically for the resistance-growth tradeoff.

Default maximum fitness costs:

`0, 0.10, 0.20, 0.30, 0.40`.

Each cost is evaluated both with no antibiotic and with a default treatment dose `A0=1`.

Outputs include:

- final biomass, survivors, Shannon diversity, and weighted IC50 versus resistance cost
- the implemented growth-multiplier curves
- species-level antibiotic benefit versus resistance under different fitness costs

Example:

```powershell
python run_resistance_cost_analysis.py --costs 0 0.1 0.2 0.3 0.4 --dose 1
```

## 6. `run_growth_death_analysis.py`

Plots resource-supported growth, antibiotic killing, net per-capita growth, and biomass for representative low-, median-, and high-resistance species.

## 7. `run_resource_interaction.py`

Tests initial resource availability by antibiotic dose.

## 8. `run_sensitivity_analysis.py`

Exploratory global sensitivity screen. It now varies:

- initial antibiotic
- antibiotic decay rate
- maximum antibiotic kill rate
- Hill coefficient
- initial resource
- leakage fraction
- **maximum resistance fitness cost**

This is a rank-correlation screen, not a formal Sobol analysis.

## Recommended analysis order

1. Control and species design
2. Dose response
3. Antibiotic decay
4. Resistance selection
5. Resistance fitness-cost tradeoff
6. Growth-versus-death mechanism
7. Resource interaction
8. Global sensitivity

## Important interpretation limits

The parameter values are exploratory normalized values. The IC50 range, antibiotic dose, decay rate, kill rate, and resistance fitness cost are not fitted to a particular organism or drug. Biological conclusions should therefore be interpreted as model behavior until those parameters are calibrated to data.
