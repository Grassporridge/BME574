# Pre-existing resistance and post-antibiotic community recovery

A minimal model of the barcoded *E. coli* community experiment: grow the community, treat it
with carbenicillin, dilute into fresh drug-free medium, let it recover, and measure how much
of the community came back.

Everything is in `resilience_model.py`. A notebook imports it, changes whatever it wants to
vary, and plots.

```python
from resilience_model import run, sweep, grid, BASELINE, check

check()                                   # once per session; 15/15 should pass

r = run()                                 # the baseline experiment
r.metrics                                 # one row of readouts
r.timeseries                              # per-strain abundance at 4, 6, 24, 48 h
r.protection                              # per-strain protection score at 24 h

r = run(growth_mean=0.2, dose=70)         # anything not named keeps its baseline value

df = sweep("n_resistant", [1, 2, 4, 8, 16, 32])
df = grid("gamma", [1, 2, 4], "vmax", [1e3, 1e4], cache="panelC.csv")
```

Units: **time in hours, biomass in OD, drug in µg/ml of carbenicillin.** Every drug number is
the number the bench protocol uses — `dose=50` is 50 µg/ml.

---

## 1. What the model does

**The community.** 57 strains, matching the barcoded Keio library. One of them (k30, index 29)
carries *bla*TEM-1 on pUC19; the rest are susceptible. `n_resistant` can make more than one
strain a carrier, for the resistant-fraction question.

**Growth.** All strains compete for one shared nutrient with Monod kinetics. Carrying capacity
is not a parameter — it is however much food the medium supplies.

**Killing.** β-lactam killing is proportional to each cell's *current* growth rate, so a
starving cell is tolerant without being resistant. This is the single most important feature of
the model and it is well established experimentally (Tuomanen et al. 1986, below).

**Resistance** is two things: the carrier is immune to the drug, and it degrades the drug out
of the shared medium with saturating (Michaelis–Menten) kinetics. The enzyme is expressed
constantly from a high-copy plasmid, so its cost is a constant reduction in the carrier's
growth rate, paid with or without drug.

**The protocol**, in hours:

| hours | step |
|---|---|
| 0 – 4 | pre-growth, all strains at equal biomass |
| 4 | carbenicillin added at `dose` |
| 4 – 24 | treatment |
| 24 | sample, then 1:100 transfer into fresh drug-free medium |
| 24 – 48 | recovery |

Samples at 4, 6, 24 and 48 h. The transfer is what makes loss irreversible: a strain below 100
cells cannot survive a 100-fold dilution, and a strain that reaches one cell (10⁻⁹ OD) is
removed and never returns.

### Equations

For strain *i*, nutrient *R*, drug *a*:

```
g_i  = mu_i · R / (K + R)                         growth
kill_i = gamma · x_i/(1+x_i) · g_i,  x_i = (a/IC50_i)^h      (0 for carriers)
dN_i/dt = (g_i − death − kill_i) · N_i
dR/dt   = − Σ N_i g_i                             1 OD of biomass consumes 1 unit of food
da/dt   = − vmax · B_R · a/(Km + a)               B_R = biomass of live carriers
```

`mu_i` is reduced by `resistance_cost` for carriers. Biomass is integrated in log space so
abundances spanning many orders of magnitude stay accurate.

---

## 2. Assumptions, and why

Each of these is a deliberate simplification. They are listed with the reason and, where it
matters, what would break if the assumption is wrong.

| Assumption | Why it is reasonable here | Risk if wrong |
|---|---|---|
| **Killing is proportional to growth rate** | Measured directly for β-lactams on *E. coli* across generation times from 0.7 to 12 h (Tuomanen 1986) | Low. This is the best-supported assumption in the model. |
| **One shared nutrient, no cross-feeding** | The strains are near-isogenic single-gene deletions of one parent eating the same defined medium | Moderate. Removes a source of strain-to-strain variation — see the open issue in §6. |
| **Strains differ only in growth rate and IC50** | Keio mutants differ in little else that matters over 48 h | Moderate, same reason. |
| **Carriers are immune to the drug** | *bla*TEM-1 on pUC19 raises the MIC 100–1000×, far above the 50–70 µg/ml used | Low at these doses. The model cannot represent killing the resistant member. |
| **Degradation saturates (Michaelis–Menten)** | TEM-1's Km for carbenicillin is 14 µM ≈ 5.3 µg/ml, so at 50–70 µg/ml the enzyme runs at 90–93% of maximum for most of the treatment | **High if removed.** First-order degradation erases the dose dependence (see §4). |
| **The enzyme cost is a constant growth penalty** | pUC19 is high-copy and *bla*TEM-1 is constitutive, so the carrier pays whether or not drug is present | Low; the cost is small anyway. |
| **No free enzyme pool from lysed carriers** | Carriers are immune, so they are not lysed by the drug | Low. |
| **No spontaneous drug decay** | Carbenicillin is stable over 48 h relative to enzymatic removal | Low. |
| **One constant death rate, fed or starving** | Starvation death only shapes the slow decline after food runs out, which no readout examines | Low. |
| **Deterministic, well mixed** | Keeps the model interpretable | **Known limit.** No replicate spread, no persister subpopulation, no spatial structure — which is where private protection matters most. |
| **A strain below 10⁻⁴ relative abundance is not counted** | Stands in for the barcode-sequencing detection limit | Replace with the real limit if known. |

### If needed, to be added

- **PA14 invasion** (manuscript Figure 3). This is the one manuscript readout the model does not
  produce. Adding it back means one extra population and one readout.
- **β-lactamase inhibitor** (clavulanate, planned Figure 4).
- **Two named media.** There is no "M9CA" or "LB" switch. A medium is just a growth rate and a
  food supply, so set `growth_mean` and `food` directly. Roughly: M9CA ≈ `growth_mean=0.55,
  food=1.0`; LB ≈ `growth_mean=1.40, food=2.5`. Deliberately *not* baked in — in the earlier
  version an extra 90× killing factor was hidden inside the "LB" preset, and it caused a
  conclusion to be misattributed to growth rate when it came from the killing rate.

---

## 3. Readouts

Everything the manuscript reports, plus two extras. `run()` returns all of them.

| Readout | In `metrics` | Definition |
|---|---|---|
| Absolute abundance | `timeseries.absolute` | OD per strain at 4, 6, 24, 48 h |
| Relative abundance | `timeseries.relative` | fraction of the community |
| Shannon diversity | `shannon_treated`, `shannon_control` | H = −Σ p ln p over **all** members, as in the manuscript |
| Richness | `richness_treated`, `richness_control` | strains above the detection limit |
| Biomass resilience | `biomass_resilience` | treated ÷ untreated biomass at 48 h |
| Shannon resilience | `shannon_resilience` | treated ÷ untreated H at 48 h |
| Compositional resilience | `compositional_resilience` | 1 − Bray–Curtis, treated against untreated |
| Protection score | `protection` (per strain) | log2[(treated/untreated with resistance) ÷ (treated/untreated without)] at 24 h |

Two extras not in the manuscript, useful for the modelling questions:

| Readout | Meaning |
|---|---|
| `resistant_share` | what fraction of the recovered community is the carrier(s) |
| `protected_*` | the same metrics computed over the **non-carrier** members only |

The second matters. The manuscript's Shannon includes k30, so a k30 takeover shows up as a loss
of diversity — which is correct, and is the manuscript's central caution. But when the question
is "did the resistance protect the others?", mixing the protector into the diversity number
confounds two things. Both are reported; `shannon_*` is the one that matches her figures.

Every readout is measured against a **matched untreated control**: the same strains, the same
schedule, no drug. `run()` produces four simulations — resistant and all-sensitive communities,
each treated and untreated — which is the manuscript's own design and is what the protection
score needs.

---

## 4. Parameters, baselines and tested ranges

`BASELINE` holds the defaults; `TESTED_RANGES` holds the spans below. These are the ranges over
which the model has been exercised and behaves sensibly — outside them the solver may struggle
or the biology stops making sense.

| Parameter | Baseline | Tested range | Notes |
|---|---|---|---|
| `growth_mean` | 0.55 /h | 0.02 – 1.8 | ~0.02–0.2 is the range described as natural; 0.55 ≈ M9CA, 1.40 ≈ LB |
| `growth_cv` | 0.08 | 0.0 – 0.6 | Keio near-isogenic ≈ 0.05–0.15; above ~0.3 the untreated control is itself uneven, so diversity *ratios* stop being interpretable |
| `food` | 1.0 OD-eq | 0.2 – 5.0 | sets carrying capacity; measurable as the plateau OD |
| `half_saturation` | 0.05 OD-eq | 0.01 – 0.5 | Monod constant |
| `death_rate` | 0.004 /h | 0.0005 – 0.02 | constant, fed or starving |
| `n_resistant` | 1 | 1 – 57 | the experiment uses 1 of 57 = 1.8% |
| `focal_multiplier` | 1.0 | 1 – 30 | k30's starting share is m/(56+m): 1× = 1.75%, 10× = 15.2% |
| `resistance_cost` | 0.05 | 0.0 – 0.5 | plasmid-borne enzyme realistically 0–0.15 |
| `dose` | 50 µg/ml | 5 – 500 | the protocol uses 50 and 70 |
| `gamma` | 2.5 | 0.3 – 30 | max killing as a multiple of current growth; below 1 it cannot drive a growing population down |
| `hill` | 2.0 | 1 – 4 | β-lactam dose–response is typically 1–2 |
| `susceptible_ic50` | 5.0 µg/ml | 1 – 20 | median MIC of the susceptible strains — **worth measuring** |
| `ic50_cv` | 0.70 | 0.1 – 1.2 | **the one parameter tuned rather than measured** (see §6) |
| `vmax` | 10 000 µg/ml/h/OD | 10 – 1e5 | not independently measurable; chosen so the baseline sits in the regime the experiment reports |
| `km` | 5.3 µg/ml | 0.5 – 50 | TEM-1 measured value; see §4 |
| `initial_biomass` | 0.005 OD | 1e-4 – 0.05 | |
| `dilution` | 100 | 10 – 1000 | |
| `pregrowth/treatment/recovery` | 4 / 20 / 24 h | — | set by the protocol |
| `cells_per_OD` | 1e9 | — | sets the one-cell extinction floor |
| `detection_rel` | 1e-4 | 1e-5 – 1e-3 | barcode detection limit |

**What the baseline produces:** the all-sensitive community is wiped out; the resistant community
recovers its biomass fully and keeps all 57 members, but the carrier grows to 44% of the
community and composition shifts (compositional resilience 0.57). Protection scores are strongly
positive. This is the manuscript's qualitative result.

---

## 5. Sanity checks

`check()` runs 15 checks and prints a table. All pass at the baseline. They cover: mass balance,
no negative state, dose 0 reproducing the untreated control, killing requiring growth, the
transfer being an exact 1:100 split, extinction being irreversible, carriers surviving a lethal
dose, more carriers leaving less drug, degradation saturating above Km and being first-order
below it, resilience of an untreated run against itself, degradation rather than competition
driving the rescue, Shannon bounded by ln(57), solver convergence, and the direction of the main
result holding across four seeds.

Run it once per session. It takes a few seconds and is not run on import.

---

## References

- Tuomanen E, Cozens R, Tosch W, Zak O, Tomasz A (1986). *The rate of killing of Escherichia coli
  by β-lactam antibiotics is strictly proportional to the rate of bacterial growth.*
  **Microbiology** 132(5):1297–1304.
  <https://www.microbiologyresearch.org/content/journal/micro/10.1099/00221287-132-5-1297>
