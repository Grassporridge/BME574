"""
ma2024_model.py

Reproduction of the sensitive/resistant + beta-lactamase (Bla) ODE model from:

  Ma HR, Xu HZ, Kim K, Anderson DJ, You L (2024).
  "Private benefit of beta-lactamase dictates selection dynamics of
  combination antibiotic treatment." Nature Communications 15:8337.
  https://doi.org/10.1038/s41467-024-52711-w

The five-state ODE system (Supplementary Eqs. S1-S10, identical to the
Methods section) and the two selection criteria (Eqs. 1-2 / S13) are
implemented below.

ALL parameter values, initial conditions and randomization ranges are taken
from the paper's Supplementary Information, Section 2 ("Parameters for the
base model", table of "Initial Value" and "Randomized Range") and Section 1
(initial conditions n_s(0) = n_r(0) = 0.2, s(0) = 4, 1 < a(0) < 100, b(0) = 0).
With these values the Fig. 2C dose (a0 = 2.15, i = 10) reproduces the
published time course: both populations dip, then regrow to S ~ 2.2,
R ~ 1.7 by t ~ 15-20, with f_R falling from 0.5 to ~0.43.

Usage:
    python ma2024_model.py
This runs two things and saves plots to the current directory:
  1. A single time-course simulation (analogous to Fig. 2C), plus the
     resistant-fraction inset (analogous to the Fig. 2C inset).
  2. A randomized-parameter check of the simplified selection criterion
     (Eq. 2) against simulated selection outcomes (analogous to Fig. 3C),
     to sanity-check that the code and the criterion agree.
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass, replace
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FIGURE_DIR = PROJECT_ROOT / "results" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Parameters
# ---------------------------------------------------------------------------

@dataclass
class Params:
    # Base ("Initial Value") parameters, Supplementary Information Section 2.
    # Randomized ranges used by the paper for simulated strains are in
    # SI_RANDOMIZED_RANGES below.
    alpha: float = 0.95       # growth-rate factor of resistant cells (burden = 1 - alpha)
    beta_min: float = 0.9     # baseline lysis-rate factor of resistant cells (private benefit = 1 - beta_min)
    c: float = 0.7            # inhibitor's intracellular efficacy
    gamma: float = 1.38       # maximum lysis coefficient (l/g at saturating antibiotic)
    h_a: float = 3.0          # Hill coefficient, antibiotic effect on lysis
    h_i: float = 2.0          # Hill coefficient, inhibitor effect
    kappa_b: float = 0.35     # antibiotic degradation rate by free Bla
    d_a: float = 0.02         # basal antibiotic decay rate
    d_b: float = 1.0          # inhibitor-dependent inactivation rate of free Bla
    xi: float = 0.8           # nutrient recycled per unit lysed biomass
    phi_max: float = 1.0      # antibiotic degradation rate by living resistant cells


# Supplementary Information Section 2, "Randomized Range" column (uniform).
# d_a has no randomized range in the paper (held at 0.02); the inhibitor dose
# i is a treatment variable spanning 0.1-10, not a strain parameter.
SI_RANDOMIZED_RANGES = {
    "alpha": (0.75, 1.0),
    "beta_min": (0.0, 1.0),
    "xi": (0.0, 1.0),
    "kappa_b": (0.0, 1.0),
    "d_b": (1.0, 10.0),
    "gamma": (1.1, 1.4),
    "h_a": (1.0, 5.0),
    "h_i": (1.0, 5.0),
    "phi_max": (0.0, 5.0),
    "c": (0.0, 1.0),
}

# Supplementary Information Section 1: dose ranges and initial conditions.
SI_DOSE_RANGES = {"a0": (1.0, 100.0), "i": (0.1, 10.0)}
SI_INITIAL_CONDITIONS = {"ns0": 0.2, "nr0": 0.2, "s0": 4.0, "b0": 0.0}


def default_y0(a0: float) -> np.ndarray:
    """[n_s, n_r, s, a, b] at t = 0, per SI Section 1."""
    ic = SI_INITIAL_CONDITIONS
    return np.array([ic["ns0"], ic["nr0"], ic["s0"], a0, ic["b0"]])


DEFAULT_PARAMS = Params()


# ---------------------------------------------------------------------------
# 2. Closure relations (paper's Eqs. M6-M10 in the group doc)
# ---------------------------------------------------------------------------

def hill(x: np.ndarray, h: float) -> np.ndarray:
    """Standard Hill function x^h / (1 + x^h), used for both the antibiotic
    and inhibitor dose-response terms."""
    xh = np.power(np.maximum(x, 0.0), h)
    return xh / (1.0 + xh)


def growth_rate(s: float) -> float:
    """g = s / (1+s)  -- Monod-type growth on the shared nutrient pool."""
    return s / (1.0 + s)


def lysis_rate(s: float, a: float, p: Params) -> float:
    """l = gamma * hill(a, h_a) * g(s)."""
    return p.gamma * hill(a, p.h_a) * growth_rate(s)


def inhibitor_effect(i: float, p: Params) -> float:
    """iota = hill(i, h_i)."""
    return hill(i, p.h_i)


def beta_of(iota: float, p: Params) -> float:
    """beta = beta_min + c*(1-beta_min)*iota
    (private-benefit factor for the resistant population; beta -> 1 as the
    inhibitor saturates and c -> 1, meaning private protection is fully
    suppressed)."""
    return p.beta_min + p.c * (1.0 - p.beta_min) * iota


def phi_of(iota: float, p: Params) -> float:
    """phi = phi_max * (1 - c*iota)
    (antibiotic-degradation rate by living resistant cells; phi -> 0 as the
    inhibitor saturates and c -> 1)."""
    return p.phi_max * (1.0 - p.c * iota)


# ---------------------------------------------------------------------------
# 3. The five-ODE system (paper's Eqs. M1-M5)
# ---------------------------------------------------------------------------

def rhs(tau: float, y: np.ndarray, p: Params, i_dose: float) -> np.ndarray:
    ns, nr, s, a, b = y
    g = growth_rate(s)
    l = lysis_rate(s, a, p)
    iota = inhibitor_effect(i_dose, p)
    beta = beta_of(iota, p)
    phi = phi_of(iota, p)

    dns = (g - l) * ns
    dnr = (p.alpha * g - beta * l) * nr
    ds = (p.xi * l - g) * ns + (p.xi * beta * l - p.alpha * g) * nr
    da = -p.kappa_b * b * a - phi * nr * a - p.d_a * a
    db = beta * l * nr - p.d_b * iota * b

    return np.array([dns, dnr, ds, da, db])


def simulate(
    p: Params,
    a0: float,
    i_dose: float,
    y0: np.ndarray | None = None,
    t_span: tuple[float, float] = (0.0, 20.0),
    n_points: int = 400,
):
    """Integrate the ODE system for a single constant dose (a0, i_dose)."""
    if y0 is None:
        y0 = default_y0(a0)   # SI initial conditions
    else:
        y0 = y0.copy()
        y0[3] = a0

    t_eval = np.linspace(*t_span, n_points)
    sol = solve_ivp(
        rhs, t_span, y0, args=(p, i_dose), t_eval=t_eval,
        method="LSODA", rtol=1e-8, atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")
    return sol


# ---------------------------------------------------------------------------
# 4. Selection criteria (paper's Eqs. 1 and 2)
# ---------------------------------------------------------------------------

def criterion_general(p: Params, s: float, a: float, i_dose: float) -> bool:
    """Eq. 1: resistant cells enriched iff 1-beta > (1-alpha)/(l/g).
    l/g is evaluated at the given (s, a) state."""
    g = growth_rate(s)
    l = lysis_rate(s, a, p)
    if g == 0:
        return False
    iota = inhibitor_effect(i_dose, p)
    beta = beta_of(iota, p)
    return (1.0 - beta) > (1.0 - p.alpha) / (l / g)


def criterion_simplified(p: Params) -> bool:
    """Eq. 2, the saturating-dose limit: (1-c)(1-beta_min) > (1-alpha)/gamma."""
    return (1.0 - p.c) * (1.0 - p.beta_min) > (1.0 - p.alpha) / p.gamma


# ---------------------------------------------------------------------------
# 5. Reproduction 1: a single time course, analogous to Fig. 2C
# ---------------------------------------------------------------------------

def reproduce_fig2c(p: Params = DEFAULT_PARAMS, a0: float = 2.15, i_dose: float = 10.0):
    sol = simulate(p, a0=a0, i_dose=i_dose, t_span=(0.0, 20.0))  # Fig. 2C x-axis: 0-20
    ns, nr, s, a, b = sol.y
    fr = nr / (ns + nr)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.plot(sol.t, ns, label="S (sensitive)", color="tab:green")
    ax.plot(sol.t, nr, label="R (resistant)", color="tab:gray")
    ax.plot(sol.t, a, label="A (antibiotic)", color="tab:red", linestyle="--")
    ax.plot(sol.t, b, label="B (free Bla)", color="black", linestyle=":")
    ax.set_xlabel("Time (dimensionless, tau)")
    ax.set_ylabel("Concentration / density")
    ax.set_title(f"Time course (a0={a0}, i={i_dose})")
    ax.legend(fontsize=8)

    ax2 = axes[1]
    ax2.plot(sol.t, fr, color="tab:purple")
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("Time (dimensionless, tau)")
    ax2.set_ylabel("Resistant fraction f_R")
    ax2.set_title("Resistant fraction over time")

    fig.tight_layout()
    output = FIGURE_DIR / "fig2c_reproduction.png"
    fig.savefig(output, dpi=150)
    print(f"Saved {output}")
    return sol


# ---------------------------------------------------------------------------
# 6. Reproduction 2: randomized-parameter check of the criterion, Fig. 3C style
# ---------------------------------------------------------------------------

def resistant_advantage(sol, p: Params, i_dose: float) -> float:
    """Integral of d/dtau ln(nr/ns) = -(1-alpha)*g + (1-beta)*l over the run.

    At a saturating, near-lethal dose the raw densities ns, nr underflow to
    numerical zero, so their ratio becomes float noise. The log-ratio's
    *rate of change* only depends on g(t) and l(t) (via s(t), a(t)), which
    stay well-behaved even as the populations vanish, so integrating that
    rate gives a robust answer to "did the resistant share grow or shrink
    over the run" without ever dividing two near-zero numbers.
    """
    s_t, a_t = sol.y[2], sol.y[3]
    iota = inhibitor_effect(i_dose, p)
    beta = beta_of(iota, p)
    g_t = growth_rate(s_t)
    l_t = p.gamma * hill(a_t, p.h_a) * g_t
    rate = -(1.0 - p.alpha) * g_t + (1.0 - beta) * l_t
    trapezoid = getattr(np, "trapezoid", None) or np.trapz
    return float(trapezoid(rate, sol.t))


def reproduce_fig3c(n_strains: int = 300, seed: int = 0,
                    a0: float = 100.0, i_dose: float = 10.0):
    rng = np.random.default_rng(seed)

    # Randomize the three parameters the paper found to matter most, within the
    # SI's randomized ranges, holding everything else at the SI base values.
    alphas = rng.uniform(*SI_RANDOMIZED_RANGES["alpha"], n_strains)
    beta_mins = rng.uniform(*SI_RANDOMIZED_RANGES["beta_min"], n_strains)
    cs = rng.uniform(*SI_RANDOMIZED_RANGES["c"], n_strains)

    x_vals = []  # (1-alpha)/gamma
    y_vals = []  # (1-c)(1-beta_min)
    outcome = []  # simulated: True if resistant cells had a net log-ratio advantage

    for alpha, beta_min, c in zip(alphas, beta_mins, cs):
        p = replace(DEFAULT_PARAMS, alpha=alpha, beta_min=beta_min, c=c)
        # Highest doses in the paper's range (a0 = 100, i = 10), the closest the
        # paper's dose grid gets to the saturating limit Eq. 2 assumes. With
        # d_a = 0.02 the antibiotic takes a few hundred time units to clear, so
        # the run is long enough for the populations to regrow and settle.
        sol = simulate(p, a0=a0, i_dose=i_dose, t_span=(0.0, 600.0), n_points=1200)
        advantage = resistant_advantage(sol, p, i_dose=i_dose)
        outcome.append(advantage > 0)

        x_vals.append((1.0 - alpha) / p.gamma)
        y_vals.append((1.0 - c) * (1.0 - beta_min))

    x_vals = np.array(x_vals)
    y_vals = np.array(y_vals)
    outcome = np.array(outcome)

    # Criterion prediction: resistant-favored iff y > x (Eq. 2).
    predicted = y_vals > x_vals
    agreement = np.mean(predicted == outcome)
    false_neg = np.sum(predicted & ~outcome)   # criterion says R, sim says S
    false_pos = np.sum(~predicted & outcome)   # criterion says S, sim says R
    print(f"Criterion (Eq. 2) vs simulated outcome agreement: {agreement:.1%} "
          f"over {n_strains} randomized strains")
    print(f"  'criterion says resistant, sim says sensitive': {false_neg} cases")
    print(f"  'criterion says sensitive, sim says resistant': {false_pos} cases")

    fig, ax = plt.subplots(figsize=(5.5, 5))
    colors = np.where(outcome, "tab:red", "tab:blue")
    ax.scatter(x_vals, y_vals, c=colors, s=18, alpha=0.7,
               label=None)
    lims = [0, max(x_vals.max(), y_vals.max()) * 1.05]
    ax.plot(lims, lims, "k-", linewidth=1, label="y = x (criterion boundary)")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel(r"$(1-\alpha)/\gamma$")
    ax.set_ylabel(r"$(1-c)(1-\beta_{min})$")
    ax.set_title(f"Criterion check ({agreement:.0%} agreement)\n"
                 "red = resistant fraction grew, blue = it shrank")
    ax.legend(fontsize=8)
    fig.tight_layout()
    output = FIGURE_DIR / "fig3c_reproduction.png"
    fig.savefig(output, dpi=150)
    print(f"Saved {output}")
    return x_vals, y_vals, outcome


# ---------------------------------------------------------------------------
# 7. Community diversity readout (new addition for the group's revised goal)
# ---------------------------------------------------------------------------

def shannon_diversity(*densities: float) -> float:
    """Shannon diversity index H = -sum(p_i * log(p_i)) over the given
    population densities (any number of taxa). Returns 0 if only one taxon
    is present, and 0 if total density is 0."""
    total = sum(densities)
    if total <= 0:
        return 0.0
    fracs = [d / total for d in densities if d > 0]
    return float(-sum(f * np.log(f) for f in fracs))


def diversity_over_time(sol) -> np.ndarray:
    """Shannon diversity of the S/R pair at each saved time point of a
    simulate() solution. With only two taxa this ranges from 0 (one taxon
    dominant) to log(2) approx 0.693 (equal split)."""
    ns, nr = sol.y[0], sol.y[1]
    return np.array([shannon_diversity(s, r) for s, r in zip(ns, nr)])


if __name__ == "__main__":
    print("Reproducing Fig. 2C-style time course...")
    sol = reproduce_fig2c()

    print("\nReproducing Fig. 3C-style criterion check...")
    reproduce_fig3c()

    print("\nShannon diversity at start and end of the Fig. 2C run:")
    div = diversity_over_time(sol)
    print(f"  H(t=0)  = {div[0]:.3f}")
    print(f"  H(t=end) = {div[-1]:.3f}")
