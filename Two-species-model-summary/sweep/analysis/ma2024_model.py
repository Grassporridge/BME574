"""
ma2024_model.py -- the five-state sensitive/resistant + beta-lactamase model of

    Ma HR, Xu HZ, Kim K, Anderson DJ, You L (2024). "Private benefit of
    beta-lactamase dictates selection dynamics of combination antibiotic
    treatment." Nature Communications 15:8337.

Equations are Supplementary Eqs. S1-S10 (identical to the Methods). Every
parameter value, initial condition and range below is taken from the paper's
Supplementary Information, Sections 1-2.

State (dimensionless):  n_s  sensitive cells      n_r  resistant cells
                        s    nutrient             a    antibiotic
                        b    free (extracellular) beta-lactamase

    dn_s/dt = (g - l) n_s                                           (S1)
    dn_r/dt = (alpha g - beta l) n_r                                (S2)
    ds/dt   = (xi l - g) n_s + (xi beta l - alpha g) n_r            (S3)
    da/dt   = -kappa_b b a - phi n_r a - d_a a                      (S4)
    db/dt   = beta l n_r - d_b iota b                               (S5)
    g = s/(1+s)                      l = gamma a^h_a/(1+a^h_a) g    (S6, S7)
    iota = i^h_i/(1+i^h_i)                                          (S8)
    beta = beta_min + c (1-beta_min) iota                           (S9)
    phi  = phi_max (1 - c iota)                                     (S10)

This file is the plain, linear-coordinate reference implementation. The sweeps
use src/analysis/sweep_engine.py, which integrates the same equations in
log-population coordinates and is checked against this file in notebook 00.
"""

from __future__ import annotations
from dataclasses import dataclass, fields
import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Params:
    """Base ("Initial Value") parameters, SI Section 2."""
    alpha: float = 0.95      # growth-rate factor of resistant cells (burden = 1 - alpha)
    beta_min: float = 0.9    # lysis-rate factor of resistant cells with no inhibitor
                             # (private benefit = 1 - beta_min)
    c: float = 0.7           # inhibitor's intracellular efficacy
    gamma: float = 1.38      # maximum lysis coefficient (l/g at saturating antibiotic)
    h_a: float = 3.0         # Hill coefficient of the antibiotic
    h_i: float = 2.0         # Hill coefficient of the inhibitor
    kappa_b: float = 0.35    # antibiotic degradation by free beta-lactamase
    d_a: float = 0.02        # basal antibiotic decay
    d_b: float = 1.0         # inhibitor-dependent inactivation of free beta-lactamase
    xi: float = 0.8          # nutrient recycled per unit of lysed biomass
    phi_max: float = 1.0     # antibiotic degradation by living resistant cells


DEFAULT_PARAMS = Params()
PARAM_NAMES = tuple(f.name for f in fields(Params))

# SI Section 2, "Randomized Range" column (uniform). d_a has no randomized range
# in the paper; the group sweeps it over 0-0.1 (0 to 5x its base value), which is
# the group's own choice and is flagged wherever it is used.
SI_RANDOMIZED_RANGES = {
    "alpha": (0.75, 1.0),
    "beta_min": (0.0, 1.0),
    "c": (0.0, 1.0),
    "gamma": (1.1, 1.4),
    "h_a": (1.0, 5.0),
    "h_i": (1.0, 5.0),
    "kappa_b": (0.0, 1.0),
    "d_b": (1.0, 10.0),
    "xi": (0.0, 1.0),
    "phi_max": (0.0, 5.0),
}
GROUP_CHOSEN_RANGES = {"d_a": (0.0, 0.10)}

# SI Section 1: dose ranges and initial conditions.
SI_DOSE_RANGES = {"a0": (1.0, 100.0), "i": (0.1, 10.0)}
SI_INITIAL_CONDITIONS = {"ns0": 0.2, "nr0": 0.2, "s0": 4.0, "b0": 0.0}

# The paper's Fig. 2C dose: the reference treatment for every sweep.
BASE_A0 = 2.15
BASE_I = 10.0


def default_y0(a0: float) -> np.ndarray:
    """[n_s, n_r, s, a, b] at t = 0 (SI Section 1)."""
    ic = SI_INITIAL_CONDITIONS
    return np.array([ic["ns0"], ic["nr0"], ic["s0"], float(a0), ic["b0"]])


def hill(x, h):
    """x^h / (1 + x^h), with x clipped at 0."""
    xh = np.power(np.maximum(x, 0.0), h)
    return xh / (1.0 + xh)


def growth_rate(s):
    return s / (1.0 + s)


def inhibitor_effect(i_dose: float, p: Params) -> float:
    return float(hill(i_dose, p.h_i))


def beta_of(iota: float, p: Params) -> float:
    return p.beta_min + p.c * (1.0 - p.beta_min) * iota


def phi_of(iota: float, p: Params) -> float:
    return p.phi_max * (1.0 - p.c * iota)


def rhs(tau: float, y: np.ndarray, p: Params, i_dose: float) -> np.ndarray:
    """Right-hand side of S1-S5 in the original (linear) coordinates."""
    ns, nr, s, a, b = y
    g = growth_rate(s)
    l = p.gamma * hill(a, p.h_a) * g
    iota = inhibitor_effect(i_dose, p)
    beta = beta_of(iota, p)
    phi = phi_of(iota, p)
    return np.array([
        (g - l) * ns,
        (p.alpha * g - beta * l) * nr,
        (p.xi * l - g) * ns + (p.xi * beta * l - p.alpha * g) * nr,
        -p.kappa_b * b * a - phi * nr * a - p.d_a * a,
        beta * l * nr - p.d_b * iota * b,
    ])


def simulate(p: Params, a0: float, i_dose: float, y0: np.ndarray | None = None,
             t_span=(0.0, 20.0), n_points: int = 400,
             method: str = "LSODA", rtol: float = 1e-8, atol: float = 1e-10):
    """Integrate S1-S5 for one dose in linear coordinates."""
    if y0 is None:
        y0 = default_y0(a0)
    else:
        y0 = np.array(y0, dtype=float).copy()
        y0[3] = a0
    sol = solve_ivp(rhs, t_span, y0, args=(p, i_dose),
                    t_eval=np.linspace(*t_span, n_points),
                    method=method, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")
    return sol


def criterion_general(p: Params, s: float, a: float, i_dose: float) -> bool:
    """Eq. 1 / S13: resistant cells enriched iff 1 - beta > (1 - alpha)/(l/g)."""
    ratio = p.gamma * float(hill(a, p.h_a))          # l/g
    if ratio <= 0:
        return False
    beta = beta_of(inhibitor_effect(i_dose, p), p)
    return (1.0 - beta) > (1.0 - p.alpha) / ratio


def criterion_simplified(p: Params) -> bool:
    """Eq. 2 (saturating antibiotic and inhibitor):
    (1 - c)(1 - beta_min) > (1 - alpha)/gamma."""
    return (1.0 - p.c) * (1.0 - p.beta_min) > (1.0 - p.alpha) / p.gamma
