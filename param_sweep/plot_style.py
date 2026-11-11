"""
plot_style.py

Shared color choices for the sweep notebooks, following the dataviz skill's
rules even though these are matplotlib/seaborn static charts rather than an
HTML artifact:
  - sequential (magnitude: f_R, peak loss, diversity, times) -> a single
    perceptually-uniform ramp (viridis), never a rainbow map like jet/turbo.
  - diverging (signed quantities, e.g. resistant advantage or a criterion
    margin that crosses zero) -> a two-hue-plus-neutral-midpoint map (RdBu),
    centered explicitly at zero via a TwoSlopeNorm rather than relying on
    the data's own min/max to land the midpoint on zero by luck.
  - categorical (a handful of named outcome classes) -> the Okabe-Ito
    palette, a standard colorblind-safe qualitative set, assigned in a fixed
    order (never re-cycled if a category is dropped).
"""

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

SEQUENTIAL_CMAP = "viridis"

# Okabe-Ito qualitative palette, fixed assignment order.
CATEGORICAL = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "vermillion": "#D55E00",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "sky_blue": "#56B4E9",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
}

DIVERGING_CMAP = "RdBu_r"  # red = positive (e.g. resistant favored), blue = negative


def diverging_norm(values):
    """A TwoSlopeNorm centered at zero, robust to all-positive or all-negative
    inputs (falls back to a symmetric range around 0 in those edge cases)."""
    vmin = min(values.min(), -1e-9)
    vmax = max(values.max(), 1e-9)
    return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)


def new_fig(figsize=(6, 5)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return fig, ax
