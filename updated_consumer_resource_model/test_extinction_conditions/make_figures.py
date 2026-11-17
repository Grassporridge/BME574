"""Figures and summary tables for test_extinction_conditions/ (run run_scans.py first)."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_scans import CARRIERS, DOSES, EXPOSURES, FIXED, VMAXES

OUT = Path(__file__).resolve().parent / "extinction_results"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
mpl.rcParams.update({"figure.dpi": 120, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25,
                     "grid.linewidth": 0.6, "lines.linewidth": 2.0})
BLUE, ORANGE, GREEN, VERMILION, PURPLE, SKY, GREY = (
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#7F7F7F")
KCOL = dict(zip(CARRIERS, [BLUE, ORANGE, GREEN, VERMILION, PURPLE]))

grid = {f: pd.read_csv(OUT / f"grid_{f}.csv") for f in ("one_cell", "1e-12")}
sens = {f: pd.read_csv(OUT / f"all_sensitive_{f}.csv") for f in ("one_cell", "1e-12")}
VPOS = {v: i for i, v in enumerate(VMAXES)}          # vmax axis is categorical (it includes 0)


def vmax_axis(ax):
    ax.set_xticks(range(len(VMAXES)))
    ax.set_xticklabels([f"{v:g}" for v in VMAXES], rotation=90, fontsize=7)
    ax.axvline(VPOS[FIXED["vmax"]], color="k", lw=0.6, ls=":")
    ax.set_xlabel("vmax, degradation rate (µg/ml/h per OD of carrier)")


def exposure_axis(ax):
    ax.set_xlabel("drug exposure time (h)")
    ax.axvline(FIXED["exposure_h"], color="k", lw=0.6, ls=":")


def dose_axis(ax):
    ax.set_xscale("log"); ax.set_xlabel("dose (µg/ml)")
    ax.set_xticks([5, 10, 20, 50, 100, 200, 500]); ax.set_xticklabels(["5", "10", "20", "50", "100", "200", "500"])
    ax.axvline(FIXED["dose"], color="k", lw=0.6, ls=":")


def one_factor_row(axs, g, k, s=None, frac=False):
    """Three panels: extinction vs vmax, vs dose, vs exposure, others at FIXED (bold) + context."""
    col = "extinct_frac" if frac else "extinct"
    n_sus = 57 - k
    d, t, v = FIXED["dose"], FIXED["exposure_h"], FIXED["vmax"]
    gk = g[g.n_resistant == k]
    # (a) vmax
    ax = axs[0]
    for dd, c, lw in [(20.0, SKY, 1.1), (d, "k", 2.4), (140.0, ORANGE, 1.1), (500.0, VERMILION, 1.1)]:
        q = gk[(gk.dose == dd) & (gk.exposure_h == t)].sort_values("vmax")
        ax.plot(q.vmax.map(VPOS), q[col], "-o", ms=3, lw=lw, color=c,
                label=f"{dd:g} µg/ml" + (" (fixed)" if dd == d else ""))
    vmax_axis(ax); ax.legend(frameon=False, fontsize=7, title=f"exposure {t:g} h", title_fontsize=7)
    # (b) dose
    ax = axs[1]
    for vv, c, lw in [(0.0, VERMILION, 1.1), (30.0, ORANGE, 1.1), (100.0, GREEN, 1.1),
                      (300.0, SKY, 1.1), (1000.0, PURPLE, 1.1), (v, "k", 2.4)]:
        q = gk[(gk.vmax == vv) & (gk.exposure_h == t)].sort_values("dose")
        ax.plot(q.dose, q[col], "-o", ms=3, lw=lw, color=c,
                label=f"vmax {vv:g}" + (" (fixed)" if vv == v else ""))
    if s is not None:
        q = s[s.exposure_h == t].sort_values("dose")
        ax.plot(q.dose, q[col] if frac else q.extinct * n_sus / 57, "--", color=GREY, lw=1.2,
                label="no carrier at all" + ("" if frac else " (scaled)"))
    dose_axis(ax); ax.legend(frameon=False, fontsize=7, title=f"exposure {t:g} h", title_fontsize=7)
    # (c) exposure
    ax = axs[2]
    for vv, c, lw in [(0.0, VERMILION, 1.1), (100.0, GREEN, 1.1), (300.0, SKY, 1.1), (v, "k", 2.4)]:
        q = gk[(gk.vmax == vv) & (gk.dose == d)].sort_values("exposure_h")
        ax.plot(q.exposure_h, q[col], "-o", ms=3, lw=lw, color=c,
                label=f"vmax {vv:g}" + (" (fixed)" if vv == v else ""))
    if s is not None:
        q = s[s.dose == d].sort_values("exposure_h")
        ax.plot(q.exposure_h, q[col] if frac else q.extinct * n_sus / 57, "--", color=GREY, lw=1.2,
                label="no carrier at all")
    exposure_axis(ax); ax.legend(frameon=False, fontsize=7, title=f"dose {d:g} µg/ml", title_fontsize=7)
    for ax in axs:
        ax.set_ylim(-0.03 * (1 if frac else n_sus), (1.05 if frac else n_sus + 2))
    axs[0].set_ylabel("fraction of susceptible strains extinct" if frac
                      else f"susceptible strains extinct (of {n_sus})")


# ------------------------------------------------------------------ figure 1: one carrier, one factor
g, s = grid["one_cell"], sens["one_cell"]
fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
one_factor_row(axs, g, 1, s)
axs[0].set_title("vary degradation rate\n(dose 50 µg/ml, exposure 20 h)", fontsize=9)
axs[1].set_title("vary dose\n(vmax 10⁴, exposure 20 h; thin lines = slower degraders)", fontsize=9)
axs[2].set_title("vary exposure time\n(dose 50 µg/ml, vmax 10⁴; thin lines = slower degraders)", fontsize=9)
fig.suptitle("One resistant strain (k30) among 57: susceptible strains permanently lost "
             "(below one cell) by 48 h. Bold = one factor varied, the other two fixed (dotted).",
             y=1.04, fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "01_one_carrier_one_factor.png", dpi=150, bbox_inches="tight")


# ------------------------------------------------------------------ figure 2: one carrier, all three
def heat(ax, q, title, frac=False):
    col = "extinct_frac" if frac else "extinct"
    piv = q.pivot_table(index="vmax", columns="dose", values=col).reindex(index=VMAXES[::-1], columns=DOSES)
    im = ax.imshow(piv.values, aspect="auto", cmap="Reds", vmin=0, vmax=1 if frac else 56)
    ax.set_xticks(range(len(DOSES))); ax.set_xticklabels([f"{x:g}" for x in DOSES], rotation=90, fontsize=6)
    ax.set_yticks(range(len(VMAXES))); ax.set_yticklabels([f"{x:g}" for x in VMAXES[::-1]], fontsize=6)
    ax.grid(False); ax.set_title(title, fontsize=8.5)
    # outline the cells with zero extinction boundary: mark protected cells with a dot
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if piv.values[i, j] == 0:
                ax.plot(j, i, ".", color=GREEN, ms=2.2)
    return im


SHOW_T = [2.0, 6.0, 10.0, 12.0, 14.0, 20.0]
fig, axs = plt.subplots(2, 4, figsize=(19, 9.4), layout="constrained")
gk = g[g.n_resistant == 1]
for ax, t in zip([axs[0, 0], axs[0, 1], axs[0, 2], axs[1, 0], axs[1, 1], axs[1, 2]], SHOW_T):
    im = heat(ax, gk[gk.exposure_h == t], f"exposure {t:g} h")
    ax.set_xlabel("dose (µg/ml)"); ax.set_ylabel("vmax")
axs[0, 3].axis("off")
fig.colorbar(im, ax=axs[0, 3], fraction=0.5, label="susceptible strains extinct (of 56)\ngreen dot = none")
ax = axs[1, 3]
tcol = dict(zip(EXPOSURES, plt.cm.viridis(np.linspace(0, 0.95, len(EXPOSURES)))))
bound_rows = []
for t in EXPOSURES:
    q = gk[gk.exposure_h == t]
    mv = []
    for dd in DOSES:
        safe = q[(q.dose == dd) & (q.extinct == 0)].vmax
        mv.append(safe.min() if len(safe) else np.nan)
        bound_rows.append(dict(n_resistant=1, exposure_h=t, dose=dd, min_safe_vmax=mv[-1]))
    if t in (6.0, 10.0, 12.0, 14.0, 20.0):
        ax.plot(DOSES, np.maximum(mv, 1), "-o", ms=3, color=tcol[t], label=f"{t:g} h")
ax.set(xscale="log", yscale="log", xlabel="dose (µg/ml)", ylim=(0.8, 5e4),
       ylabel="lowest vmax with no extinction")
ax.set_title("degradation needed to lose no strain\n(1 = none needed; 14 h lies under 20 h)", fontsize=8.5)
ax.axhline(1e4, color="k", lw=0.6, ls=":"); ax.text(5.5, 1.25e4, "model default", fontsize=7)
ax.axhspan(3e4, 5e4, color=GREY, alpha=0.2, lw=0); ax.text(5.5, 3.4e4, "above physical ceiling", fontsize=7)
ax.legend(frameon=False, fontsize=7, title="exposure", title_fontsize=7)
fig.suptitle("One resistant strain: all three factors together")
fig.savefig(FIG / "02_one_carrier_all_three.png", dpi=150, bbox_inches="tight")

# ------------------------------------------------------------------ figure 3: repeat fig 1 per carrier count
fig, axs = plt.subplots(len(CARRIERS), 3, figsize=(15, 3.3 * len(CARRIERS)))
for row, k in zip(axs, CARRIERS):
    one_factor_row(row, g, k, None, frac=False)
    row[0].set_ylabel(f"{k} carrier{'s' if k > 1 else ''}:\nextinct (of {57 - k})")
    for ax in row[1:]:
        ax.set_ylabel("")
for ax, t in zip(axs[0], ["vary degradation rate (dose 50, 20 h)", "vary dose (vmax 10⁴, 20 h)",
                          "vary exposure (dose 50, vmax 10⁴)"]):
    ax.set_title(t, fontsize=9)
fig.suptitle("Figure 1 repeated for 1, 2, 5, 10 and 20 resistant strains (bold = factor varied, "
             "others fixed; thin = context)", y=1.0)
fig.tight_layout(); fig.savefig(FIG / "03_carriers_one_factor.png", dpi=150, bbox_inches="tight")

# ------------------------------------------------------------------ figure 4: all three, per carrier count
HT = [4.0, 10.0, 12.0, 14.0, 20.0]
fig, axs = plt.subplots(len(CARRIERS), len(HT), figsize=(19, 3.4 * len(CARRIERS)), layout="constrained")
for i, k in enumerate(CARRIERS):
    gk = g[g.n_resistant == k]
    for j, t in enumerate(HT):
        im = heat(axs[i, j], gk[gk.exposure_h == t], f"{k} carrier{'s' if k > 1 else ''}, exposure {t:g} h",
                  frac=True)
        if j == 0: axs[i, j].set_ylabel("vmax")
        if i == len(CARRIERS) - 1: axs[i, j].set_xlabel("dose (µg/ml)")
fig.colorbar(im, ax=axs, shrink=0.4, label="fraction of susceptible strains extinct; green dot = none")
fig.suptitle("All three factors, for each number of resistant strains")
fig.savefig(FIG / "04_carriers_all_three.png", dpi=150, bbox_inches="tight")

for k in CARRIERS[1:]:
    gk = g[g.n_resistant == k]
    for t in EXPOSURES:
        q = gk[gk.exposure_h == t]
        for dd in DOSES:
            safe = q[(q.dose == dd) & (q.extinct == 0)].vmax
            bound_rows.append(dict(n_resistant=k, exposure_h=t, dose=dd,
                                   min_safe_vmax=safe.min() if len(safe) else np.nan))
bound = pd.DataFrame(bound_rows)
bound.to_csv(OUT / "min_safe_vmax.csv", index=False)

# ------------------------------------------------------------------ figure 5: carrier-count summary
fig, axs = plt.subplots(1, 3, figsize=(15, 3.9))
for k in CARRIERS:
    q = bound[(bound.n_resistant == k) & (bound.exposure_h == 20.0)].sort_values("dose")
    axs[0].plot(q.dose, np.maximum(q.min_safe_vmax, 1), "-o", ms=3, color=KCOL[k], label=f"{k}")
    q = g[(g.n_resistant == k) & (g.vmax == 0) & (g.dose == 50.0)].sort_values("exposure_h")
    axs[1].plot(q.exposure_h, q.extinct_frac, "-o", ms=3, color=KCOL[k], label=f"{k}")
    q = g[(g.n_resistant == k) & (g.vmax == 100.0) & (g.dose == 200.0)].sort_values("exposure_h")
    axs[2].plot(q.exposure_h, q.extinct_frac, "-o", ms=3, color=KCOL[k], label=f"{k}")
q = s[s.dose == 50.0].sort_values("exposure_h")
axs[1].plot(q.exposure_h, q.extinct_frac, "--", color=GREY, lw=1.2, label="0 (no carrier)")
axs[0].set(xscale="log", yscale="log", xlabel="dose (µg/ml)", ylabel="lowest vmax per carrier with no extinction",
           ylim=(0.8, 5e4))
axs[0].set_title("degradation needed (exposure 20 h)\n1 = none needed", fontsize=9)
axs[1].set(xlabel="drug exposure time (h)", ylabel="fraction of susceptibles extinct", ylim=(-0.03, 1.05))
axs[1].set_title("non-degrading resistant strains (vmax 0), 50 µg/ml", fontsize=9)
axs[2].set(xlabel="drug exposure time (h)", ylim=(-0.03, 1.05))
axs[2].set_title("slow degraders (vmax 100), 200 µg/ml", fontsize=9)
for ax in axs:
    ax.legend(frameon=False, fontsize=7, title="resistant strains", title_fontsize=7)
fig.tight_layout(); fig.savefig(FIG / "05_carrier_count_summary.png", dpi=150, bbox_inches="tight")

# ------------------------------------------------------------------ figure 6: extinction threshold
fig, axs = plt.subplots(1, 3, figsize=(15, 3.9))
for fl, ls in [("one_cell", "-"), ("1e-12", "--")]:
    gg = grid[fl]
    for k in (1, 5, 20):
        q = gg[(gg.n_resistant == k) & (gg.exposure_h == 20.0)]
        mv = [q[(q.dose == dd) & (q.extinct == 0)].vmax.min() for dd in DOSES]
        axs[0].plot(DOSES, np.maximum(mv, 1), ls, marker="o", ms=3, color=KCOL[k],
                    label=f"{k} carrier{'s' if k > 1 else ''}, floor {'1 cell' if fl == 'one_cell' else '1e-12 OD'}")
        q = gg[(gg.n_resistant == k) & (gg.vmax == 0) & (gg.dose == 50.0)].sort_values("exposure_h")
        axs[1].plot(q.exposure_h, q.extinct_frac, ls, marker="o", ms=3, color=KCOL[k])
    q = sens[fl][sens[fl].exposure_h == 20.0].sort_values("dose")
    axs[2].plot(q.dose, q.extinct_frac, ls, marker="o", ms=3, color="k",
                label=f"floor {'1 cell (1e-9 OD)' if fl == 'one_cell' else '1e-12 OD'}")
axs[0].set(xscale="log", yscale="log", xlabel="dose (µg/ml)", ylabel="lowest vmax with no extinction",
           ylim=(0.8, 5e4)); axs[0].set_title("exposure 20 h; solid = 1 cell, dashed = 1e-12 OD", fontsize=9)
axs[0].legend(frameon=False, fontsize=6.5)
axs[1].set(xlabel="drug exposure time (h)", ylabel="fraction extinct", ylim=(-0.03, 1.05))
axs[1].set_title("vmax 0, 50 µg/ml (colours as left)", fontsize=9)
axs[2].set(xscale="log", xlabel="dose (µg/ml)", ylim=(-0.03, 1.05))
axs[2].set_title("no resistant strain at all, exposure 20 h", fontsize=9); axs[2].legend(frameon=False, fontsize=7)
fig.suptitle("Sensitivity to the extinction threshold", y=1.03)
fig.tight_layout(); fig.savefig(FIG / "06_threshold_sensitivity.png", dpi=150, bbox_inches="tight")

# ------------------------------------------------------------------ printed summary
pd.set_option("display.width", 250)
print("lowest safe vmax per carrier at exposure 20 h (one-cell floor):")
print(bound[bound.exposure_h == 20.0].pivot_table(index="n_resistant", columns="dose",
                                                  values="min_safe_vmax").to_string())
print("\nlowest safe vmax, 1 carrier, by exposure (rows) and dose (cols):")
print(bound[bound.n_resistant == 1].pivot_table(index="exposure_h", columns="dose",
                                                values="min_safe_vmax").to_string())
gg = grid["1e-12"]
rows = []
for k in CARRIERS:
    q = gg[(gg.n_resistant == k) & (gg.exposure_h == 20.0)]
    rows.append({dd: q[(q.dose == dd) & (q.extinct == 0)].vmax.min() for dd in DOSES})
print("\nsame, 1e-12 floor (rows = carriers):")
print(pd.DataFrame(rows, index=CARRIERS).to_string())
print("\nfraction extinct at vmax 0, 50 ug/ml, by exposure (rows = carriers):")
print(g[(g.vmax == 0) & (g.dose == 50.0)].pivot_table(index="n_resistant", columns="exposure_h",
                                                       values="extinct_frac").round(2).to_string())
print("\nall-sensitive community, extinct (of 57), by dose (rows) x exposure:")
print(s.pivot_table(index="dose", columns="exposure_h", values="extinct").to_string())
fixed = g[(g.dose == 50.0) & (g.exposure_h == 20.0) & (g.vmax == 10000.0)]
print("\nat the fixed protocol values (50 ug/ml, 20 h, vmax 1e4):")
print(fixed[["n_resistant", "extinct", "undetected_alive", "carrier_share", "clearance_h"]].to_string(index=False))
