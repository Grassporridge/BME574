"""Matched initial-composition experiment: equal, fast/slow enriched, ten random mixtures."""
from pathlib import Path
from dataclasses import replace
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import model as model_module
import analysis_tools as analysis_module
from model import Config
def display(x):
    print(x.to_string() if hasattr(x,"to_string") else x)
def run(config=Config(),output='initial_composition_experiment'):
    # Use the aligned "both" model: variable affinity + resistance growth cost.
    both_config = replace(
        config,
        variable_affinity=True,
        resistance_cost_max=.20,
        structured_community=False,
    )
    
    template = model_module.generate(both_config)
    traits = template["table"].copy()
    species = traits["species"].to_numpy()
    growth = traits["mu_max"].to_numpy()
    n = len(species)
    total_initial_biomass = both_config.initial_biomass
    
    # Starting percentages:
    # Equal = 4% each for 25 species.
    # Growth-biased communities retain every species, but favor fast/slow growers.
    # Random communities use different lognormal abundance mixtures.
    rng = np.random.default_rng(2026)
    growth_standardized = (
        (growth - growth.mean()) / max(growth.std(), 1e-12)
    )
    
    starting_shares = {
        "Equal": np.ones(n),
        "Fast enriched": np.exp(1.5 * growth_standardized),
        "Slow enriched": np.exp(-1.5 * growth_standardized),
    }
    
    for k in range(10):
        starting_shares[f"Random {k + 1:02d}"] = rng.lognormal(
            mean=0, sigma=1.5, size=n
        )
    
    starting_shares = {
        name: weights / weights.sum()
        for name, weights in starting_shares.items()
    }
    
    results = {}
    summary_rows = []
    species_rows = []
    
    for name, shares in starting_shares.items():
        # Same seed gives exactly the same traits in every experiment.
        community = model_module.generate(both_config)
    
        # Only change initial biomass allocation.
        initial = total_initial_biomass * shares
        community["initial"] = initial.copy()
        community["table"]["initial_biomass"] = initial
    
        result = model_module.simulate(community)
        diversity = analysis_module.diversity(result)
    
        final = result["N"][-1]
        final_total = final.sum()
        final_shares = (
            final / final_total
            if final_total > 0
            else np.full(n, np.nan)
        )
        survived = np.isnan(result["extinction"])
    
        results[name] = {
            "initial_shares": shares,
            "final_shares": final_shares,
            "survived": survived,
            "simulation": result,
            "diversity": diversity,
        }
    
        summary_rows.append({
            "starting_community": name,
            "initial_Shannon": diversity["shannon"].iloc[0],
            "final_Shannon": diversity["shannon"].iloc[-1],
            "survivors": int(survived.sum()),
            "final_total_biomass": final_total,
            "initial_weighted_growth": shares @ growth,
            "final_weighted_growth": (
                final_shares @ growth if final_total > 0 else np.nan
            ),
            "initial_weighted_resistance":
                shares @ traits["resistance"].to_numpy(),
            "final_weighted_resistance": (
                final_shares @ traits["resistance"].to_numpy()
                if final_total > 0 else np.nan
            ),
            "initial_weighted_resistance_cost":
                shares @ traits["resistance_fitness_cost"].to_numpy(),
            "final_weighted_resistance_cost": (
                final_shares @ traits["resistance_fitness_cost"].to_numpy()
                if final_total > 0 else np.nan
            ),
        })
    
        species_table = traits.copy()
        species_table["starting_community"] = name
        species_table["initial_percent"] = 100 * shares
        species_table["final_percent"] = 100 * final_shares
        species_table["final_biomass"] = final
        species_table["survived"] = survived
        species_table["extinction_time"] = result["extinction"]
        species_rows.append(species_table)
    
    summary = pd.DataFrame(summary_rows).set_index("starting_community")
    species_outcomes = pd.concat(species_rows, ignore_index=True)
    
    display(summary.round(4))
    
    # Save tables and figures.
    save_folder = Path(output)
    save_folder.mkdir(exist_ok=True,parents=True)
    summary.to_csv(save_folder / "community_summary.csv")
    species_outcomes.to_csv(
        save_folder / "species_traits_and_outcomes.csv", index=False
    )
    
    def caption(text):
        print("\n" + text + "\n")
    
    def finish_figure(fig, filename):
        fig.tight_layout()
        fig.savefig(save_folder / filename, dpi=180, bbox_inches="tight")
        plt.close(fig)
    
    names = list(results)
    x = np.arange(len(names))
    
    # 1. Shannon trajectories and final outcomes.
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    
    for name in names:
        d = results[name]["diversity"]
        highlighted = name in ["Equal", "Fast enriched", "Slow enriched"]
        axes[0].plot(
            d["time"], d["shannon"],
            linewidth=2 if highlighted else 0.8,
            alpha=1 if highlighted else 0.35,
            label=name if highlighted else None,
        )
    
    axes[0].set(xlabel="Time", ylabel="Shannon diversity H")
    axes[0].legend()
    
    axes[1].bar(x, summary["final_Shannon"])
    axes[1].set(ylabel="Final Shannon diversity")
    
    axes[2].bar(x, summary["survivors"])
    axes[2].set(ylabel="Species surviving at final time", ylim=(0, n + 1))
    
    for ax in axes[1:]:
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=65, ha="right")
    
    finish_figure(fig, "01_shannon_and_survival.png")
    caption(
        "Figure 1. Starting proportions can alter Shannon diversity and final composition; zero mortality preserves species. "
        "All experiments use identical traits, total starting biomass, and food. "
        "Unequal communities already have lower initial Shannon, so interpret "
        "the trajectories together with their starting values."
    )
    
    # 2. Starting and final community composition.
    initial_matrix = np.array([
        results[name]["initial_shares"] for name in names
    ])
    final_matrix = np.array([
        results[name]["final_shares"] for name in names
    ])
    
    colors = plt.cm.turbo(np.linspace(0, 1, n))
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    
    for ax, matrix, title in zip(
        axes,
        [initial_matrix, final_matrix],
        ["Starting composition", "Final composition"],
    ):
        bottom = np.zeros(len(names))
        for i in range(n):
            values = 100 * matrix[:, i]
            ax.bar(x, values, bottom=bottom, color=colors[i],
                   label=species[i])
            bottom += np.nan_to_num(values)
        ax.set(title=title, ylabel="Living biomass (%)", ylim=(0, 100))
    
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(names, rotation=65, ha="right")
    fig.legend(
        handles=axes[0].get_legend_handles_labels()[0],
        labels=species.tolist(),
        loc="upper center", bbox_to_anchor=(0.5, 1.04),
        ncol=13, fontsize=8,
    )
    finish_figure(fig, "02_initial_final_composition.png")
    caption(
        "Figure 2. Each color represents the same species in both panels. "
        "Compare which starting mixtures lead to different final dominants. "
        "Percentages describe the remaining living biomass, not its absolute "
        "amount; dominance reflects relative growth, not mortality in this control."
    )
    
    # 3. Biomass-weighted community traits.
    trait_columns = [
        ("growth", "Maximum growth rate"),
        ("resistance", "Resistance score"),
        ("resistance_cost", "Resistance fitness cost"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    for ax, (column, label) in zip(axes, trait_columns):
        ax.bar(
            x - 0.18, summary[f"initial_weighted_{column}"],
            width=0.36, label="Initial",
        )
        ax.bar(
            x + 0.18, summary[f"final_weighted_{column}"],
            width=0.36, label="Final",
        )
        ax.set(ylabel=f"Biomass-weighted {label.lower()}")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=65, ha="right")
    
    axes[0].legend()
    finish_figure(fig, "03_community_trait_enrichment.png")
    caption(
        "Figure 3. Community traits are weighted by species' biomass shares. "
        "A higher final weighted growth capacity indicates enrichment in fast "
        "growers, even if food depletion prevents them from growing at that "
        "rate. Resistance scores and growth costs are fixed traits; "
        "enrichment reflects resource competition without antibiotic killing."
    )
    
    # 4. Which species survived? Order species by growth capacity.
    order = np.argsort(growth)
    survival_matrix = np.array([
        results[name]["survived"] for name in names
    ], dtype=int)[:, order]
    
    fig, ax = plt.subplots(figsize=(13, 7))
    image = ax.imshow(
        survival_matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1
    )
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels(
        [f"{species[i]}\n{growth[i]:.2f}" for i in order],
        rotation=90, fontsize=8,
    )
    ax.set_yticks(x)
    ax.set_yticklabels(names)
    ax.set(
        xlabel="Species ordered from slow to fast; number = maximum growth",
        ylabel="Starting community",
    )
    bar = fig.colorbar(image, ax=ax, ticks=[0, 1])
    bar.ax.set_yticklabels(["Extinct", "Survived"])
    finish_figure(fig, "04_species_survival.png")
    caption(
        "Figure 4. Dark cells indicate species above the extinction cutoff "
        "at the final time. Species are ordered by maximum growth capacity. "
        "All species persist in this mortality-free control; "
        "the panel is a consistency check rather than evidence of killing."
    )
    
    # 5. Are fast or slow growers enriched overall?
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    axes[0].scatter(
        summary["initial_weighted_growth"],
        summary["final_Shannon"], s=60,
    )
    axes[0].set(
        xlabel="Initial biomass-weighted maximum growth",
        ylabel="Final Shannon diversity",
    )
    
    axes[1].scatter(
        summary["initial_weighted_growth"],
        summary["final_weighted_growth"], s=60,
    )
    low = min(
        summary["initial_weighted_growth"].min(),
        summary["final_weighted_growth"].min(),
    )
    high = max(
        summary["initial_weighted_growth"].max(),
        summary["final_weighted_growth"].max(),
    )
    axes[1].plot([low, high], [low, high], "k--", linewidth=1)
    axes[1].set(
        xlabel="Initial biomass-weighted maximum growth",
        ylabel="Final biomass-weighted maximum growth",
    )
    
    for name in ["Equal", "Fast enriched", "Slow enriched"]:
        row = summary.loc[name]
        axes[0].annotate(
            name, (row["initial_weighted_growth"], row["final_Shannon"]),
            xytext=(4, 4), textcoords="offset points", fontsize=8,
        )
        axes[1].annotate(
            name,
            (row["initial_weighted_growth"], row["final_weighted_growth"]),
            xytext=(4, 4), textcoords="offset points", fontsize=8,
        )
    
    finish_figure(fig, "05_growth_enrichment_and_diversity.png")
    caption(
        "Figure 5. Left: starting growth-trait enrichment versus final Shannon. "
        "Right: points above the dashed line indicate enrichment in faster "
        "growers over the batch; points below indicate enrichment in slower "
        "growers. These comparisons use one fixed set of species, so they "
        "describe this community rather than a universal ecological rule."
    )
    
    print(f"Saved tables and figures to: {save_folder.resolve()}")

if __name__=='__main__':
    run()
