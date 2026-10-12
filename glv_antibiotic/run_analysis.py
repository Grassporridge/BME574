"""Run the gLV antibiotic analysis suite and build a PDF report."""
from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from model import Config, metrics, sample_community, simulate, species_outcomes, stability_diagnostics

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
FIG = OUT / "figures"
REPORTS = OUT / "reports"


def savefig(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_run(c, p, label, **overrides):
    result = simulate(c, p, **overrides)
    metric_args = {
        "initial_biomass": overrides.get("initial_biomass"),
        "dose": overrides.get("antibiotic_initial"),
        "decay": overrides.get("antibiotic_decay"),
        "kmax": overrides.get("antibiotic_kmax"),
        "hill": overrides.get("antibiotic_hill"),
        "cost": overrides.get("resistance_cost_max"),
    }
    row = metrics(c, p, result, **metric_args)
    row["scenario"] = label
    return result, row


def assign_resistance(community, rho_values):
    """Return a community copy with the same rho values assigned to new species."""
    assigned = dict(community)
    rho_values = np.asarray(rho_values, dtype=float).copy()
    lo = community["ic50_min"]
    hi = community["ic50_max"]
    new_ic50 = lo * (hi / lo) ** rho_values
    table = community["species"].copy(deep=True)
    table["resistance_rho"] = rho_values
    table["IC50"] = new_ic50
    table["fitness_cost_at_max"] = community["resistance_cost_max"] * rho_values
    assigned["resistance"] = rho_values
    assigned["IC50"] = new_ic50
    assigned["species"] = table
    return assigned


def resistance_assignment_sweep(c, community, n_permutations=200):
    """Compare rho-to-species assignments while holding ecology and rho multiset fixed."""
    rho_values = np.sort(community["resistance"].copy())
    growth_order = np.argsort(community["r"])
    assignments = [("seeded_assignment", community["resistance"].copy())]
    high_on_high = np.empty(c.n_species)
    high_on_high[growth_order] = rho_values
    assignments.append(("high_rho_on_high_growth", high_on_high))
    high_on_low = np.empty(c.n_species)
    high_on_low[growth_order] = rho_values[::-1]
    assignments.append(("high_rho_on_low_growth", high_on_low))
    abundance_order = np.argsort(community["x_star"])
    high_on_high_abundance = np.empty(c.n_species)
    high_on_high_abundance[abundance_order] = rho_values
    assignments.append(("high_rho_on_high_baseline_abundance", high_on_high_abundance))
    high_on_low_abundance = np.empty(c.n_species)
    high_on_low_abundance[abundance_order] = rho_values[::-1]
    assignments.append(("high_rho_on_low_baseline_abundance", high_on_low_abundance))

    rng = np.random.default_rng(c.seed + 980001)
    random_assignments = [rng.permutation(rho_values) for _ in range(n_permutations)]
    records, named_species = [], []
    for label, rho in assignments + [(f"random_{i:03d}", values) for i, values in enumerate(random_assignments)]:
        pp = assign_resistance(community, rho)
        control, control_metrics = save_run(c, pp, f"{label}:control", antibiotic_initial=0)
        treated, treated_metrics = save_run(c, pp, f"{label}:treated", antibiotic_initial=1)
        rec = {"assignment": label, "kind": "random_permutation" if label.startswith("random_") else "structured",
               "spearman_rho_vs_growth": spearmanr(rho, community["r"]).statistic,
               "spearman_rho_vs_baseline_abundance": spearmanr(rho, community["x_star"]).statistic,
               "initial_mean_rho": float(np.mean(rho)),
               "control_final_biomass": control_metrics["final_biomass"],
               "treated_final_biomass": treated_metrics["final_biomass"],
               "treated_control_biomass_ratio": treated_metrics["final_biomass"] / max(control_metrics["final_biomass"], 1e-300),
               "control_survivors": control_metrics["survivors"],
               "treated_survivors": treated_metrics["survivors"],
               "control_shannon": control_metrics["final_shannon"],
               "treated_shannon": treated_metrics["final_shannon"],
               "control_weighted_rho": control_metrics["final_weighted_resistance"],
               "treated_weighted_rho": treated_metrics["final_weighted_resistance"],
               "selection_shift_weighted_rho": treated_metrics["final_weighted_resistance"] - control_metrics["final_weighted_resistance"],
               "drug_effect_delta_biomass": treated_metrics["final_biomass"] - control_metrics["final_biomass"]}
        records.append(rec)
        if label in {"seeded_assignment", "high_rho_on_high_growth", "high_rho_on_low_growth",
                     "high_rho_on_high_baseline_abundance", "high_rho_on_low_baseline_abundance", "random_000"}:
            named_species.append(pd.DataFrame({"species": community["species"]["species"],
                                               "intrinsic_growth_r": community["r"],
                                               "baseline_equilibrium_abundance": community["x_star"],
                                               "resistance_rho": rho,
                                               "IC50": pp["IC50"],
                                               "assignment": label,
                                               "control_final_abundance": control["X"][-1],
                                               "treated_final_abundance": treated["X"][-1]}))
    return pd.DataFrame(records), pd.concat(named_species, ignore_index=True)


def make_report(summary, sensitivity, diagnostics, report_path, assignments=None, ensemble=None):
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont

    # ReportLab's built-in Helvetica/Times fonts lack many Greek, subscript, and
    # superscript glyphs used here. Embed Matplotlib's bundled DejaVu Sans so the
    # PDF renders these characters consistently on systems without extra fonts.
    font_dir = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    pdfmetrics.registerFont(TTFont("ReportDejaVu", str(font_dir / "DejaVuSans.ttf")))
    pdfmetrics.registerFont(TTFont("ReportDejaVu-Bold", str(font_dir / "DejaVuSans-Bold.ttf")))
    pdfmetrics.registerFont(TTFont("ReportDejaVu-Oblique", str(font_dir / "DejaVuSans-Oblique.ttf")))
    pdfmetrics.registerFont(TTFont("ReportDejaVu-BoldOblique", str(font_dir / "DejaVuSans-BoldOblique.ttf")))
    pdfmetrics.registerFontFamily("ReportDejaVu", normal="ReportDejaVu", bold="ReportDejaVu-Bold",
                                  italic="ReportDejaVu-Oblique", boldItalic="ReportDejaVu-BoldOblique")

    styles = getSampleStyleSheet()
    for style_name in ["Normal", "BodyText"]:
        styles[style_name].fontName = "ReportDejaVu"
    for style_name in ["Heading1", "Heading2", "Title"]:
        styles[style_name].fontName = "ReportDejaVu-Bold"
    styles.add(ParagraphStyle(name="ReportTableCell", parent=styles["BodyText"],
                              fontName="ReportDejaVu", fontSize=7.5, leading=9))
    styles.add(ParagraphStyle(name="ReportTableHeader", parent=styles["BodyText"],
                              fontName="ReportDejaVu-Bold", fontSize=7.5, leading=9,
                              textColor=colors.white))
    styles.add(ParagraphStyle(name="ReportTitle", parent=styles["Title"], fontName="ReportDejaVu-Bold",
                              alignment=TA_CENTER, fontSize=18, leading=22))
    doc = SimpleDocTemplate(str(report_path), pagesize=letter, rightMargin=0.65*inch,
                            leftMargin=0.65*inch, topMargin=0.55*inch, bottomMargin=0.55*inch)
    story = [Paragraph("Random diagonally stable gLV community with antibiotic exposure", styles["ReportTitle"]),
             Spacer(1, 10), Paragraph("Model equations and parameter meanings", styles["Heading1"]),
             Paragraph("For species i=1,...,n:<br/><br/>"
                       "dX<sub>i</sub>/dt = X<sub>i</sub> [ r<sub>i</sub>(1 - c rho<sub>i</sub>) + "
                       "Σ<sub>j=1</sub><super>n</super> A<sub>ij</sub>X<sub>j</sub> - "
                       "k<sub>max</sub>A(t)<super>h</super>/(IC50<sub>i</sub><super>h</super> + A(t)<super>h</super>) ].<br/><br/>"
                       "X<sub>i</sub> is abundance; r<sub>i</sub> is intrinsic growth; A<sub>ij</sub> is the per-capita effect of species j on i; "
                       "rho<sub>i</sub> in [0,1] is resistance; c is the maximum fractional growth cost at rho=1; and the Hill term is direct antibiotic killing.", styles["BodyText"]),
             Spacer(1, 5), Paragraph("Antibiotic dynamics: dA/dt = -k<sub>deg</sub>A; A(0)=A<sub>0</sub>.<br/>"
                       "Resistance mapping: IC50<sub>i</sub> = IC50<sub>min</sub>(IC50<sub>max</sub>/IC50<sub>min</sub>)<super>rho<sub>i</sub></super>.<br/><br/>"
                       "There is no explicit resource, beta-lactamase, or inhibitor state. The ecological interaction matrix is randomly sampled, "
                       "then constructed as A=W⁻¹(K−P), with W diagonal positive, K skew-symmetric, and P symmetric positive definite. "
                       "Thus WA+A^T W = -2P is negative definite, certifying diagonal stability. A positive target equilibrium x* is sampled and r=-Ax*, "
                       "so the no-drug, zero-cost system has a feasible equilibrium by construction.", styles["BodyText"]),
             Spacer(1, 7), Paragraph("What was run", styles["Heading1"]),
             Paragraph("Matched no-drug controls; antibiotic dose and decay sweeps; species-level resistance selection relative to matched control; "
                       "resistance-cost tradeoff; representative growth-versus-killing trajectories; initial-biomass-by-dose and interaction-strength-by-dose "
                       "factorials; a controlled resistance-to-species assignment experiment with fixed ecology and a fixed rho multiset; a 40-community rank-correlation sensitivity screen; and a multi-seed robustness check. Since gLV has no explicit resource "
                       "pool, initial biomass and interaction strength replace CRM resource availability/leakage as context factors.", styles["BodyText"]),
             Spacer(1, 8), Paragraph("Stability and interpretation checks", styles["Heading1"])]
    for k, v in diagnostics.items():
        story.append(Paragraph(f"{k}: {v:.4g}", styles["BodyText"]))
    story.extend([Spacer(1, 8), Paragraph("Scenario ranges", styles["Heading1"])])
    range_rows = [
        ["Parameter", "Values / range", "Interpretation"],
        ["Community", "n=25; interaction scale 0.25 baseline", "Fixed community for matched sweeps; separate stable matrices for ensemble analyses"],
        ["Drug dose A₀", "0, 0.1, 0.25, 0.5, 1, 2", "Initial antibiotic concentration (normalized units)"],
        ["Drug decay kdeg", "0, 0.005, 0.01, 0.02, 0.05, 0.1", "First-order clearance; 0 means no clearance"],
        ["Resistance cost c", "0, 0.1, 0.2, 0.3, 0.4", "Growth cost at resistance score ρ=1"],
        ["IC₅₀", "0.05–2.0; log-interpolated by ρ", "Species-specific antibiotic susceptibility trait"],
        ["Initial biomass", "0.0025, 0.01, 0.04, 0.16", "Density-dependent gLV context (CRM resource analog)"],
        ["Sensitivity screen", "40 draws; kmax 0.15–0.8; h 1–4", "Exploratory Spearman screen; other ranges listed above"],
    ]
    range_rows = [[Paragraph(str(cell), styles["ReportTableHeader"] if i == 0 else styles["ReportTableCell"])
                   for cell in row] for i, row in enumerate(range_rows)]
    range_table = Table(range_rows, colWidths=[1.05*inch, 1.9*inch, 4.05*inch], repeatRows=1)
    range_table.setStyle(TableStyle([("BACKGROUND", (0,0),(-1,0), colors.HexColor("#243b53")),
                                     ("TEXTCOLOR", (0,0),(-1,0), colors.white),
                                     ("GRID", (0,0),(-1,-1), .3, colors.grey),
                                     ("VALIGN", (0,0),(-1,-1), "TOP"),
                                     ("FONTSIZE", (0,0),(-1,-1), 7.5), ("FONTNAME", (0,0),(-1,-1), "ReportDejaVu"),
                                     ("FONTNAME", (0,0),(-1,0), "ReportDejaVu-Bold"),
                                     ("ROWBACKGROUNDS", (0,1),(-1,-1), [colors.white, colors.HexColor("#eef3f8")])]))
    story.append(range_table)
    key = summary[summary.scenario.isin(["control", "dose=1", "dose=2", "decay=0", "decay=0.1"])][
        ["scenario", "final_biomass", "survivors", "final_shannon", "final_weighted_resistance"]]
    if not key.empty:
        data = [list(key.columns)] + key.round(4).values.tolist()
        tab = Table(data, repeatRows=1)
        tab.setStyle(TableStyle([("BACKGROUND", (0,0),(-1,0), colors.HexColor("#243b53")),
                                 ("TEXTCOLOR", (0,0),(-1,0), colors.white),
                                 ("GRID", (0,0),(-1,-1), .35, colors.grey),
                                 ("FONTSIZE", (0,0),(-1,-1), 8), ("FONTNAME", (0,0),(-1,-1), "ReportDejaVu"),
                                 ("FONTNAME", (0,0),(-1,0), "ReportDejaVu-Bold"),
                                 ("ROWBACKGROUNDS", (0,1),(-1,-1), [colors.white, colors.HexColor("#eef3f8")])]))
        story.extend([Spacer(1, 7), tab])
    dose_rows = summary[summary.scenario.isin(["control", "dose=2"])].drop_duplicates("scenario").set_index("scenario")
    if "control" in dose_rows.index and "dose=2" in dose_rows.index:
        untreated, high = dose_rows.loc["control"], dose_rows.loc["dose=2"]
        biomass_change = 100 * (high.final_biomass / untreated.final_biomass - 1)
        story.extend([Spacer(1, 8), Paragraph("Run-specific findings", styles["Heading1"]),
                      Paragraph(f"In the fixed-seed dose sweep, final biomass changed from {untreated.final_biomass:.3g} without drug to {high.final_biomass:.3g} at A₀=2 ({biomass_change:+.1f}%). "
                                f"Survivors changed from {int(untreated.survivors)} to {int(high.survivors)}, while biomass-weighted resistance changed from {untreated.final_weighted_resistance:.3f} to {high.final_weighted_resistance:.3f}. "
                                "This is consistent with preferential retention of initially more resistant community members; resistance is assigned at initialization and does not evolve here.", styles["BodyText"])])
    decay_rows = summary[summary.scenario.isin(["decay=0", "decay=0.1"])].drop_duplicates("scenario").set_index("scenario")
    if "decay=0" in decay_rows.index and "decay=0.1" in decay_rows.index:
        story.append(Paragraph(f"At A₀=1, the no-decay and fastest-decay cases ended at {decay_rows.loc['decay=0'].final_biomass:.3g} and {decay_rows.loc['decay=0.1'].final_biomass:.3g} total biomass, respectively. This captures exposure-duration effects as well as peak dose.", styles["BodyText"]))
    rankcorr = sensitivity.drop(columns=["sample"], errors="ignore").corr(method="spearman")["final_biomass"].drop("final_biomass").abs().sort_values(ascending=False)
    top_corr = ", ".join(f"{name} (|ρ|={val:.2f})" for name, val in rankcorr.head(3).items())
    story.append(Paragraph(f"The three largest absolute rank correlations with final biomass in the 40-draw sensitivity screen were {top_corr}. They are exploratory associations across joint parameter/community draws, not isolated causal effects.", styles["BodyText"]))
    cost_rows = summary[summary.scenario.astype(str).str.startswith("cost=")]
    biomass_rows = summary[summary.scenario.astype(str).str.startswith("biomass=") & np.isclose(summary.dose.fillna(-1), 1)]
    interaction_rows = summary[summary.scenario.astype(str).str.startswith("interaction=") & np.isclose(summary.dose.fillna(-1), 1)]
    if not cost_rows.empty:
        cbase = cost_rows[np.isclose(cost_rows.resistance_cost, 0) & np.isclose(cost_rows.dose, 0)].iloc[0]
        cmax = cost_rows[np.isclose(cost_rows.resistance_cost, .4) & np.isclose(cost_rows.dose, 0)].iloc[0]
        story.append(Paragraph(f"The untreated resistance-cost sweep changed final biomass from {cbase.final_biomass:.3g} at c=0 to {cmax.final_biomass:.3g} at c=0.4. This is the expected direct fitness-cost effect and is why antibiotic response should be interpreted against a matched control at each cost.", styles["BodyText"]))
    if not biomass_rows.empty and not interaction_rows.empty:
        story.append(Paragraph(f"At A₀=1, initial-biomass scenarios spanned final biomass {biomass_rows.final_biomass.min():.3g}–{biomass_rows.final_biomass.max():.3g}; interaction-strength scenarios spanned {interaction_rows.final_biomass.min():.3g}–{interaction_rows.final_biomass.max():.3g}. The latter use separately sampled stable matrices, so this spread includes the associated change in feasible growth vector, not only interaction magnitude.", styles["BodyText"]))
    if ensemble is not None:
        e1 = ensemble[np.isclose(ensemble.dose, 1)].final_biomass
        if len(e1):
            story.append(Paragraph(f"Across the ten independently sampled communities at A₀=1, final biomass median was {e1.median():.3g} (range {e1.min():.3g}–{e1.max():.3g}); these seeds vary both ecology and resistance assignment.", styles["BodyText"]))
    if assignments is not None:
        random_rows = assignments[assignments.kind == "random_permutation"]
        aligned = assignments[assignments.assignment == "high_rho_on_high_growth"].iloc[0]
        anti = assignments[assignments.assignment == "high_rho_on_low_growth"].iloc[0]
        if not random_rows.empty:
            q05, med, q95 = random_rows.treated_control_biomass_ratio.quantile([.05, .5, .95])
            story.append(Paragraph(f"Resistance-placement experiment: the exact same 25 rho values and gLV community were retained while their species assignment was changed. "
                                   f"Across 200 random assignments, treated/untreated final-biomass ratio had median {med:.3f} (5th–95th percentile {q05:.3f}–{q95:.3f}). "
                                   f"Assigning the highest rho values to highest-growth species gave ratio {aligned.treated_control_biomass_ratio:.3f}; assigning them to lowest-growth species gave {anti.treated_control_biomass_ratio:.3f}. "
                                   "Because resistance cost is attached to rho, each assignment was compared to its own no-antibiotic control; this isolates antibiotic response from assignment-dependent baseline cost effects.", styles["BodyText"]))
    story.extend([Spacer(1, 8), Paragraph("Results and figures", styles["Heading1"]),
                  Paragraph("The results are exploratory model behavior, not fitted or validated biological estimates. "
                            "A is held fixed within matched sweeps. The sensitivity analysis is a Spearman rank screen, not a variance-based global sensitivity analysis. "
                            "Diagonal stability establishes a strong property for the autonomous gLV interaction system; it does not by itself prove stability of the full "
                            "time-varying antibiotic system. All numerical conclusions should be read with the stated parameterization and finite simulation horizon in mind.", styles["BodyText"])])
    for name, caption in [
        ("01_control_dynamics", "No-antibiotic reference and randomized community design."),
        ("01b_growth_parameter_heatmap", "Sampled intrinsic growth-rate parameter for each species."),
        ("01c_interaction_parameter_heatmap", "Sampled gLV interaction matrix; rows are focal species and columns are interacting species."),
        ("02_dose_response", "Dose response: total abundance, diversity, survivors, and resistance enrichment."),
        ("02b_dose_endpoints", "Dose-response endpoint metrics for abundance, survivors, diversity, and resistance."),
        ("03_dose_species_heatmap", "Species abundance ratios relative to matched no-antibiotic control."),
        ("04_decay_response", "Drug persistence changes the duration of killing pressure."),
        ("05_resistance_cost", "Resistance fitness-cost tradeoff under control and treatment."),
        ("05b_growth_penalty", "The explicit resistance-associated reduction in intrinsic growth."),
        ("06_growth_kill_mechanism", "Representative intrinsic growth and direct antibiotic killing terms."),
        ("07_biomass_dose_factorial", "Initial community biomass-by-dose factorial."),
        ("08_interaction_dose_factorial", "Interaction-strength-by-dose factorial across stable communities."),
        ("09_sensitivity", "Rank correlations across randomly sampled stable communities and parameter ranges."),
        ("10_ensemble_robustness", "Variation across independently sampled community matrices."),
        ("11_resistance_assignment_effect", "How rho-to-species assignment relates to antibiotic response while ecology and the rho multiset are held fixed."),
        ("12_resistance_assignment_map", "The exact resistance values assigned to each species in the structured placements and example random assignment."),
        ("13_assignment_outcomes", "Matched treated/control species outcomes under distinct resistance placements."),
    ]:
        path = FIG / f"{name}.png"
        if path.exists():
            story.extend([PageBreak(), Paragraph(caption, styles["Heading2"]), Spacer(1, 5),
                          Image(str(path), width=7.05*inch, height=4.65*inch, kind="proportional")])
    story.append(PageBreak())
    story.extend([Paragraph("How to read the main patterns", styles["Heading1"]),
                  Paragraph("A decrease in total abundance after exposure reflects the balance between growth, density-dependent interactions, and antibiotic killing. "
                            "A higher final biomass-weighted IC₅₀ or resistance score than in the matched control is consistent with ecological selection for resistant species, "
                            "but does not mean resistance arose during the simulation. Resistance is assigned at initialization. A cost can counter selection by reducing intrinsic growth. "
                            "The biomass-by-dose factorial tests whether initial density changes the trajectory through gLV feedback; the interaction-strength factorial tests how the "
                            "sampled community structure modifies treatment outcomes. Ensemble variability shows whether a result is specific to the fixed seed.", styles["BodyText"]),
                  Spacer(1, 10), Paragraph("Supporting CSV files, full per-species outcomes, sampled matrices, model settings, and every plot in PNG/PDF formats are saved alongside this report in outputs/.", styles["BodyText"])])
    doc.build(story)


def write_final_observations(out_dir, dose_df, decay_df, cost_df, biomass_df,
                             interaction_df, sensitivity_df, ensemble_df,
                             assignments_df, diagnostics):
    """Write a data-backed narrative summary of all completed gLV analyses."""
    def metric_at(frame, column, value, by="dose"):
        row = frame[np.isclose(frame[by].astype(float), float(value))].iloc[0]
        return row

    dose0, dose1, dose2 = (metric_at(dose_df, "dose", x) for x in (0, 1, 2))
    decay0, decayfast = (metric_at(decay_df, "decay", x, by="decay") for x in (0, .1))
    c0 = cost_df[np.isclose(cost_df.resistance_cost, 0) & np.isclose(cost_df.dose, 0)].iloc[0]
    c4 = cost_df[np.isclose(cost_df.resistance_cost, .4) & np.isclose(cost_df.dose, 0)].iloc[0]
    t0 = cost_df[np.isclose(cost_df.resistance_cost, 0) & np.isclose(cost_df.dose, 1)].iloc[0]
    t4 = cost_df[np.isclose(cost_df.resistance_cost, .4) & np.isclose(cost_df.dose, 1)].iloc[0]
    b0 = biomass_df[np.isclose(biomass_df.dose, 1)].sort_values("initial_biomass")
    i0 = interaction_df[np.isclose(interaction_df.dose, 1)].sort_values("interaction_strength")
    rho_corr = sensitivity_df.drop(columns=["sample"], errors="ignore").corr(method="spearman")["final_biomass"].drop("final_biomass").sort_values(key=np.abs, ascending=False)
    rand = assignments_df[assignments_df.kind == "random_permutation"]
    top_rand = rand.nlargest(1, "treated_control_biomass_ratio").iloc[0]
    low_rand = rand.nsmallest(1, "treated_control_biomass_ratio").iloc[0]
    align_hi = assignments_df[assignments_df.assignment == "high_rho_on_high_growth"].iloc[0]
    align_lo = assignments_df[assignments_df.assignment == "high_rho_on_low_growth"].iloc[0]
    abundance_hi = assignments_df[assignments_df.assignment == "high_rho_on_high_baseline_abundance"].iloc[0]
    abundance_lo = assignments_df[assignments_df.assignment == "high_rho_on_low_baseline_abundance"].iloc[0]
    ens = ensemble_df.groupby("dose").final_biomass.agg(["median", "min", "max"]).reset_index()

    lines = [
        "# Final observations: antibiotic-treated gLV system", "",
        "## Scope and model interpretation", "",
        "This report summarizes a 25-species gLV community with a randomly sampled, diagonally stable interaction matrix and direct species-specific Hill killing. Antibiotic concentration decays by first-order kinetics. There are no explicit resource, beta-lactamase, inhibitor, or evolving-resistance states. Resistance score rho is assigned at initialization; its fitness cost is included in intrinsic growth.", "",
        "The sampled matrix passed the diagonal-stability certificate: maximum eigenvalue of WA + A^T W = "
        f"{diagnostics['max_eigenvalue_symmetric_part']:.4g} (must be negative), minimum eigenvalue of P = {diagnostics['min_eigenvalue_P']:.4g}. "
        "The target equilibrium residual was " f"{diagnostics['max_equilibrium_residual']:.2g}. This certifies the autonomous gLV interaction matrix, not the full antibiotic-forced dynamics.", "",
        "## Main observations", "",
        f"- **Dose response (one fixed community and resistance assignment):** final biomass was {dose0.final_biomass:.3f} without antibiotic, {dose1.final_biomass:.3f} at A0=1, and {dose2.final_biomass:.3f} at A0=2. At dose 2, {int(dose2.survivors)} of 25 species remained above the numerical extinction threshold versus {int(dose0.survivors)} in control; Shannon diversity shifted from {dose0.final_shannon:.3f} to {dose2.final_shannon:.3f}.",
        f"- **Resistance enrichment in the dose series:** biomass-weighted rho was {dose0.final_weighted_resistance:.3f} in control and {dose2.final_weighted_resistance:.3f} at dose 2. This indicates ecological enrichment of initially more resistant members in this sampled community, not evolution of resistance.",
        f"- **Drug persistence:** with A0=1, final biomass was {decay0.final_biomass:.3f} at kdeg=0 and {decayfast.final_biomass:.3f} at kdeg=0.1. Faster drug decay shortens the killing period; this is a matched comparison with other settings fixed.",
        f"- **Resistance cost:** in the untreated control, raising the maximum cost from 0 to 0.4 changed final biomass {c0.final_biomass:.3f} -> {c4.final_biomass:.3f}. Under A0=1, the corresponding values were {t0.final_biomass:.3f} -> {t4.final_biomass:.3f}. Interpret treated values alongside their same-cost controls because cost acts even without antibiotic.",
        f"- **Initial biomass x dose:** at A0=1, final biomass ranged from {b0.final_biomass.min():.3f} to {b0.final_biomass.max():.3f} across the four starting biomass levels. In gLV this is a density-feedback context test, not a resource-availability test.",
        f"- **Interaction-strength x dose:** at A0=1, the four separately sampled stable communities produced final biomass from {i0.final_biomass.min():.3f} to {i0.final_biomass.max():.3f} over interaction scales 0.1–0.55. Since each scale uses a different sampled matrix and associated feasible growth vector, this is an ensemble comparison, not a one-parameter perturbation of a single A.",
        f"- **Sensitivity screen:** 40 joint parameter/community draws gave the largest absolute Spearman associations with final biomass for " + ", ".join(f"{idx} (rho={val:.2f})" for idx, val in rho_corr.head(3).items()) + ". These are exploratory screening signals, not causal effects or formal variance-based sensitivity indices.",
        f"- **Across-community robustness:** across ten independently sampled communities, final-biomass median (range) by dose was " + "; ".join(f"A0={row.dose:g}: {row['median']:.3f} ({row['min']:.3f}–{row['max']:.3f})" for _, row in ens.iterrows()) + ". These communities also resample resistance assignment, so this spread includes both ecological and trait-assignment differences.",
        "", "## Which species receive higher initial resistance?", "",
        "To isolate placement, this experiment held A, r, x*, the multiset of 25 rho values, dose, decay, and cost fixed while changing only which species received each rho. Each assignment was run both untreated and treated, so antibiotic response is measured relative to that assignment's own baseline and the resistance cost is not mistaken for drug benefit.",
        f"- Pairing the highest rho values with the highest intrinsic-growth species gave a treated/control final-biomass ratio of {align_hi.treated_control_biomass_ratio:.3f}; pairing high rho with the lowest-growth species gave {align_lo.treated_control_biomass_ratio:.3f}.",
        f"- Pairing high rho with the highest versus lowest baseline-equilibrium-abundance species gave treated/control ratios {abundance_hi.treated_control_biomass_ratio:.3f} and {abundance_lo.treated_control_biomass_ratio:.3f}, respectively.",
        f"- Across {len(rand)} random permutations, the treated/control biomass ratio ranged from {rand.treated_control_biomass_ratio.min():.3f} to {rand.treated_control_biomass_ratio.max():.3f}, with median {rand.treated_control_biomass_ratio.median():.3f} and 5th–95th percentiles {rand.treated_control_biomass_ratio.quantile(.05):.3f}–{rand.treated_control_biomass_ratio.quantile(.95):.3f}. The lowest-ratio placement was {low_rand.assignment} ({low_rand.treated_control_biomass_ratio:.3f}); the highest was {top_rand.assignment} ({top_rand.treated_control_biomass_ratio:.3f}).",
        "- The assignment-specific response can differ because species vary in intrinsic growth, equilibrium abundance, and interaction-network position. The structured high/low-growth and high/low-equilibrium-abundance placements ask whether resistance is allocated to different ecological roles; the random-permutation distribution estimates variation across arbitrary allocations.",
        "", "## Limitations and interpretation", "",
        "All parameter values are exploratory normalized choices, not fitted biological estimates. Numerical values of kmax and time scales are not directly comparable between CRM and gLV without scale calibration. Results depend on one primary sampled community except for the explicitly labeled ensemble screens. Resistance is a pre-existing species trait, not within-species adaptation; no mutation, horizontal transfer, or changing resistance is modeled. The random-permutation experiment holds the rho multiset fixed but includes the associated resistance-growth cost, and its treated/control matching is therefore essential."
    ]
    (out_dir / "reports" / "gLV_Final_Observations.md").write_text("\n".join(lines) + "\n")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    c = Config()
    p = sample_community(c)
    diagnostics = stability_diagnostics(p)
    if diagnostics["max_eigenvalue_symmetric_part"] >= 1e-9:
        raise RuntimeError("Diagonal stability certificate check failed")
    (OUT / "config.json").write_text(json.dumps({"config": c.__dict__, "stability_diagnostics": diagnostics}, indent=2))
    pd.DataFrame(p["A"]).to_csv(OUT / "interaction_matrix_A.csv", index=False, header=False)
    pd.DataFrame(p["W"]).to_csv(OUT / "stability_certificate_W.csv", index=False, header=False)
    p["species"].to_csv(OUT / "sampled_species_traits.csv", index=False)
    pd.DataFrame(p["x_star"]).to_csv(OUT / "baseline_equilibrium_xstar.csv", index=False, header=["abundance"])
    rows, species_rows = [], []

    # Matched control and dose response.
    doses = [0, .1, .25, .5, 1, 2]
    dose_runs = {}
    for dose in doses:
        res, row = save_run(c, p, "control" if dose == 0 else f"dose={dose:g}", antibiotic_initial=dose)
        rows.append(row); dose_runs[dose] = res
        species_rows.append(species_outcomes(p, res, row["scenario"]).assign(dose=dose))
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for dose, res in dose_runs.items():
        ax.plot(res["t"], res["X"].sum(axis=1), label=f"A₀={dose:g}")
    ax.set(xlabel="Time", ylabel="Total abundance", title="Matched gLV community: antibiotic dose response"); ax.legend(ncol=2)
    savefig(fig, "02_dose_response")
    dose_summary = pd.DataFrame(rows)
    fig, axs = plt.subplots(2, 2, figsize=(9, 7))
    for ax, col, label in zip(axs.flat, ["final_biomass", "survivors", "final_shannon", "final_weighted_resistance"],
                              ["Final biomass", "Surviving species", "Shannon diversity", "Biomass-weighted resistance"]):
        ax.plot(doses, dose_summary.dose, alpha=0) if False else None
        ax.plot(dose_summary.dose, dose_summary[col], marker="o"); ax.set(xlabel="Initial antibiotic A₀", ylabel=label)
    fig.suptitle("Dose-response endpoints"); fig.tight_layout(); savefig(fig, "02b_dose_endpoints")
    sframe = pd.concat(species_rows, ignore_index=True)
    ctrl = sframe[sframe.dose == 0].set_index("species").final_abundance
    nonzero = [d for d in doses if d > 0]
    ordered = p["species"].sort_values("IC50").species.tolist()
    heat = []
    for d in nonzero:
        v = sframe[sframe.dose == d].set_index("species").final_abundance / ctrl
        heat.append(np.log10(np.maximum(v.reindex(ordered).values, 1e-12)))
    fig, ax = plt.subplots(figsize=(11, 4.8)); im = ax.imshow(heat, aspect="auto", cmap="coolwarm", vmin=-6, vmax=2)
    ax.set(yticks=range(len(nonzero)), yticklabels=[f"A₀={d:g}" for d in nonzero], xticks=range(c.n_species), xticklabels=ordered,
           xlabel="Species ordered by IC₅₀", title="Final abundance / matched control (log₁₀)"); ax.tick_params(axis="x", rotation=90, labelsize=7)
    fig.colorbar(im, ax=ax, label="log₁₀ ratio"); savefig(fig, "03_dose_species_heatmap")
    sframe.to_csv(OUT / "dose_species_outcomes.csv", index=False)

    # Decay sweep at dose=1.
    decay_rows, decay_runs = [], {}
    for decay in [0, .005, .01, .02, .05, .1]:
        res, row = save_run(c, p, f"decay={decay:g}", antibiotic_initial=1, antibiotic_decay=decay)
        decay_rows.append(row); decay_runs[decay] = res
    pd.DataFrame(decay_rows).to_csv(OUT / "decay_summary.csv", index=False)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for d, res in decay_runs.items():
        axs[0].plot(res["t"], res["X"].sum(axis=1), label=f"kdeg={d:g}"); axs[1].plot(res["t"], res["drug"], label=f"kdeg={d:g}")
    axs[0].set(xlabel="Time", ylabel="Total abundance", title="Community trajectories"); axs[1].set(xlabel="Time", ylabel="Antibiotic", title="Exposure trajectories")
    axs[0].legend(fontsize=7); axs[1].legend(fontsize=7); fig.tight_layout(); savefig(fig, "04_decay_response")

    # Resistance-cost matched control and treated comparisons.
    costs = [0, .1, .2, .3, .4]
    cost_rows = []
    for cost in costs:
        for dose in [0, 1]:
            res, row = save_run(c, p, f"cost={cost:g},dose={dose:g}", antibiotic_initial=dose, resistance_cost_max=cost)
            cost_rows.append(row)
    cost_df = pd.DataFrame(cost_rows); cost_df.to_csv(OUT / "resistance_cost_summary.csv", index=False)
    fig, axs = plt.subplots(1, 2, figsize=(9, 4))
    for dose, group in cost_df.groupby("dose"):
        axs[0].plot(group.resistance_cost, group.final_biomass, marker="o", label=f"A₀={dose:g}")
        axs[1].plot(group.resistance_cost, group.final_weighted_resistance, marker="o", label=f"A₀={dose:g}")
    axs[0].set(xlabel="Maximum resistance cost", ylabel="Final biomass"); axs[1].set(xlabel="Maximum resistance cost", ylabel="Weighted resistance")
    for ax in axs: ax.legend()
    fig.suptitle("Resistance fitness-cost tradeoff"); fig.tight_layout(); savefig(fig, "05_resistance_cost")
    rho = p["resistance"]
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for cost in [0, .2, .4]:
        ax.scatter(rho, p["r"] * (1-cost*rho), s=20, label=f"cost={cost:g}")
    ax.set(xlabel="Resistance score ρ", ylabel="Effective intrinsic growth rᵢ(1−cρᵢ)", title="Mechanistic growth penalty"); ax.legend()
    savefig(fig, "05b_growth_penalty")

    # Resistance selection relative to matched controls.
    ctrl_species = species_outcomes(p, dose_runs[0], "control").set_index("species")
    treated_species = species_outcomes(p, dose_runs[1], "dose=0.1").set_index("species")
    treated_species["treatment_control_log10_ratio"] = np.log10(np.maximum(treated_species.final_abundance, 1e-300) / np.maximum(ctrl_species.final_abundance, 1e-300))
    treated_species.reset_index().to_csv(OUT / "resistance_selection_species.csv", index=False)
    corr = spearmanr(p["resistance"], treated_species.treatment_control_log10_ratio)
    pd.DataFrame([{"dose": 0.1, "spearman_rho": corr.statistic, "p_value": corr.pvalue,
                   "control_weighted_resistance": rows[0]["final_weighted_resistance"],
                   "treated_weighted_resistance": rows[1]["final_weighted_resistance"]}]).to_csv(OUT / "resistance_selection_summary.csv", index=False)

    # Growth-vs-killing mechanism for low/median/high resistance representatives.
    reps = np.argsort(p["resistance"])[[0, c.n_species//2, -1]]
    mech_rows = []
    res = dose_runs[1]
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    for i in reps:
        kill = p["IC50"][i]
        kt = c.antibiotic_kmax * res["drug"]**c.antibiotic_hill / (kill**c.antibiotic_hill + res["drug"]**c.antibiotic_hill)
        axs[0].plot(res["t"], res["X"][:, i], label=f"S{i+1:02d} ρ={rho[i]:.2f}")
        axs[1].plot(res["t"], res["r_eff"][i] - kt, label=f"S{i+1:02d} net low-density rate")
        mech_rows.append({"species": f"S{i+1:02d}", "resistance": rho[i], "IC50": p["IC50"][i], "intrinsic_growth": p["r"][i]})
    axs[0].set(xlabel="Time", ylabel="Abundance", title="Representative species"); axs[1].axhline(0, color="black", lw=.8); axs[1].set(xlabel="Time", ylabel="rᵢ,eff − antibiotic kill", title="Low-density net per-capita rate")
    axs[0].legend(fontsize=7); axs[1].legend(fontsize=7); fig.tight_layout(); savefig(fig, "06_growth_kill_mechanism")
    pd.DataFrame(mech_rows).to_csv(OUT / "mechanism_representatives.csv", index=False)

    # Initial biomass-by-dose factorial: resource-level analog for density-dependent gLV.
    factorial = []
    for biomass in [.0025, .01, .04, .16]:
        for dose in [0, .25, .5, 1, 2]:
            rr, row = save_run(c, p, f"biomass={biomass:g},dose={dose:g}", initial_biomass=biomass, antibiotic_initial=dose)
            factorial.append(row)
    fac_df = pd.DataFrame(factorial); fac_df.to_csv(OUT / "biomass_dose_factorial.csv", index=False)
    pivot = fac_df.pivot(index="initial_biomass", columns="dose", values="final_biomass")
    fig, ax = plt.subplots(figsize=(7, 4.8)); im = ax.imshow(np.log10(np.maximum(pivot.values, 1e-15)), aspect="auto", cmap="viridis")
    ax.set(xticks=range(len(pivot.columns)), xticklabels=[f"{x:g}" for x in pivot.columns], yticks=range(len(pivot.index)), yticklabels=[f"{x:g}" for x in pivot.index], xlabel="Initial antibiotic A₀", ylabel="Initial total biomass", title="Final biomass, log₁₀ scale")
    fig.colorbar(im, ax=ax, label="log₁₀ final biomass"); savefig(fig, "07_biomass_dose_factorial")

    # Interaction-strength x dose factorial; every A is stable by construction.
    int_rows = []
    for strength in [.1, .25, .4, .55]:
        cc = Config(**{**c.__dict__, "interaction_strength": strength, "seed": 100 + int(strength*100)})
        pp = sample_community(cc)
        for dose in [0, .25, .5, 1, 2]:
            rr, row = save_run(cc, pp, f"interaction={strength:g},dose={dose:g}", antibiotic_initial=dose)
            row["interaction_strength"] = strength; int_rows.append(row)
    int_df = pd.DataFrame(int_rows); int_df.to_csv(OUT / "interaction_dose_factorial.csv", index=False)
    piv = int_df.pivot(index="interaction_strength", columns="dose", values="final_biomass")
    fig, ax = plt.subplots(figsize=(7, 4.8)); im = ax.imshow(np.log10(np.maximum(piv.values, 1e-15)), aspect="auto", cmap="magma")
    ax.set(xticks=range(len(piv.columns)), xticklabels=[f"{x:g}" for x in piv.columns], yticks=range(len(piv.index)), yticklabels=[f"{x:g}" for x in piv.index], xlabel="Initial antibiotic A₀", ylabel="Interaction scale", title="Final biomass across stable sampled communities")
    fig.colorbar(im, ax=ax, label="log₁₀ final biomass"); savefig(fig, "08_interaction_dose_factorial")

    # Global rank screen: 40 independently sampled, feasible stable communities.
    sens = []
    srng = np.random.default_rng(7301)
    for j in range(40):
        strength = float(srng.uniform(.08, .55)); dose = float(srng.uniform(0, 2))
        decay = float(srng.uniform(0, .1)); kmax = float(srng.uniform(.15, .8))
        hill = float(srng.uniform(1, 4)); cost = float(srng.uniform(0, .4))
        biomass = float(np.exp(srng.uniform(np.log(.0025), np.log(.16))))
        cc = Config(**{**c.__dict__, "seed": 1000+j, "interaction_strength": strength,
                       "antibiotic_initial": dose, "antibiotic_decay": decay,
                       "antibiotic_kmax": kmax, "antibiotic_hill": hill,
                       "resistance_cost_max": cost, "initial_biomass": biomass})
        pp = sample_community(cc)
        rr, mm = save_run(cc, pp, f"sample={j}")
        sens.append({"sample": j, "interaction_strength": strength, "dose": dose, "decay": decay,
                     "kmax": kmax, "hill": hill, "resistance_cost": cost,
                     "initial_biomass": biomass, **{k: mm[k] for k in ["final_biomass", "survivors", "final_shannon", "final_weighted_resistance"]}})
    sensdf = pd.DataFrame(sens); sensdf.to_csv(OUT / "sensitivity_runs.csv", index=False)
    inputs = ["interaction_strength", "dose", "decay", "kmax", "hill", "resistance_cost", "initial_biomass"]
    outputs = ["final_biomass", "survivors", "final_shannon", "final_weighted_resistance"]
    corrmat = pd.DataFrame(index=inputs, columns=outputs, dtype=float)
    for x in inputs:
        for y in outputs: corrmat.loc[x,y] = spearmanr(sensdf[x], sensdf[y]).statistic
    corrmat.to_csv(OUT / "sensitivity_spearman_rho.csv")
    fig, ax = plt.subplots(figsize=(8, 5)); im = ax.imshow(corrmat.values, vmin=-1, vmax=1, cmap="coolwarm", aspect="auto")
    ax.set(xticks=range(len(outputs)), xticklabels=outputs, yticks=range(len(inputs)), yticklabels=inputs, title="Spearman rank correlations (n=40 sampled communities)")
    for i in range(len(inputs)):
        for j in range(len(outputs)): ax.text(j, i, f"{corrmat.iloc[i,j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="Spearman ρ"); savefig(fig, "09_sensitivity")

    # A dose series over ten additional seeds quantifies community-to-community spread.
    ensemble = []
    for seed in range(10, 20):
        cc = Config(**{**c.__dict__, "seed": seed})
        pp = sample_community(cc)
        for dose in [0, .5, 1, 2]:
            rr, mm = save_run(cc, pp, f"seed={seed},dose={dose:g}", antibiotic_initial=dose)
            mm["seed"] = seed; ensemble.append(mm)
    ens = pd.DataFrame(ensemble); ens.to_csv(OUT / "ensemble_robustness.csv", index=False)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for dose, group in ens.groupby("dose"):
        ax.scatter(np.full(len(group), dose), group.final_biomass, alpha=.7)
        ax.plot([dose-.12, dose+.12], [group.final_biomass.median()]*2, lw=2)
    ax.set(xlabel="Initial antibiotic A₀", ylabel="Final biomass", title="Across 10 independently sampled stable communities"); savefig(fig, "10_ensemble_robustness")

    # Baseline/control design diagnostic panel.
    control = dose_runs[0]
    fig, axs = plt.subplots(1, 2, figsize=(9, 4))
    axs[0].plot(control["t"], control["X"].sum(axis=1), color="#2878b5"); axs[0].set(xlabel="Time", ylabel="Total abundance", title="No-antibiotic reference")
    axs[1].scatter(p["resistance"], p["IC50"], c=p["r"], cmap="viridis"); axs[1].set(xlabel="Resistance score ρ", ylabel="IC₅₀", title="Sampled species traits")
    fig.tight_layout(); savefig(fig, "01_control_dynamics")

    # Show the sampled gLV parameters directly, not just their dynamical outcomes.
    fig, ax = plt.subplots(figsize=(4.8, 7.2))
    growth_image = ax.imshow(p["r"][:, None], aspect="auto", cmap="viridis")
    ax.set(xticks=[0], xticklabels=["rᵢ"], yticks=np.arange(c.n_species),
           yticklabels=p["species"]["species"], ylabel="Species",
           title="Sampled intrinsic growth rates")
    fig.colorbar(growth_image, ax=ax, label="Intrinsic growth rate rᵢ")
    savefig(fig, "01b_growth_parameter_heatmap")

    fig, ax = plt.subplots(figsize=(9, 8))
    bound = float(np.max(np.abs(p["A"])))
    interaction_image = ax.imshow(p["A"], cmap="coolwarm", vmin=-bound, vmax=bound, aspect="equal")
    ax.set(xticks=np.arange(c.n_species), xticklabels=p["species"]["species"],
           yticks=np.arange(c.n_species), yticklabels=p["species"]["species"],
           xlabel="Interacting species j", ylabel="Focal species i",
           title="Sampled gLV interaction matrix A")
    ax.tick_params(axis="x", rotation=90, labelsize=7); ax.tick_params(axis="y", labelsize=7)
    fig.colorbar(interaction_image, ax=ax, label="Aᵢⱼ (negative = inhibition; positive = facilitation)")
    savefig(fig, "01c_interaction_parameter_heatmap")

    # Controlled rho-placement study: identical ecology and identical rho multiset.
    assignment_df, assignment_species = resistance_assignment_sweep(c, p, n_permutations=200)
    assignment_df.to_csv(OUT / "resistance_assignment_permutations.csv", index=False)
    random_assignments = assignment_df[assignment_df.kind == "random_permutation"]
    structured_assignments = assignment_df[assignment_df.kind == "structured"]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    sc = ax.scatter(random_assignments.spearman_rho_vs_growth,
                    random_assignments.treated_control_biomass_ratio,
                    c=random_assignments.selection_shift_weighted_rho,
                    cmap="viridis", alpha=.65, edgecolor="none", label="Random rho permutations")
    for _, row in structured_assignments.iterrows():
        ax.scatter(row.spearman_rho_vs_growth, row.treated_control_biomass_ratio,
                   marker="*", s=130, edgecolor="black", zorder=3)
        short = {
            "seeded_assignment": "seeded",
            "high_rho_on_high_growth": "ρ→high growth",
            "high_rho_on_low_growth": "ρ→low growth",
            "high_rho_on_high_baseline_abundance": "ρ→high abundance",
            "high_rho_on_low_baseline_abundance": "ρ→low abundance",
        }[row.assignment]
        align = "right" if row.spearman_rho_vs_growth > .75 else "left"
        offset = (-6, 5) if align == "right" else (5, 5)
        ax.annotate(short, (row.spearman_rho_vs_growth, row.treated_control_biomass_ratio),
                    xytext=offset, textcoords="offset points", fontsize=8, ha=align)
    ax.axvline(0, color="grey", lw=.8); ax.set(xlabel="Spearman correlation: assigned ρ vs intrinsic growth",
        ylabel="Treated / matched-control final biomass", title="Effect of which species receive high initial resistance")
    ax.set_xlim(-1.2, 1.2)
    fig.colorbar(sc, ax=ax, label="Treatment shift in biomass-weighted ρ")
    savefig(fig, "11_resistance_assignment_effect")

    placement_order = ["seeded_assignment", "high_rho_on_high_growth", "high_rho_on_low_growth",
                       "high_rho_on_high_baseline_abundance", "high_rho_on_low_baseline_abundance", "random_000"]
    species_order = p["species"].sort_values("intrinsic_growth_r")["species"].tolist()
    rho_map = assignment_species.pivot(index="assignment", columns="species", values="resistance_rho").reindex(index=placement_order, columns=species_order)
    fig, ax = plt.subplots(figsize=(12, 4.6))
    im = ax.imshow(rho_map.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set(yticks=np.arange(len(rho_map.index)), yticklabels=[x.replace("high_rho_on_", "high ρ on ").replace("high_baseline_abundance", "high abundance").replace("low_baseline_abundance", "low abundance").replace("growth", "growth") for x in rho_map.index],
           xticks=np.arange(c.n_species), xticklabels=species_order, xlabel="Species ordered by intrinsic growth rate",
           title="Fixed resistance values reassigned across species")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    fig.colorbar(im, ax=ax, label="Assigned resistance score ρ")
    savefig(fig, "12_resistance_assignment_map")

    assignment_species["treatment_control_log10_ratio"] = np.log10(
        np.maximum(assignment_species.treated_final_abundance, 1e-300) /
        np.maximum(assignment_species.control_final_abundance, 1e-300))
    assignment_species.to_csv(OUT / "resistance_assignment_species_outcomes.csv", index=False)
    outcome_map = assignment_species.pivot(index="assignment", columns="species", values="treatment_control_log10_ratio").reindex(index=placement_order, columns=species_order)
    fig, ax = plt.subplots(figsize=(12, 4.6))
    im = ax.imshow(outcome_map.values, aspect="auto", cmap="coolwarm", vmin=-6, vmax=2)
    ax.set(yticks=np.arange(len(outcome_map.index)), yticklabels=[x.replace("high_rho_on_", "high ρ on ").replace("high_baseline_abundance", "high abundance").replace("low_baseline_abundance", "low abundance") for x in outcome_map.index],
           xticks=np.arange(c.n_species), xticklabels=species_order, xlabel="Species ordered by intrinsic growth rate",
           title="Species-level antibiotic response by resistance placement")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    fig.colorbar(im, ax=ax, label="log₁₀(treated abundance / matched control)")
    savefig(fig, "13_assignment_outcomes")

    all_summary = pd.concat([pd.DataFrame(rows), pd.DataFrame(decay_rows), cost_df,
                             fac_df, int_df], ignore_index=True, sort=False)
    all_summary.to_csv(OUT / "scenario_summary.csv", index=False)
    (OUT / "stability_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    pd.DataFrame(rows).to_csv(OUT / "dose_summary.csv", index=False)
    write_final_observations(OUT, pd.DataFrame(rows), pd.DataFrame(decay_rows), cost_df,
                             fac_df, int_df, sensdf, ens, assignment_df, diagnostics)
    make_report(all_summary, sensdf, diagnostics, REPORTS / "gLV_Antibiotic_Analysis_Report.pdf", assignment_df, ens)
    print(f"Completed analysis. Outputs: {OUT}")
    print(f"Report: {REPORTS / 'gLV_Antibiotic_Analysis_Report.pdf'}")


if __name__ == "__main__":
    main()
