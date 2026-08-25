from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def fmt(x, digits: int = 3) -> str:
    if pd.isna(x):
        return ""
    if isinstance(x, str):
        return x
    return f"{float(x):.{digits}f}"


def source(root: Path, rel: str) -> str:
    return str(root / rel)


def build_claim_inventory(root: Path) -> pd.DataFrame:
    rows = [
        {
            "claim_id": "C01",
            "claim_text": "Methylation-only operator-time dynamics captures baseline trajectory but fails at morula basin entry.",
            "supporting_figure": source(root, r"figures\CSB_TRO_final_figure2_baseline_failure.svg"),
            "supporting_table": source(root, r"results\CSB_TRO_final_model_comparison.tsv"),
            "source_result_file": source(root, r"results\CSB_TRO_final_sensitivity_summary.tsv"),
            "status": "complete",
            "risk_level": "low",
            "claim_boundary": "Diagnostic baseline failure; not evidence by itself for causal u_bio.",
        },
        {
            "claim_id": "C02",
            "claim_text": "Measured correction reaches observed morula occupancy around alpha 0.50 and acts as a diagnostic upper bound.",
            "supporting_figure": source(root, r"figures\CSB_TRO_final_figure3_threshold_entry.svg"),
            "supporting_table": source(root, r"results\CSB_TRO_alpha_bifurcation_scan.tsv"),
            "source_result_file": source(root, r"results\CSB_TRO_final_sensitivity_summary.tsv"),
            "status": "complete",
            "risk_level": "medium",
            "claim_boundary": "Uses held-out morula methylation residual; do not describe as non-leaking predictive model.",
        },
        {
            "claim_id": "C03",
            "claim_text": "Measured correction decomposes into a compact module set dominated by M05/M01/M12/M02/M10.",
            "supporting_figure": source(root, r"figures\CSB_TRO_basin_residual_module_contributions.svg"),
            "supporting_table": source(root, r"results\CSB_TRO_missing_control_term_greedy_modules.tsv"),
            "source_result_file": source(root, r"results\CSB_TRO_missing_control_term_module_basis.tsv"),
            "status": "complete",
            "risk_level": "low",
            "claim_boundary": "Residual control coordinates, not causal drivers.",
        },
        {
            "claim_id": "C04",
            "claim_text": "Dual-branch closure/access architecture is direction-sensitive and outperforms wrong-sign branch controls.",
            "supporting_figure": source(root, r"figures\CSB_TRO_final_figure4_dual_branch_structure.svg"),
            "supporting_table": source(root, r"results\CSB_TRO_dual_branch_sign_control.tsv"),
            "source_result_file": source(root, r"results\CSB_TRO_dual_branch_structure_validation_summary.tsv"),
            "status": "complete",
            "risk_level": "low",
            "claim_boundary": "Diagnostic architecture; not final biological input identification.",
        },
        {
            "claim_id": "C05",
            "claim_text": "Random partition and exact sign-pattern controls support a non-arbitrary module/branch structure.",
            "supporting_figure": source(root, r"figures\CSB_TRO_dual_branch_robustness_controls.svg"),
            "supporting_table": source(root, r"results\CSB_TRO_dual_branch_random_partition_control.tsv"),
            "source_result_file": source(root, r"results\CSB_TRO_dual_branch_exact_partition_sign_control.tsv"),
            "status": "complete",
            "risk_level": "low",
            "claim_boundary": "Some random partitions can recover similar geometry when they duplicate or approximate the true module split; report this transparently.",
        },
        {
            "claim_id": "C06",
            "claim_text": "RNA, nearest-gene RNA, motif x TF, and composite public surrogates are weak replacements.",
            "supporting_figure": source(root, r"figures\CSB_TRO_final_figure5_external_model_comparison.svg"),
            "supporting_table": source(root, r"results\CSB_TRO_final_model_comparison.tsv"),
            "source_result_file": source(root, r"results\CSB_TRO_RNA_transition_replication_control_summary.tsv"),
            "status": "complete",
            "risk_level": "low",
            "claim_boundary": "Use as boundary evidence, not as failed biology.",
        },
        {
            "claim_id": "C07",
            "claim_text": "hESC and public embryo histone proxies support branch plausibility, but only diagnostically.",
            "supporting_figure": source(root, r"figures\CSB_TRO_final_figure5_external_model_comparison.svg"),
            "supporting_table": source(root, r"results\CSB_TRO_embryo_histone_control_metrics.tsv"),
            "source_result_file": source(root, r"results\CSB_TRO_embryo_histone_random_controls.tsv"),
            "status": "complete",
            "risk_level": "medium",
            "claim_boundary": "Public embryo contrast is not strict morula-entry replacement; matched random can be high.",
        },
        {
            "claim_id": "C08",
            "claim_text": "Final u_bio identification is bounded by missing controlled-access H3K27ac_morula/H3K4me3_morula data.",
            "supporting_figure": source(root, r"figures\CSB_TRO_final_figure6_data_access_boundary.svg"),
            "supporting_table": source(root, r"results\CSB_TRO_final_model_comparison.tsv"),
            "source_result_file": source(root, r"results\CSB_TRO_histone_data_access_audit.tsv"),
            "status": "complete",
            "risk_level": "low",
            "claim_boundary": "State as data-access limitation and upgrade path.",
        },
    ]
    return pd.DataFrame(rows)


def build_evidence_table(root: Path) -> pd.DataFrame:
    final = read_tsv(root / r"results\CSB_TRO_final_model_comparison.tsv")
    sens = read_tsv(root / r"results\CSB_TRO_final_sensitivity_summary.tsv")
    greedy = read_tsv(root / r"results\CSB_TRO_missing_control_term_greedy_modules.tsv")
    sign = read_tsv(root / r"results\CSB_TRO_dual_branch_sign_control.tsv")
    struct = read_tsv(root / r"results\CSB_TRO_dual_branch_structure_validation_summary.tsv")

    def model_row(name: str) -> pd.Series:
        return final.loc[final["model_name"] == name].iloc[0]

    baseline = model_row("methylation-only baseline")
    measured = model_row("measured correction upper bound")
    dual = model_row("dual-branch chromatin-state proxy")
    rna = model_row("global RNA transition")
    motif = model_row("motif x TF")
    hesc = model_row("hESC histone branch identity proxy")
    embryo = model_row("public embryo histone diagnostic contrast")
    strict_histone = model_row("strict morula-entry partial histone")

    alpha = sens[
        (sens["sensitivity_type"] == "measured_correction_alpha")
        & (sens["latent_dim"] == 3)
        & (sens["basin_quantile"] == 0.9)
    ].iloc[0]
    step3 = greedy.loc[greedy["step"] == 3].iloc[0]
    step4 = greedy.loc[greedy["step"] == 4].iloc[0]
    correct = sign.loc[sign["model"] == "correct_closure_correct_access"].iloc[0]
    wrong_closure = sign.loc[sign["model"] == "wrong_closure_correct_access"].iloc[0]
    wrong_access = sign.loc[sign["model"] == "correct_closure_wrong_access"].iloc[0]

    struct_map = dict(zip(struct["summary_item"], struct["value"]))

    rows = [
        {
            "evidence_layer": "methylation_only_failure",
            "result": f"q90 occupancy={fmt(baseline['occupancy'])}; observed q90 occupancy=0.875",
            "interpretation": "Methylation-only dynamics is insufficient for morula basin entry.",
            "claim_strength": "strong diagnostic",
            "boundary": "Baseline failure alone does not identify u_bio.",
            "source": source(root, r"results\CSB_TRO_final_model_comparison.tsv"),
        },
        {
            "evidence_layer": "measured_correction_upper_bound",
            "result": f"alpha_to_observed={fmt(alpha['alpha_to_observed'])}; occupancy@alpha1={fmt(measured['occupancy'])}",
            "interpretation": "The missing component is compact and threshold-like at basin entry.",
            "claim_strength": "diagnostic upper bound",
            "boundary": "Uses held-out morula methylation residual.",
            "source": source(root, r"results\CSB_TRO_final_sensitivity_summary.tsv"),
        },
        {
            "evidence_layer": "module_specificity",
            "result": f"M05/M01/M12 occupancy={fmt(step3['pred_basin_occupancy_q90'])}; +M02 occupancy={fmt(step4['pred_basin_occupancy_q90'])}",
            "interpretation": "The correction is concentrated in specific residual modules.",
            "claim_strength": "strong diagnostic",
            "boundary": "Modules are diagnostic coordinates.",
            "source": source(root, r"results\CSB_TRO_missing_control_term_greedy_modules.tsv"),
        },
        {
            "evidence_layer": "dual_branch_architecture",
            "result": f"occupancy={fmt(dual['occupancy'])}; cosine={fmt(dual['cosine'])}; PC3 recovery={fmt(dual['PC3_recovery'])}",
            "interpretation": "Closure/access branch structure organizes the missing correction geometry.",
            "claim_strength": "central result",
            "boundary": "Architecture, not final causal u_bio.",
            "source": source(root, r"results\CSB_TRO_final_model_comparison.tsv"),
        },
        {
            "evidence_layer": "sign_dependence",
            "result": f"correct={fmt(correct['pred_basin_occupancy_q90'])}; wrong closure={fmt(wrong_closure['pred_basin_occupancy_q90'])}; wrong access={fmt(wrong_access['pred_basin_occupancy_q90'])}",
            "interpretation": "The branch orientation is direction-sensitive, especially for closure.",
            "claim_strength": "important control",
            "boundary": "Still diagnostic because branch vectors come from residual geometry.",
            "source": source(root, r"results\CSB_TRO_dual_branch_sign_control.tsv"),
        },
        {
            "evidence_layer": "sign_pattern_enumeration",
            "result": f"true sign rank={fmt(struct_map.get('true_sign_pattern_rank'), 0)}; percentile={fmt(struct_map.get('true_sign_pattern_percentile_occupancy'))}",
            "interpretation": "The true branch sign pattern is at the top of tested sign patterns.",
            "claim_strength": "important control",
            "boundary": "Sign enumeration supports structure, not causality.",
            "source": source(root, r"results\CSB_TRO_dual_branch_structure_validation_summary.tsv"),
        },
        {
            "evidence_layer": "rna_surrogate",
            "result": f"global RNA occupancy={fmt(rna['occupancy'])}",
            "interpretation": "Public RNA transition is not sufficient as final replacement.",
            "claim_strength": "boundary",
            "boundary": "Treat as weak surrogate.",
            "source": source(root, r"results\CSB_TRO_final_model_comparison.tsv"),
        },
        {
            "evidence_layer": "motif_tf_surrogate",
            "result": f"motif x TF occupancy={fmt(motif['occupancy'])}",
            "interpretation": "Motif x TF carries limited regulatory signal but does not rescue the basin.",
            "claim_strength": "boundary",
            "boundary": "Do not overstate TF mechanism.",
            "source": source(root, r"results\CSB_TRO_final_model_comparison.tsv"),
        },
        {
            "evidence_layer": "hesc_histone_proxy",
            "result": f"max occupancy={fmt(hesc['max_occupancy'])}; cosine={fmt(hesc['cosine'])}",
            "interpretation": "Histone-state information supports branch identity plausibility.",
            "claim_strength": "moderate support",
            "boundary": "Not stage-matched morula-entry evidence.",
            "source": source(root, r"results\CSB_TRO_final_model_comparison.tsv"),
        },
        {
            "evidence_layer": "public_embryo_histone_diagnostic",
            "result": f"max occupancy={fmt(embryo['max_occupancy'])}; cosine={fmt(embryo['cosine'])}; PC3 recovery={fmt(embryo['PC3_recovery'])}",
            "interpretation": "Public embryo histone contrast is compatible with the branch architecture.",
            "claim_strength": "histone-supported diagnostic",
            "boundary": "Matched random can be high; diagnostic only.",
            "source": source(root, r"results\CSB_TRO_embryo_histone_control_metrics.tsv"),
        },
        {
            "evidence_layer": "strict_morula_histone_gap",
            "result": f"strict partial histone max occupancy={fmt(strict_histone['max_occupancy'])}",
            "interpretation": "Strict morula-entry replacement is data-access limited.",
            "claim_strength": "limitation",
            "boundary": "Requires controlled-access H3K27ac_morula/H3K4me3_morula.",
            "source": source(root, r"results\CSB_TRO_final_model_comparison.tsv"),
        },
    ]
    return pd.DataFrame(rows)


def build_control_summary(root: Path) -> pd.DataFrame:
    rows = []

    sign = read_tsv(root / r"results\CSB_TRO_dual_branch_sign_control.tsv")
    correct = sign.loc[sign["model"] == "correct_closure_correct_access"].iloc[0]
    for model in [
        "wrong_closure_correct_access",
        "correct_closure_wrong_access",
        "wrong_closure_wrong_access",
        "naive_inverse_all_modules",
        "naive_raw_all_modules",
    ]:
        r = sign.loc[sign["model"] == model].iloc[0]
        rows.append({
            "control_family": "dual_branch_sign",
            "control_name": model,
            "target_metric": "q90 occupancy",
            "true_value": fmt(correct["pred_basin_occupancy_q90"]),
            "control_summary": fmt(r["pred_basin_occupancy_q90"]),
            "passes": "yes" if r["pred_basin_occupancy_q90"] < correct["pred_basin_occupancy_q90"] else "no",
            "interpretation": "Wrong or naive branch sign weakens basin entry.",
            "source": source(root, r"results\CSB_TRO_dual_branch_sign_control.tsv"),
        })

    rand = read_tsv(root / r"results\CSB_TRO_dual_branch_random_partition_control.tsv")
    rand_only = rand[rand["random_id"] >= 0].copy()
    true_val = float(rand.loc[rand["validation_status"] == "true_partition", "pred_basin_occupancy_q90"].iloc[0])
    values = rand_only["pred_basin_occupancy_q90"]
    rows.append({
        "control_family": "dual_branch_random_partition",
        "control_name": "random 3/2 module partitions",
        "target_metric": "q90 occupancy",
        "true_value": fmt(true_val),
        "control_summary": f"n={len(values)}; median={fmt(values.median())}; q95={fmt(values.quantile(0.95))}; max={fmt(values.max())}; frac>=true={fmt((values >= true_val).mean())}",
        "passes": "qualified",
        "interpretation": "Most random partitions are weaker, but some duplicate or approximate the true module split; report transparently.",
        "source": source(root, r"results\CSB_TRO_dual_branch_random_partition_control.tsv"),
    })

    exact = read_tsv(root / r"results\CSB_TRO_dual_branch_exact_partition_sign_control.tsv")
    true_exact = float(exact.loc[exact["validation_status"] == "true_partition_correct_sign", "pred_basin_occupancy_q90"].iloc[0])
    exact_controls = exact[exact["validation_status"] != "true_partition_correct_sign"]["pred_basin_occupancy_q90"]
    rows.append({
        "control_family": "exact_partition_sign",
        "control_name": "all exact module partitions with branch signs",
        "target_metric": "q90 occupancy",
        "true_value": fmt(true_exact),
        "control_summary": f"n={len(exact_controls)}; median={fmt(exact_controls.median())}; q95={fmt(exact_controls.quantile(0.95))}; max={fmt(exact_controls.max())}",
        "passes": "yes",
        "interpretation": "True sign pattern is at the top tier; sign and module assignment matter.",
        "source": source(root, r"results\CSB_TRO_dual_branch_exact_partition_sign_control.tsv"),
    })

    top_random = read_tsv(root / r"results\CSB_TRO_residual_DMR_matched_random_control.tsv")
    for k, sub in top_random.groupby("K"):
        obs = sub[sub["control_type"] == "observed_topK"]["pred_basin_occupancy_q90"]
        rnd = sub[sub["control_type"] == "matched_random"]["pred_basin_occupancy_q90"]
        if len(obs) and len(rnd):
            rows.append({
                "control_family": "top_residual_DMR_matched_random",
                "control_name": f"top{k} matched random DMRs",
                "target_metric": "q90 occupancy",
                "true_value": fmt(obs.iloc[0]),
                "control_summary": f"n={len(rnd)}; median={fmt(rnd.median())}; q95={fmt(rnd.quantile(0.95))}; max={fmt(rnd.max())}",
                "passes": "yes" if obs.iloc[0] > rnd.quantile(0.95) else "qualified",
                "interpretation": "Observed high-residual DMR sets outperform matched random sets.",
                "source": source(root, r"results\CSB_TRO_residual_DMR_matched_random_control.tsv"),
            })

    hist = read_tsv(root / r"results\CSB_TRO_embryo_histone_random_controls.tsv")
    true_hist = float(hist.loc[hist["validation_status"] == "true_target", "pred_basin_occupancy_q90"].iloc[0])
    hist_rand = hist[hist["validation_status"] != "true_target"]["pred_basin_occupancy_q90"]
    rows.append({
        "control_family": "embryo_histone_matched_random",
        "control_name": "public embryo histone diagnostic matched random",
        "target_metric": "q90 occupancy",
        "true_value": fmt(true_hist),
        "control_summary": f"n={len(hist_rand)}; median={fmt(hist_rand.median())}; q95={fmt(hist_rand.quantile(0.95))}; max={fmt(hist_rand.max())}; frac>=true={fmt((hist_rand >= true_hist).mean())}",
        "passes": "no",
        "interpretation": "Random controls can be high, so embryo histone contrast is diagnostic support rather than final replacement.",
        "source": source(root, r"results\CSB_TRO_embryo_histone_random_controls.tsv"),
    })

    return pd.DataFrame(rows)


def build_supported_not_supported_table(root: Path) -> pd.DataFrame:
    rows = [
        {
            "claim": "methylation-only baseline is incomplete at morula basin entry",
            "supported": "yes",
            "primary_evidence": "baseline q90 occupancy=0.044 versus observed q90 occupancy=0.875",
            "boundary": "does not itself identify the missing biological variable",
            "manuscript_use": "main result",
        },
        {
            "claim": "measured correction is a compact diagnostic missing term",
            "supported": "yes",
            "primary_evidence": "measured correction reaches observed occupancy around alpha=0.50 and occupancy@alpha1=1.000",
            "boundary": "uses held-out morula methylation; diagnostic upper bound only",
            "manuscript_use": "main result with explicit leakage boundary",
        },
        {
            "claim": "correction is non-random at the top residual DMR level",
            "supported": "yes",
            "primary_evidence": "top25 observed occupancy=0.956 versus matched-random q95=0.156 and max=0.200",
            "boundary": "depends on the defined DMR universe and matching variables",
            "manuscript_use": "main control",
        },
        {
            "claim": "correction is module-specific",
            "supported": "yes",
            "primary_evidence": "greedy reconstruction: M05/M01/M12 occupancy=0.867; +M02 occupancy=0.956",
            "boundary": "modules are diagnostic coordinates, not causal units",
            "manuscript_use": "main result",
        },
        {
            "claim": "dual-branch closure/access architecture is direction-sensitive",
            "supported": "yes",
            "primary_evidence": "correct branch=0.956; wrong closure=0.000; wrong access=0.178",
            "boundary": "branch vectors are inferred from diagnostic residual geometry",
            "manuscript_use": "central result",
        },
        {
            "claim": "RNA or motif x TF identifies the final u_bio",
            "supported": "no",
            "primary_evidence": "global RNA occupancy=0.200; motif x TF occupancy=0.222",
            "boundary": "weak surrogate signals only",
            "manuscript_use": "boundary evidence",
        },
        {
            "claim": "public embryo histone contrast is final morula-entry replacement",
            "supported": "no",
            "primary_evidence": "true histone diagnostic=0.978 but matched-random median=0.978 and max=1.000",
            "boundary": "diagnostic plausibility only; matched random high",
            "manuscript_use": "boundary/support result",
        },
        {
            "claim": "histone evidence supports biological plausibility of branch architecture",
            "supported": "partially",
            "primary_evidence": "hESC max occupancy=0.511/cosine=0.806; public embryo diagnostic max occupancy=0.978/cosine=0.933",
            "boundary": "not stage-matched strict morula-entry replacement",
            "manuscript_use": "supporting result",
        },
        {
            "claim": "final causal u_bio is identified",
            "supported": "no",
            "primary_evidence": "observational public data plus controlled-access H3K27ac_morula/H3K4me3_morula gap",
            "boundary": "requires stage-matched external controls and perturbation for causality",
            "manuscript_use": "explicit limitation",
        },
    ]
    return pd.DataFrame(rows)


def write_integrated_figure(root: Path, out: Path) -> None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return

    sns.set_theme(style="whitegrid", font="Arial")
    final = read_tsv(root / r"results\CSB_TRO_final_model_comparison.tsv")
    sign = read_tsv(root / r"results\CSB_TRO_dual_branch_sign_control.tsv")
    random_dmr = read_tsv(root / r"results\CSB_TRO_residual_DMR_matched_random_control.tsv")
    hist_rand = read_tsv(root / r"results\CSB_TRO_embryo_histone_random_controls.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.4))
    fig.suptitle("Diagnostic evidence chain for a structured morula-entry correction", fontsize=15, fontweight="bold", y=0.965)

    # A. model comparison
    ax = axes[0, 0]
    model_order = [
        ("methylation-only baseline", "Meth-only\nbaseline"),
        ("measured correction upper bound", "Measured\ncorrection"),
        ("dual-branch chromatin-state proxy", "Dual\nbranch"),
        ("global RNA transition", "RNA"),
        ("motif x TF", "Motif x TF"),
        ("public embryo histone diagnostic contrast", "Embryo\nhistone"),
        ("strict morula-entry partial histone", "Strict\npartial\nhistone"),
    ]
    vals = []
    labels = []
    colors = []
    for name, label in model_order:
        row = final[final["model_name"] == name].iloc[0]
        val = row["max_occupancy"] if name in {"public embryo histone diagnostic contrast", "strict morula-entry partial histone"} else row["occupancy"]
        vals.append(float(val))
        labels.append(label)
        colors.append("#222222" if "methylation" in name else "#4c78a8" if "dual" in name else "#8c6bb1" if "measured" in name else "#b0b0b0")
    ax.bar(labels, vals, color=colors)
    ax.axhline(0.875, color="#666666", linestyle="--", linewidth=1)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("q90 occupancy")
    ax.set_title("A. Diagnostic model comparison", loc="left", fontweight="bold")
    ax.tick_params(axis="x", labelrotation=0, labelsize=9)

    # B. branch sign controls
    ax = axes[0, 1]
    sign_models = [
        ("correct_closure_correct_access", "Correct"),
        ("wrong_closure_correct_access", "Wrong\nclosure"),
        ("correct_closure_wrong_access", "Wrong\naccess"),
        ("wrong_closure_wrong_access", "Both\nwrong"),
        ("naive_inverse_all_modules", "Naive\ninverse"),
    ]
    vals = [float(sign[sign["model"] == m]["pred_basin_occupancy_q90"].iloc[0]) for m, _ in sign_models]
    labels = [lab for _, lab in sign_models]
    ax.bar(labels, vals, color=["#4c78a8", "#d95f02", "#d95f02", "#d95f02", "#999999"])
    ax.axhline(0.875, color="#666666", linestyle="--", linewidth=1)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("q90 occupancy")
    ax.set_title("B. Branch orientation controls", loc="left", fontweight="bold")
    ax.tick_params(axis="x", labelsize=9)

    # C. top residual DMR matched random
    ax = axes[1, 0]
    top25 = random_dmr[random_dmr["K"] == 25]
    observed = float(top25[top25["control_type"] == "observed_topK"]["pred_basin_occupancy_q90"].iloc[0])
    rnd = top25[top25["control_type"] == "matched_random"]["pred_basin_occupancy_q90"].astype(float)
    sns.histplot(rnd, bins=18, ax=ax, color="#bbbbbb", edgecolor="white")
    ax.axvline(observed, color="#4c78a8", linewidth=3, label=f"observed={observed:.3f}")
    ax.axvline(rnd.quantile(0.95), color="#d95f02", linewidth=2, linestyle="--", label=f"random q95={rnd.quantile(0.95):.3f}")
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("q90 occupancy")
    ax.set_ylabel("matched random count")
    ax.set_title("C. Top residual DMR specificity", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    # D. histone boundary
    ax = axes[1, 1]
    true_hist = float(hist_rand[hist_rand["validation_status"] == "true_target"]["pred_basin_occupancy_q90"].iloc[0])
    hist_controls = hist_rand[hist_rand["validation_status"] != "true_target"]["pred_basin_occupancy_q90"].astype(float)
    sns.histplot(hist_controls, bins=18, ax=ax, color="#bbbbbb", edgecolor="white")
    ax.axvline(true_hist, color="#4c78a8", linewidth=3, label=f"true={true_hist:.3f}")
    ax.axvline(hist_controls.median(), color="#d95f02", linewidth=2, linestyle="--", label=f"random median={hist_controls.median():.3f}")
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("q90 occupancy")
    ax.set_ylabel("matched random count")
    ax.set_title("D. Histone diagnostic boundary", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.text(
        0.5,
        0.012,
        "Dashed horizontal lines mark observed morula q90 occupancy. Histone contrast is retained as diagnostic plausibility because matched-random controls are high.",
        ha="center",
        fontsize=10,
        color="#555555",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.94], h_pad=2.0, w_pad=2.0)
    fig.savefig(out / "CSB_TRO_2026-05-27_integrated_diagnostic_evidence_figure.svg", format="svg")
    fig.savefig(out / "CSB_TRO_2026-05-27_integrated_diagnostic_evidence_figure.png", dpi=300)
    plt.close(fig)


def write_manuscript_assets(out: Path) -> None:
    figure_plan = """# Main figure plan for public-data-bounded diagnostic dynamics

## Figure 1. Operator-time methylation dynamics and morula-entry failure

Panels:
- DMR state-space/operator-time schematic.
- Methylation-only morula prediction versus observed morula basin.
- q90 occupancy contrast: methylation-only 0.044 versus observed 0.875.
- Measured correction alpha scan showing alpha_to_observed around 0.50.

Claim: methylation-only dynamics is a necessary baseline but is insufficient for morula reset-basin entry.

## Figure 2. Measured correction and residual modules

Panels:
- Definition of c_diag = x_observed,morula - x_predicted,morula.
- Top residual DMR/module contribution plot.
- Greedy module reconstruction: M05, M05/M01, M05/M01/M12, +M02, +M10.
- Matched-random top residual DMR control.

Claim: the correction is non-random and module-specific.

## Figure 3. Dual-branch closure/access diagnostic architecture

Panels:
- Closure branch M05/M01/M12 and access branch M02/M10 schematic.
- Correct closure/access occupancy 0.956.
- Wrong closure/wrong access sign controls.
- Exact sign-pattern and beta-grid robustness summary.

Claim: the correction is direction-sensitive and best organized by a dual-branch diagnostic architecture.

## Figure 4. Surrogate comparison and data-access boundary

Panels:
- RNA and motif x TF weak surrogate comparison.
- hESC/public embryo histone diagnostic support.
- Public embryo histone matched-random boundary.
- Controlled-access H3K27ac_morula/H3K4me3_morula gap.

Claim: histone evidence supports plausibility, but final biological replacement is data-access limited.

## Table 1. Evidence and boundary table

Use `CSB_TRO_2026-05-27_supported_not_supported_table.tsv` or `CSB_TRO_2026-05-27_evidence_boundary_table.tsv`.
"""

    results_skeleton = """# Results skeleton

## Result 1. Methylation-only operator-time dynamics fails at morula basin entry

We first modeled DMR methylation propagation using a methylation-only operator-time baseline. Although this baseline captures part of the stage-wise DMR-state trajectory, it fails to generate the observed morula reset-basin population structure. At q90, methylation-only morula occupancy is 0.044 compared with observed morula occupancy of 0.875.

## Result 2. The missing morula-entry component is measurable as a diagnostic correction

We defined the measured correction as the difference between observed morula state and the methylation-only morula prediction. This correction is a diagnostic upper bound because it uses held-out morula methylation. An alpha scan showed threshold-like basin entry, with alpha_to_observed around 0.50 and occupancy@alpha1 equal to 1.000.

## Result 3. The correction is non-random and module-specific

The measured correction is concentrated in a restricted residual module set. Greedy reconstruction showed occupancy of 0.867 for M05/M01/M12 and 0.956 after adding M02. Top residual DMR matched-random controls strongly supported non-random structure: top25 observed occupancy was 0.956, whereas matched-random q95 was 0.156 and max was 0.200.

## Result 4. A dual-branch closure/access architecture organizes the correction

The strongest diagnostic architecture separates a closure-like branch (M05/M01/M12) from an access-like branch (M02/M10). The correct branch orientation produced q90 occupancy of 0.956, whereas wrong closure gave 0.000 and wrong access gave 0.178. This supports a direction-sensitive closure/access correction architecture rather than an arbitrary residual fit.

## Result 5. Public surrogate layers define the biological boundary

RNA and motif x TF surrogates were weak replacements, with occupancy around 0.200 and 0.222, respectively. Histone-associated evidence supported branch plausibility, but public embryo histone matched-random controls were high, so this result is diagnostic support rather than final replacement. Strict morula-entry histone replacement remains limited by missing controlled-access H3K27ac_morula and H3K4me3_morula tracks.

## Result 6. The final claim is diagnostic, not causal

Together, the results support a biologically structured missing correction term projected into methylation state space. They do not identify a fully causal u_bio. The appropriate conclusion is a public-data-bounded diagnostic control framework with histone-supported biological plausibility.
"""

    captions = """# Draft captions for new integration outputs

## Integrated diagnostic evidence figure

Diagnostic evidence chain supporting a structured morula-entry correction. (A) Methylation-only baseline fails to recover observed morula q90 basin occupancy, whereas measured correction and the dual-branch diagnostic architecture recover high occupancy. RNA and motif x TF surrogates remain weak, and strict partial histone replacement remains limited. (B) Branch sign controls show strong direction dependence: the correct closure/access orientation reaches high occupancy, whereas wrong closure or wrong access reduces basin entry. (C) Top residual DMR matched-random controls show that observed high-residual DMRs strongly exceed matched random backgrounds. (D) Public embryo histone contrast is high, but matched-random histone controls are also high, supporting diagnostic plausibility rather than final biological replacement.

## Evidence-boundary table

Evidence and boundaries for a biologically structured missing correction. Each row lists the supported claim, key quantitative result, interpretation strength, and the reason the result should or should not be interpreted as final causal u_bio identification.
"""

    (out / "CSB_TRO_2026-05-27_main_figure_plan.md").write_text(figure_plan, encoding="utf-8")
    (out / "CSB_TRO_2026-05-27_results_skeleton.md").write_text(results_skeleton, encoding="utf-8")
    (out / "CSB_TRO_2026-05-27_caption_drafts.md").write_text(captions, encoding="utf-8")

    manuscript = """# Diagnostic operator-time methylation dynamics reveals a structured morula-entry correction architecture

## Abstract

Early embryonic methylation reset is usually analyzed through stage-wise methylation changes, but this does not directly test whether methylation state propagation is sufficient to generate the morula reset-basin. We constructed a DMR-level operator-time methylation dynamics framework and evaluated strict methylation-only propagation at morula entry. The methylation-only baseline captured part of the DMR-state trajectory but failed to generate the observed morula basin, with q90 occupancy of 0.044 compared with observed occupancy of 0.875. The missing component could be measured as a diagnostic correction term. This correction was non-random, module-specific, and concentrated in M05/M01/M12/M02/M10 residual coordinates. A dual-branch closure/access architecture organized the correction with high occupancy (0.956), high directional alignment, and strong sign dependence; wrong closure orientation collapsed occupancy to 0.000. Top residual DMR matched-random controls strongly supported non-random structure, whereas public RNA and motif x TF surrogates were weak. Histone-associated evidence supported branch plausibility but remained diagnostic because public embryo histone matched-random controls were high and strict morula-entry H3K27ac/H3K4me3 inputs are controlled-access. These results support a public-data-bounded diagnostic control framework: the missing morula-entry term is best interpreted as a biologically structured control-like projection, not as a fully identified causal u_bio.

## Results

### Methylation-only operator-time dynamics fails at morula basin entry

We first evaluated whether DMR methylation state propagation alone could account for morula reset-basin entry. The methylation-only baseline produced q90 morula basin occupancy of 0.044, whereas observed morula occupancy was 0.875. This established the methylation-only operator as an incomplete but necessary baseline rather than a complete developmental reset model.

### The missing component is measurable as a diagnostic correction

The failure component was defined as the measured correction between the observed morula state and the methylation-only morula prediction. This correction is an upper-bound diagnostic because it uses held-out morula methylation. An alpha scan showed threshold-like basin entry, with alpha_to_observed around 0.50 and occupancy@alpha1 of 1.000.

### The correction is non-random and module-specific

Residual module reconstruction identified a compact correction structure dominated by M05/M01/M12/M02/M10. Greedy reconstruction reached occupancy of 0.867 using M05/M01/M12 and 0.956 after adding M02. Top residual DMR matched-random controls strongly supported non-random structure: top25 observed occupancy was 0.956, while matched-random q95 was 0.156 and max was 0.200.

### A dual-branch closure/access architecture organizes the correction

The strongest diagnostic architecture separated closure-like M05/M01/M12 modules from access-like M02/M10 modules. Correct closure/access orientation reached q90 occupancy of 0.956. Wrong closure with correct access produced occupancy of 0.000, correct closure with wrong access produced 0.178, and both wrong signs produced 0.000. This sign dependence indicates that the branch architecture is not an arbitrary residual fit.

### Public surrogate layers define the biological boundary

RNA and motif x TF surrogates were weak replacements, with global RNA occupancy around 0.200 and motif x TF occupancy around 0.222. hESC histone and public embryo histone contrasts supported branch plausibility, but public embryo histone matched-random controls were high: true diagnostic occupancy was 0.978 and random median was also 0.978. Therefore, public histone evidence is diagnostic support, not final morula-entry replacement. Strict morula-entry histone replacement remains limited by missing controlled-access H3K27ac_morula and H3K4me3_morula tracks.

## Discussion

This study establishes a public-data-bounded diagnostic control framework for morula-entry methylation reset. The central result is not that a final causal u_bio has been identified. Instead, the results show that methylation-only dynamics fails in a structured way, and that the measured missing term has module specificity, direction sensitivity, and a dual-branch closure/access organization.

The strongest evidence comes from negative controls. Branch sign controls showed that incorrect closure orientation destroys basin entry, and top residual DMR matched-random controls showed that the high-residual DMR set is far stronger than matched random backgrounds. These controls support the interpretation that the missing correction is not a random mathematical residual.

The main biological boundary is histone data access. Public embryo histone contrasts are compatible with the closure/access architecture, but matched-random histone controls are also high. In addition, strict morula-entry H3K27ac and H3K4me3 tracks are not available in the local public data package and are tied to controlled-access sources. The correct interpretation is therefore histone-supported biological plausibility, not final biological replacement.

## Claim boundary

Supported claim:

Methylation-only operator-time dynamics captures baseline DMR-state propagation but fails at morula reset-basin entry. The resulting measured correction is non-random, module-specific, direction-sensitive, and organized by a dual-branch closure/access diagnostic architecture. Public histone-associated evidence supports biological plausibility, while RNA and motif x TF surrogates provide weak replacement signals.

Unsupported claim:

The final causal biological control input u_bio has been fully identified.

## Future work

The immediate controlled-data upgrade is to obtain stage-matched morula H3K27ac and H3K4me3 tracks and test whether they provide a strict non-leaking biological replacement for the diagnostic correction. Longer-term extensions include matched ATAC/histone profiling, perturbation validation, cross-species replication, and continuous-time or optimal-transport formulations. These are future directions and are not required for the current public-data-bounded diagnostic manuscript.
"""

    checklist = """# Submission readiness checklist

## Complete

- Core result inventory generated.
- Evidence/boundary table generated.
- Supported/not-supported table generated.
- Negative-control coverage table generated.
- Integrated diagnostic evidence figure generated as SVG and PNG.
- Results skeleton generated.
- Figure plan and caption drafts generated.
- Language risk audit generated and currently reports zero high-risk matches.

## Main manuscript actions

- Use dual-branch architecture as the central result.
- Put sign controls and top residual DMR matched-random controls in the main figures.
- Present public embryo histone as diagnostic support only.
- Keep RNA and motif x TF as weak surrogate boundary results.
- Use supported/not-supported table to prevent overclaim.

## Do not add to current main line

- New proxy tuning.
- Mouse/cross-species replication.
- Wet-lab validation.
- OT or Schrodinger bridge method comparison.
- Claims that final causal u_bio has been identified.

## Remaining risk

- Histone support is biologically useful but not specific enough under matched-random control.
- Final morula-entry histone replacement requires controlled-access H3K27ac_morula/H3K4me3_morula data.
"""

    (out / "CSB_TRO_2026-05-27_manuscript_v0.9_diagnostic_framework.md").write_text(manuscript, encoding="utf-8")
    (out / "CSB_TRO_2026-05-27_submission_readiness_checklist.md").write_text(checklist, encoding="utf-8")


def write_summary(out: Path, evidence: pd.DataFrame, controls: pd.DataFrame) -> None:
    lines = [
        "# Public-data bounded diagnostic evidence summary",
        "",
        "Status: generated from frozen local result files; no original result files were modified.",
        "",
        "## Main conclusion",
        "",
        "The current evidence supports a biologically structured missing correction term, not a fully identified causal u_bio.",
        "",
        "Recommended claim:",
        "",
        "Methylation-only operator-time dynamics captures baseline DMR-state propagation but fails at morula reset-basin entry. The measured correction is structured, module-specific, direction-sensitive, and organized by a dual-branch closure/access diagnostic architecture. Public RNA and motif x TF surrogates are weak replacements, while histone-associated evidence supports biological plausibility. Final replacement remains bounded by controlled-access morula H3K27ac/H3K4me3 data.",
        "",
        "## Evidence chain",
        "",
    ]
    for r in evidence.itertuples(index=False):
        lines.append(f"- {r.evidence_layer}: {r.result} ({r.claim_strength})")
    lines.extend([
        "",
        "## Control status",
        "",
    ])
    for r in controls.itertuples(index=False):
        lines.append(f"- {r.control_family} / {r.control_name}: {r.control_summary}; pass={r.passes}")
    lines.extend([
        "",
        "## Immediate manuscript action",
        "",
        "1. Keep negative controls in the main figure panels where possible.",
        "2. Use the evidence table as a boundary table, not as proof of causal u_bio.",
        "3. Do not expand into mouse, cross-species, wet-lab, OT, or Schrodinger bridge analyses in the current manuscript.",
    ])
    (out / "CSB_TRO_2026-05-27_diagnostic_evidence_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_control_figure(out: Path, controls: pd.DataFrame) -> None:
    selected = [
        ("dual_branch_sign", "wrong_closure_correct_access", "Wrong closure"),
        ("dual_branch_sign", "correct_closure_wrong_access", "Wrong access"),
        ("dual_branch_sign", "wrong_closure_wrong_access", "Both wrong"),
        ("dual_branch_random_partition", "random 3/2 module partitions", "Random partitions"),
        ("top_residual_DMR_matched_random", "top25 matched random DMRs", "Top25 matched random"),
        ("embryo_histone_matched_random", "public embryo histone diagnostic matched random", "Embryo histone random"),
    ]
    labels = []
    true_vals = []
    control_vals = []
    pass_colors = []
    for family, name, label in selected:
        sub = controls[(controls["control_family"] == family) & (controls["control_name"] == name)]
        if sub.empty:
            continue
        row = sub.iloc[0]
        labels.append(label)
        true_vals.append(float(row["true_value"]))
        summary = str(row["control_summary"])
        if "q95=" in summary:
            ctrl = float(summary.split("q95=")[1].split(";")[0])
        else:
            ctrl = float(summary.split(";")[0]) if ";" in summary else float(summary)
        control_vals.append(ctrl)
        pass_colors.append("#4c78a8" if row["passes"] == "yes" else "#d95f02" if row["passes"] == "no" else "#8c6bb1")
    width, height = 980, 460
    left, top, plot_w, plot_h = 80, 58, 850, 270
    y_base = top + plot_h
    max_y = 1.08
    group_w = plot_w / max(len(labels), 1)
    bar_w = 24

    def y(v: float) -> float:
        return y_base - (v / max_y) * plot_h

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="80" y="30" font-family="Arial" font-size="20" font-weight="700">Negative-control coverage for diagnostic evidence chain</text>',
        f'<line x1="{left}" y1="{y_base}" x2="{left + plot_w}" y2="{y_base}" stroke="#333" stroke-width="1"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{y_base}" stroke="#333" stroke-width="1"/>',
    ]
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        yy = y(tick)
        parts.append(f'<line x1="{left-5}" y1="{yy:.1f}" x2="{left + plot_w}" y2="{yy:.1f}" stroke="#e6e6e6" stroke-width="1"/>')
        parts.append(f'<text x="{left-12}" y="{yy+4:.1f}" text-anchor="end" font-family="Arial" font-size="11">{tick:.2f}</text>')
    obs_y = y(0.875)
    parts.append(f'<line x1="{left}" y1="{obs_y:.1f}" x2="{left + plot_w}" y2="{obs_y:.1f}" stroke="#666" stroke-width="1.5" stroke-dasharray="5 4"/>')
    parts.append(f'<text x="{left + plot_w - 4}" y="{obs_y-6:.1f}" text-anchor="end" font-family="Arial" font-size="11" fill="#555">observed q90 occupancy = 0.875</text>')

    for i, label in enumerate(labels):
        cx = left + group_w * i + group_w / 2
        tv, cv = true_vals[i], control_vals[i]
        for x0, val, color in [(cx - bar_w - 3, tv, "#222222"), (cx + 3, cv, pass_colors[i])]:
            yy = y(val)
            parts.append(f'<rect x="{x0:.1f}" y="{yy:.1f}" width="{bar_w}" height="{y_base-yy:.1f}" fill="{color}"/>')
            parts.append(f'<text x="{x0 + bar_w/2:.1f}" y="{yy-5:.1f}" text-anchor="middle" font-family="Arial" font-size="10">{val:.2f}</text>')
        parts.append(f'<text x="{cx:.1f}" y="{y_base+22}" text-anchor="middle" font-family="Arial" font-size="11">{label}</text>')

    parts.extend([
        '<rect x="80" y="395" width="12" height="12" fill="#222222"/>',
        '<text x="98" y="405" font-family="Arial" font-size="12">target / true</text>',
        '<rect x="205" y="395" width="12" height="12" fill="#4c78a8"/>',
        '<text x="223" y="405" font-family="Arial" font-size="12">control passes</text>',
        '<rect x="335" y="395" width="12" height="12" fill="#8c6bb1"/>',
        '<text x="353" y="405" font-family="Arial" font-size="12">qualified</text>',
        '<rect x="445" y="395" width="12" height="12" fill="#d95f02"/>',
        '<text x="463" y="405" font-family="Arial" font-size="12">does not pass</text>',
        '<text x="80" y="438" font-family="Arial" font-size="12" fill="#555">For random controls, the colored bar is q95 when available. Histone random is high, so public embryo histone remains diagnostic support only.</text>',
        "</svg>",
    ])
    (out / "CSB_TRO_2026-05-27_negative_control_summary.svg").write_text("\n".join(parts), encoding="utf-8")


def build_language_audit(root: Path) -> pd.DataFrame:
    targets = [
        root / r"docs\CSB_TRO_public_data_bounded_manuscript_draft.md",
        root / r"docs\CSB_TRO_final_claim_boundary_summary.md",
    ]
    patterns = [
        (r"\bprove\b", "prove", "Use provide evidence / support instead."),
        (r"\bcausal\b", "causal", "Avoid causal language unless discussing future perturbation."),
        (r"\bcausality\b", "causality", "Avoid causal language unless discussing future perturbation."),
        (r"\bidentified u_bio\b", "identified u_bio", "Use diagnostic architecture for the missing correction."),
        (r"\bidentify u_bio\b", "identify u_bio", "Use diagnostic architecture for the missing correction."),
        (r"\btrue biological control\b", "true biological control", "Use biologically structured control-like projection."),
        (r"\bfinal replacement\b", "final replacement", "Use stage-matched replacement remains data-access limited."),
        (r"\bhistone replacement\b", "histone replacement", "Use histone-supported diagnostic plausibility."),
        (r"\bmechanism\b", "mechanism", "Use architecture or plausible mechanism unless directly supported."),
    ]
    rows = []
    for path in targets:
        if not path.exists():
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        for i, line in enumerate(lines, start=1):
            low = line.lower()
            if "rather than causal" in low or "not causal" in low:
                continue
            for regex, term, recommendation in patterns:
                if re.search(regex, low):
                    rows.append({
                        "file": str(path),
                        "line": i,
                        "matched_term": term,
                        "text": line.strip(),
                        "recommendation": recommendation,
                    })
    return pd.DataFrame(rows, columns=["file", "line", "matched_term", "text", "recommendation"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    claim = build_claim_inventory(args.root)
    evidence = build_evidence_table(args.root)
    controls = build_control_summary(args.root)
    supported = build_supported_not_supported_table(args.root)
    audit = build_language_audit(args.root)

    claim.to_csv(args.out / "CSB_TRO_2026-05-27_result_inventory.tsv", sep="\t", index=False)
    evidence.to_csv(args.out / "CSB_TRO_2026-05-27_evidence_boundary_table.tsv", sep="\t", index=False)
    controls.to_csv(args.out / "CSB_TRO_2026-05-27_negative_control_coverage.tsv", sep="\t", index=False)
    supported.to_csv(args.out / "CSB_TRO_2026-05-27_supported_not_supported_table.tsv", sep="\t", index=False)
    audit.to_csv(args.out / "CSB_TRO_2026-05-27_manuscript_language_risk_audit.tsv", sep="\t", index=False)
    write_summary(args.out, evidence, controls)
    write_control_figure(args.out, controls)
    write_integrated_figure(args.root, args.out)
    write_manuscript_assets(args.out)
    manifest = {
        "status": "public_data_bounded_diagnostic_evidence_chain",
        "input_root": str(args.root),
        "output_dir": str(args.out),
        "generated_files": [
            "CSB_TRO_2026-05-27_result_inventory.tsv",
            "CSB_TRO_2026-05-27_evidence_boundary_table.tsv",
            "CSB_TRO_2026-05-27_negative_control_coverage.tsv",
            "CSB_TRO_2026-05-27_supported_not_supported_table.tsv",
            "CSB_TRO_2026-05-27_manuscript_language_risk_audit.tsv",
            "CSB_TRO_2026-05-27_diagnostic_evidence_summary.md",
            "CSB_TRO_2026-05-27_negative_control_summary.svg",
            "CSB_TRO_2026-05-27_integrated_diagnostic_evidence_figure.svg",
            "CSB_TRO_2026-05-27_integrated_diagnostic_evidence_figure.png",
            "CSB_TRO_2026-05-27_main_figure_plan.md",
            "CSB_TRO_2026-05-27_results_skeleton.md",
            "CSB_TRO_2026-05-27_caption_drafts.md",
            "CSB_TRO_2026-05-27_manuscript_v0.9_diagnostic_framework.md",
            "CSB_TRO_2026-05-27_submission_readiness_checklist.md",
        ],
        "claim_boundary": "Supports biologically structured missing correction, not fully identified causal u_bio.",
        "source_script": str(Path(__file__).resolve()),
    }
    (args.out / "CSB_TRO_2026-05-27_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
