#!/usr/bin/env python3
"""Summarize external transcriptomic validation for the MN manuscript.

This script compares the discovery glomerular PMN contrast (GSE108113) with
an independent glomerular MGN/control contrast from GSE99340-GPL19184. It also
checks whether the network-medicine ranking is preserved when the disease set is
rebuilt from the external transcriptomic contrast. It does not overwrite
primary discovery-cohort outputs.
"""
from __future__ import annotations

import pickle
import random
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact, pearsonr, spearmanr

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data"
RESULTS = BASE / "results"
EXTERNAL = RESULTS / "external"
TABLES = BASE / "tables"
FIGURES = BASE / "figures"
EXTERNAL.mkdir(exist_ok=True, parents=True)
TABLES.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

DISCOVERY_DE = RESULTS / "DE_MGN_vs_Normal.csv"
EXTERNAL_DE = EXTERNAL / "GSE99340_GPL19184_glom_MGN_vs_control_DE.csv"
EXTERNAL_META = EXTERNAL / "GSE99340_GPL19184_glom_MGN_vs_control_metadata.csv"
GWAS = {"PLA2R1", "NFKB1", "IRF4", "HLA-DQA1", "HLA-DRB1"}
SEED = 20260519


def clean_gene_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().replace({"": np.nan, "nan": np.nan, "NA": np.nan})


def prepare_de(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"Gene", "logFC", "P.Value", "adj.P.Val"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    df = df.copy()
    df["Gene"] = clean_gene_series(df["Gene"])
    df = df.dropna(subset=["Gene", "logFC", "P.Value", "adj.P.Val"])
    df = df.sort_values("P.Value").drop_duplicates("Gene", keep="first")
    df = df.rename(
        columns={
            "logFC": f"logFC_{label}",
            "P.Value": f"P.Value_{label}",
            "adj.P.Val": f"adj.P.Val_{label}",
        }
    )
    keep = ["Gene", f"logFC_{label}", f"P.Value_{label}", f"adj.P.Val_{label}"]
    for optional in ["Probe", "Entrez"]:
        if optional in df.columns:
            keep.append(optional)
    return df[keep]


def fmt_float(x: float, digits: int = 4) -> str:
    if pd.isna(x):
        return "NA"
    if abs(x) < 1e-3 and x != 0:
        return f"{x:.2e}"
    return f"{x:.{digits}f}"


def top_set(df: pd.DataFrame, p_col: str, n: int) -> set[str]:
    return set(df.sort_values(p_col).head(n)["Gene"])


def load_lcc_graph() -> nx.Graph:
    graph = nx.read_graphml(DATA / "mn_ppi_network.graphml").to_undirected()
    for _, _, data in graph.edges(data=True):
        data["weight"] = float(data.get("weight", 1.0))
    lcc = max(nx.connected_components(graph), key=len)
    return graph.subgraph(lcc).copy()


def mean_distance_to_disease(graph: nx.Graph, source_nodes: list[str], disease_nodes: set[str]) -> float:
    source = [node for node in source_nodes if node in graph]
    disease = [node for node in disease_nodes if node in graph]
    if len(source) < 2 or len(disease) < 2:
        return np.nan
    disease_distance = nx.multi_source_dijkstra_path_length(graph, disease, weight=None)
    vals = [disease_distance[node] for node in source if node in disease_distance]
    return float(np.mean(vals)) if vals else np.nan


def degree_bins(graph: nx.Graph, bin_width: int = 5) -> dict[int, list[str]]:
    bins: dict[int, list[str]] = {}
    for node, degree in graph.degree():
        bins.setdefault(int(degree // bin_width), []).append(node)
    return bins


def sample_degree_matched(graph: nx.Graph, targets: list[str], bins: dict[int, list[str]], rng: random.Random, bin_width: int = 5) -> list[str]:
    all_nodes = list(graph.nodes())
    sampled = []
    for target in targets:
        key = int(graph.degree(target) // bin_width)
        candidates = [node for node in bins.get(key, []) if node != target]
        if len(candidates) < 2:
            candidates = []
            for delta in range(1, 6):
                candidates.extend(node for node in bins.get(key - delta, []) if node != target)
                candidates.extend(node for node in bins.get(key + delta, []) if node != target)
                if len(candidates) >= 2:
                    break
        if not candidates:
            candidates = all_nodes
        sampled.append(rng.choice(candidates))
    return sampled


def proximity_result(
    graph: nx.Graph,
    drug_targets: set[str],
    disease_nodes: set[str],
    *,
    mode: str,
    n_perm: int,
    seed: int,
) -> dict[str, object]:
    targets = sorted(set(drug_targets) & set(graph.nodes()))
    disease = set(disease_nodes) & set(graph.nodes())
    if len(targets) < 2 or len(disease) < 2:
        return {
            "Targets in network": len(targets),
            "Disease proteins in network": len(disease),
            "Observed distance": np.nan,
            "Permutation mean distance": np.nan,
            "Proximity z-score": np.nan,
        }
    obs = mean_distance_to_disease(graph, targets, disease)
    rng = random.Random(seed)
    if mode == "size-matched":
        nodes = list(graph.nodes())
        perm = [mean_distance_to_disease(graph, rng.sample(nodes, len(targets)), disease) for _ in range(n_perm)]
    elif mode == "degree-binned":
        bins = degree_bins(graph)
        perm = [mean_distance_to_disease(graph, sample_degree_matched(graph, targets, bins, rng), disease) for _ in range(n_perm)]
    else:
        raise ValueError(f"Unknown permutation mode: {mode}")
    perm = np.asarray(perm, dtype=float)
    sd = float(np.nanstd(perm))
    z = float((obs - np.nanmean(perm)) / sd) if sd > 0 else np.nan
    return {
        "Targets in network": len(targets),
        "Disease proteins in network": len(disease),
        "Observed distance": obs,
        "Permutation mean distance": float(np.nanmean(perm)),
        "Proximity z-score": z,
    }


def external_network_proximity(external_full: pd.DataFrame) -> pd.DataFrame:
    graph = load_lcc_graph()
    with open(DATA / "drug_targets.pkl", "rb") as handle:
        drug_targets = pickle.load(handle)

    external_sig = set(external_full.loc[external_full["adj.P.Val_external"] < 0.05, "Gene"])
    external_nominal = set(external_full.loc[external_full["P.Value_external"] < 0.05, "Gene"])
    top1000 = set(external_full.sort_values("P.Value_external").head(1000)["Gene"])
    disease_sets = {
        "External FDR<0.05 DEGs": external_sig,
        "External nominal P<0.05 DEGs": external_nominal,
        "External top 1000 genes": top1000,
        "External FDR<0.05 DEGs + GWAS": external_sig | GWAS,
    }

    rows = []
    counter = 0
    for disease_name, disease_nodes in disease_sets.items():
        for drug in ["RTX", "TAC", "CTX"]:
            for mode, n_perm in [("size-matched", 1000), ("degree-binned", 300)]:
                counter += 1
                res = proximity_result(
                    graph,
                    set(drug_targets[drug]),
                    disease_nodes,
                    mode=mode,
                    n_perm=n_perm,
                    seed=SEED + counter,
                )
                rows.append(
                    {
                        "Disease set": disease_name,
                        "Drug": drug,
                        "Permutation mode": mode,
                        **res,
                        "Interpretation": "Proximal" if pd.notna(res["Proximity z-score"]) and res["Proximity z-score"] < -1.5 else "Not proximal by z < -1.5",
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(EXTERNAL / "external_validation_network_proximity.csv", index=False)
    return out


def main() -> None:
    discovery = prepare_de(DISCOVERY_DE, "discovery")
    external = prepare_de(EXTERNAL_DE, "external")
    meta = pd.read_csv(EXTERNAL_META)

    merged = discovery.merge(external, on="Gene", how="inner")
    merged["same_direction"] = np.sign(merged["logFC_discovery"]) == np.sign(merged["logFC_external"])
    merged["discovery_FDR05"] = merged["adj.P.Val_discovery"] < 0.05
    merged["external_FDR05"] = merged["adj.P.Val_external"] < 0.05
    merged["external_nominal05"] = merged["P.Value_external"] < 0.05
    merged["both_FDR05"] = merged["discovery_FDR05"] & merged["external_FDR05"]
    merged["both_FDR05_same_direction"] = merged["both_FDR05"] & merged["same_direction"]
    merged.to_csv(EXTERNAL / "external_validation_merged_gene_level.csv", index=False)

    shared_n = len(merged)
    external_gene_level_fdr = int((external["adj.P.Val_external"] < 0.05).sum())
    discovery_sig = set(merged.loc[merged["discovery_FDR05"], "Gene"])
    external_sig = set(merged.loc[merged["external_FDR05"], "Gene"])
    overlap_sig = discovery_sig & external_sig
    overlap_sig_same = set(merged.loc[merged["both_FDR05_same_direction"], "Gene"])

    pearson_all = pearsonr(merged["logFC_discovery"], merged["logFC_external"])
    spearman_all = spearmanr(merged["logFC_discovery"], merged["logFC_external"])

    disc_sig_df = merged[merged["discovery_FDR05"]].copy()
    disc_sig_same_direction = int(disc_sig_df["same_direction"].sum())
    disc_sig_same_prop = disc_sig_same_direction / len(disc_sig_df) if len(disc_sig_df) else np.nan
    disc_sig_binom = binomtest(disc_sig_same_direction, len(disc_sig_df), 0.5, alternative="greater") if len(disc_sig_df) else None

    disc_sig_ext_nom_df = disc_sig_df[disc_sig_df["external_nominal05"]].copy()
    disc_sig_ext_nom_same = int(disc_sig_ext_nom_df["same_direction"].sum())
    disc_sig_ext_nom_prop = disc_sig_ext_nom_same / len(disc_sig_ext_nom_df) if len(disc_sig_ext_nom_df) else np.nan
    disc_sig_ext_nom_binom = binomtest(disc_sig_ext_nom_same, len(disc_sig_ext_nom_df), 0.5, alternative="greater") if len(disc_sig_ext_nom_df) else None

    a = len(overlap_sig)
    b = len(discovery_sig - external_sig)
    c = len(external_sig - discovery_sig)
    d = shared_n - a - b - c
    fisher_or, fisher_p = fisher_exact([[a, b], [c, d]], alternative="greater")

    top_rows = []
    for n in [100, 300, 500, 1000]:
        n_eff = min(n, len(discovery), len(external))
        disc_top = top_set(discovery, "P.Value_discovery", n_eff)
        ext_top = top_set(external, "P.Value_external", n_eff)
        universe = set(merged["Gene"])
        disc_top &= universe
        ext_top &= universe
        overlap = disc_top & ext_top
        a_top = len(overlap)
        b_top = len(disc_top - ext_top)
        c_top = len(ext_top - disc_top)
        d_top = len(universe) - a_top - b_top - c_top
        top_or, top_p = fisher_exact([[a_top, b_top], [c_top, d_top]], alternative="greater")
        top_rows.append(
            {
                "Top N by nominal P": n_eff,
                "Discovery top genes in shared universe": len(disc_top),
                "External top genes in shared universe": len(ext_top),
                "Overlap": a_top,
                "Jaccard": a_top / len(disc_top | ext_top) if (disc_top | ext_top) else np.nan,
                "Fisher odds ratio": top_or,
                "Fisher P value": top_p,
            }
        )
    top_df = pd.DataFrame(top_rows)
    top_df.to_csv(EXTERNAL / "external_validation_top_gene_overlap.csv", index=False)

    key = merged.sort_values("P.Value_discovery").head(200).copy()
    key["replicated_nominal_same_direction"] = key["external_nominal05"] & key["same_direction"]
    key.to_csv(EXTERNAL / "external_validation_key_discovery_genes.csv", index=False)

    sample_summary = meta.groupby(["compartment", "disease_group"]).size().reset_index(name="n")
    sample_summary.to_csv(EXTERNAL / "external_validation_sample_summary.csv", index=False)

    external_full = external.rename(
        columns={
            "logFC_external": "logFC_external",
            "P.Value_external": "P.Value_external",
            "adj.P.Val_external": "adj.P.Val_external",
        }
    )
    prox_df = external_network_proximity(external_full)
    prox_focus = prox_df[
        (prox_df["Disease set"] == "External FDR<0.05 DEGs + GWAS")
        & (prox_df["Permutation mode"].isin(["size-matched", "degree-binned"]))
    ].copy()

    summary_rows = [
        {
            "Validation component": "External cohort",
            "Metric": "Dataset and contrast",
            "Value": "GSE99340-GPL19184 glomeruli, MGN vs control",
            "Interpretation": "Independent microarray glomerular validation cohort",
        },
        {
            "Validation component": "External cohort",
            "Metric": "Samples analyzed",
            "Value": f"MGN={int((meta['disease_group'] == 'MGN').sum())}; control={int((meta['disease_group'] == 'Control').sum())}",
            "Interpretation": "Filtered to glomerular compartment only",
        },
        {
            "Validation component": "Gene overlap",
            "Metric": "Shared mapped genes",
            "Value": str(shared_n),
            "Interpretation": "Genes available for cross-platform concordance testing",
        },
        {
            "Validation component": "Differential expression",
            "Metric": "External gene-level DEGs at FDR < 0.05",
            "Value": str(external_gene_level_fdr),
            "Interpretation": "Independent glomerular MGN transcriptomic signal before restricting to cross-platform shared genes",
        },
        {
            "Validation component": "Differential expression",
            "Metric": "External FDR<0.05 DEGs represented in shared discovery-external gene universe",
            "Value": str(len(external_sig)),
            "Interpretation": "External FDR-significant genes available for cross-platform overlap testing",
        },
        {
            "Validation component": "Global concordance",
            "Metric": "Pearson correlation of logFC across shared genes",
            "Value": f"r={fmt_float(pearson_all.statistic)}; P={fmt_float(pearson_all.pvalue)}",
            "Interpretation": "Positive cross-cohort effect-size concordance" if pearson_all.statistic > 0 else "No positive global concordance",
        },
        {
            "Validation component": "Global concordance",
            "Metric": "Spearman correlation of logFC across shared genes",
            "Value": f"rho={fmt_float(spearman_all.statistic)}; P={fmt_float(spearman_all.pvalue)}",
            "Interpretation": "Positive rank-order concordance" if spearman_all.statistic > 0 else "No positive rank-order concordance",
        },
        {
            "Validation component": "Discovery DEG concordance",
            "Metric": "Discovery FDR<0.05 genes with same external direction",
            "Value": f"{disc_sig_same_direction}/{len(disc_sig_df)} ({fmt_float(disc_sig_same_prop * 100, 1)}%); binomial P={fmt_float(disc_sig_binom.pvalue if disc_sig_binom else np.nan)}",
            "Interpretation": "Directionality of discovery transcriptomic seed genes in external cohort",
        },
        {
            "Validation component": "Discovery DEG concordance",
            "Metric": "Discovery FDR<0.05 and external nominal P<0.05 genes with same direction",
            "Value": f"{disc_sig_ext_nom_same}/{len(disc_sig_ext_nom_df)} ({fmt_float(disc_sig_ext_nom_prop * 100, 1)}%); binomial P={fmt_float(disc_sig_ext_nom_binom.pvalue if disc_sig_ext_nom_binom else np.nan)}",
            "Interpretation": "Directionality among discovery genes showing at least nominal external evidence",
        },
        {
            "Validation component": "DEG overlap",
            "Metric": "Discovery and external FDR<0.05 overlap",
            "Value": f"{a} genes; OR={fmt_float(fisher_or)}; Fisher P={fmt_float(fisher_p)}",
            "Interpretation": "Strict cross-platform DEG overlap was not enriched" if fisher_p > 0.05 else "Enrichment of discovery DEGs among external DEGs",
        },
        {
            "Validation component": "DEG overlap",
            "Metric": "FDR<0.05 overlap with same direction",
            "Value": f"{len(overlap_sig_same)}/{a} overlapping genes ({fmt_float((len(overlap_sig_same) / a * 100) if a else np.nan, 1)}%)",
            "Interpretation": "Direction-consistent replicated DEG subset",
        },
    ]

    for _, row in top_df.iterrows():
        summary_rows.append(
            {
                "Validation component": "Top-ranked DEG overlap",
                "Metric": f"Top {int(row['Top N by nominal P'])} overlap",
                "Value": f"overlap={int(row['Overlap'])}; Jaccard={fmt_float(row['Jaccard'])}; OR={fmt_float(row['Fisher odds ratio'])}; Fisher P={fmt_float(row['Fisher P value'])}",
                "Interpretation": "Top-ranked discovery genes are enriched among top-ranked external genes" if row["Fisher P value"] < 0.05 else "No significant top-ranked overlap enrichment",
            }
        )

    for _, row in prox_focus.iterrows():
        summary_rows.append(
            {
                "Validation component": "External-module network proximity",
                "Metric": f"{row['Drug']} proximity to external FDR<0.05 DEG+GWAS module ({row['Permutation mode']})",
                "Value": f"z={fmt_float(row['Proximity z-score'])}; observed={fmt_float(row['Observed distance'])}; perm.mean={fmt_float(row['Permutation mean distance'])}",
                "Interpretation": row["Interpretation"],
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(EXTERNAL / "external_validation_concordance.csv", index=False)
    summary.to_csv(TABLES / "Supplementary_Table_S8_External_validation.csv", index=False)

    # Compact figure for supplementary use.
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), dpi=300)
    ax = axes[0]
    sig = merged["discovery_FDR05"] | merged["external_FDR05"]
    ax.scatter(merged.loc[~sig, "logFC_discovery"], merged.loc[~sig, "logFC_external"], s=6, alpha=0.18, color="#999999", linewidths=0)
    ax.scatter(merged.loc[sig, "logFC_discovery"], merged.loc[sig, "logFC_external"], s=8, alpha=0.5, color="#1f77b4", linewidths=0)
    ax.axhline(0, color="black", lw=0.6)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("GSE108113 logFC (PMN vs control)")
    ax.set_ylabel("GSE99340 logFC (MGN vs control)")
    ax.set_title(f"Cross-cohort logFC concordance\nr={pearson_all.statistic:.2f}, rho={spearman_all.statistic:.2f}")

    ax = axes[1]
    plot_top = top_df.copy()
    ax.bar(plot_top["Top N by nominal P"].astype(str), plot_top["Overlap"], color="#4c78a8")
    ax.set_xlabel("Top N genes by nominal P")
    ax.set_ylabel("Overlap count")
    ax.set_title("Top-ranked DEG overlap")
    for i, row in plot_top.iterrows():
        ax.text(i, row["Overlap"], f"P={row['Fisher P value']:.1e}", ha="center", va="bottom", fontsize=7, rotation=45)

    ax = axes[2]
    plot_prox = prox_focus[prox_focus["Permutation mode"] == "degree-binned"].copy()
    ax.bar(plot_prox["Drug"], plot_prox["Proximity z-score"], color=["#1f77b4", "#ff7f0e", "#2ca02c"])
    ax.axhline(-1.5, color="black", ls="--", lw=0.8)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_ylabel("Degree-binned proximity z-score")
    ax.set_title("External DEG+GWAS module proximity")
    for i, (_, row) in enumerate(plot_prox.iterrows()):
        z = row["Proximity z-score"]
        ax.text(i, z, f"{z:.2f}", ha="center", va="bottom" if z >= 0 else "top", fontsize=8)

    fig.tight_layout()
    fig.savefig(FIGURES / "Supplementary_Figure_S2_External_Validation.png", bbox_inches="tight")
    plt.close(fig)

    print("External validation summary")
    print(summary.to_string(index=False))
    print(f"Saved {TABLES / 'Supplementary_Table_S8_External_validation.csv'}")
    print(f"Saved {FIGURES / 'Supplementary_Figure_S2_External_Validation.png'}")


if __name__ == "__main__":
    main()
