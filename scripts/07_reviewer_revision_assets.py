#!/usr/bin/env python3
"""Generate reviewer-facing tables, robustness checks, and merged figures."""

from pathlib import Path
import math
import pickle
import random

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


BASE = Path("/workspaces/platform/MN")
DATA = BASE / "data"
RESULTS = BASE / "results"
FIGURES = BASE / "figures"
TABLES = BASE / "tables"
SUBFIG = FIGURES / "submission"

TABLES.mkdir(exist_ok=True)
SUBFIG.mkdir(parents=True, exist_ok=True)


def read_gene_count(path):
    return len(pd.read_csv(path)["Gene"].dropna().unique())


def write_tables():
    datasets = pd.DataFrame(
        [
            {
                "Resource": "GSE108113",
                "Data type": "Bulk glomerular transcriptomics",
                "Samples used": "87 PMN; 11 living donor controls",
                "Platform/source": "Affymetrix Human Gene 2.1 ST Array; NEPTUNE",
                "Role in analysis": "Differential expression and transcriptomic disease module seed genes",
            },
            {
                "Resource": "GSE99340-GPL19184",
                "Data type": "External bulk glomerular transcriptomics",
                "Samples used": "21 MGN; 8 control glomerular samples",
                "Platform/source": "GEO SuperSeries microarray platform GPL19184",
                "Role in analysis": "External transcriptomic concordance and external-module proximity validation",
            },
            {
                "Resource": "STRING v11.5",
                "Data type": "Protein-protein interaction network",
                "Samples used": "Human proteins; combined score >= 0.4",
                "Platform/source": "STRING physical interaction network",
                "Role in analysis": "Disease module construction and network distance calculations",
            },
            {
                "Resource": "DrugBank/STITCH/Harmonizome/literature",
                "Data type": "Drug-target annotations",
                "Samples used": "RTX, TAC and CTX target sets",
                "Platform/source": "Curated pharmacology databases and targeted literature review",
                "Role in analysis": "Drug target set definition and expansion",
            },
            {
                "Resource": "GSE171458",
                "Data type": "Single-cell RNA-seq",
                "Samples used": "PMN and normal kidney samples from source dataset",
                "Platform/source": "10x Genomics scRNA-seq",
                "Role in analysis": "Cell-type mapping of drug target sets",
            },
            {
                "Resource": "RCSB PDB/PubChem",
                "Data type": "Protein structures and small-molecule structures",
                "Samples used": "Six hub proteins; TAC and CTX ligands",
                "Platform/source": "RCSB Protein Data Bank and PubChem",
                "Role in analysis": "Orthogonal structural plausibility screening by docking",
            },
        ]
    )
    datasets.to_csv(TABLES / "Table1_Datasets_and_resources.csv", index=False)

    disease = pickle.load(open(DATA / "disease_module.pkl", "rb"))
    de = pd.read_csv(RESULTS / "DE_MGN_vs_Normal.csv")
    sig = de[de["adj.P.Val"] < 0.05]
    drug_targets = pickle.load(open(DATA / "drug_targets.pkl", "rb"))

    target_rows = [
        {
            "Set": "PMN DEGs",
            "Layer/source": "Transcriptomic disease seeds",
            "Seed count": len(sig["Gene"].dropna().unique()),
            "Expanded/in-module count": len(disease["deg_genes"]),
            "Representative genes": "TP53, KEAP1, MYDGF, SRSF5",
            "Interpretation": "Active glomerular disease-state signature",
        },
        {
            "Set": "PMN GWAS genes",
            "Layer/source": "Genetic disease seeds",
            "Seed count": len(disease["gwas_genes"]),
            "Expanded/in-module count": len(disease["gwas_genes"]),
            "Representative genes": "PLA2R1, NFKB1, IRF4, HLA-DQA1, HLA-DRB1",
            "Interpretation": "Disease susceptibility core",
        },
        {
            "Set": "PMN LCC module",
            "Layer/source": "STRING-expanded disease module",
            "Seed count": len(disease["seed_genes"]),
            "Expanded/in-module count": len(disease["lcc_genes"]),
            "Representative genes": "NFKB1, TP53, AKT1, RELA, PLA2R1",
            "Interpretation": "Interactome neighborhood used for proximity analysis",
        },
        {
            "Set": "RTX targets",
            "Layer/source": "Direct, Fc receptor, complement, BCR signaling and Treg layers",
            "Seed count": read_gene_count(RESULTS / "RTX_targets.csv"),
            "Expanded/in-module count": len(drug_targets["RTX"]),
            "Representative genes": "MS4A1, FCGR3A, C1QA, SYK, BTK, FOXP3, TGFB1",
            "Interpretation": "Multi-layer biologic perturbation model",
        },
        {
            "Set": "TAC targets",
            "Layer/source": "FKBP/calcineurin/NFAT/T-cell signaling",
            "Seed count": read_gene_count(RESULTS / "TAC_targets.csv"),
            "Expanded/in-module count": len(drug_targets["TAC"]),
            "Representative genes": "FKBP1A, PPP3CA, NFATC1, IL2, JAK1",
            "Interpretation": "T-cell activation neighborhood",
        },
        {
            "Set": "CTX targets",
            "Layer/source": "Metabolism, DNA damage response and apoptosis",
            "Seed count": read_gene_count(RESULTS / "CTX_targets.csv"),
            "Expanded/in-module count": len(drug_targets["CTX"]),
            "Representative genes": "CYP2B6, CYP3A4, TP53, ATM, CASP3, BAX",
            "Interpretation": "Broad cytotoxic/apoptotic neighborhood",
        },
    ]
    pd.DataFrame(target_rows).to_csv(TABLES / "Table2_Target_sets_and_layers.csv", index=False)

    prox = pd.read_csv(RESULTS / "drug_disease_proximity.csv")
    sep = pd.read_csv(RESULTS / "drug_drug_separation.csv")
    metrics = prox.copy()
    metrics = metrics.rename(
        columns={
            "Drug": "Drug or pair",
            "TargetSet": "Disease set",
            "Z_score": "Proximity z-score",
            "Observed": "Observed distance",
            "PermMean": "Permutation mean distance",
        }
    )
    metrics["Disease set"] = metrics["Disease set"].replace({"ALL": "Integrated"})
    metrics["Metric type"] = "Drug-disease proximity"
    metrics["Separation"] = ""
    metrics["d_ab"] = ""
    metrics["d_aa"] = ""
    metrics["d_bb"] = ""
    metrics = metrics[
        [
            "Metric type",
            "Drug or pair",
            "Disease set",
            "Observed distance",
            "Permutation mean distance",
            "Proximity z-score",
            "Separation",
            "d_ab",
            "d_aa",
            "d_bb",
        ]
    ]
    sep_rows = sep.rename(columns={"Pair": "Drug or pair", "Separation": "Separation"}).copy()
    sep_rows["Metric type"] = "Drug-drug separation"
    sep_rows["Disease set"] = ""
    sep_rows["Observed distance"] = ""
    sep_rows["Permutation mean distance"] = ""
    sep_rows["Proximity z-score"] = ""
    sep_rows = sep_rows[
        [
            "Metric type",
            "Drug or pair",
            "Disease set",
            "Observed distance",
            "Permutation mean distance",
            "Proximity z-score",
            "Separation",
            "d_ab",
            "d_aa",
            "d_bb",
        ]
    ]
    pd.concat([metrics, sep_rows], ignore_index=True).to_csv(
        TABLES / "Table3_Network_metrics.csv", index=False
    )


def degree_bins(G, bin_width=5):
    bins = {}
    for node, deg in G.degree():
        key = int(deg // bin_width)
        bins.setdefault(key, []).append(node)
    return bins


def sample_degree_matched(G, targets, bins, bin_width=5):
    sampled = []
    all_nodes = list(G.nodes())
    for t in targets:
        deg = G.degree(t)
        key = int(deg // bin_width)
        candidates = bins.get(key, [])
        if len(candidates) <= 1:
            # Expand search to neighboring bins when a bin is sparse.
            candidates = []
            for delta in range(1, 5):
                candidates.extend(bins.get(key - delta, []))
                candidates.extend(bins.get(key + delta, []))
                if len(candidates) > 1:
                    break
        if not candidates:
            candidates = all_nodes
        sampled.append(random.choice(candidates))
    return sampled


def avg_min_distance(G, sources, targets):
    vals = []
    target_set = set(targets)
    for s in sources:
        lengths = nx.single_source_shortest_path_length(G, s)
        ds = [d for n, d in lengths.items() if n in target_set]
        if ds:
            vals.append(min(ds))
    return float(np.mean(vals)) if vals else math.nan


def robustness_checks():
    random.seed(11)
    np.random.seed(11)
    G = nx.read_graphml(DATA / "mn_ppi_network.graphml").to_undirected()
    lcc = max(nx.connected_components(G), key=len)
    G = G.subgraph(lcc).copy()
    all_nodes = set(G.nodes())

    disease = pickle.load(open(DATA / "disease_module.pkl", "rb"))
    drug_targets = pickle.load(open(DATA / "drug_targets.pkl", "rb"))
    de = pd.read_csv(RESULTS / "DE_MGN_vs_Normal.csv")
    deg_genes = set(de.loc[de["adj.P.Val"] < 0.05, "Gene"].dropna())
    gwas_genes = set(disease["gwas_genes"])
    disease_sets = {
        "DEGs": deg_genes & all_nodes,
        "GWAS": gwas_genes & all_nodes,
        "Integrated": (deg_genes | gwas_genes) & all_nodes,
    }

    bins = degree_bins(G)
    rows = []
    n_perm = 300

    for drug in ["RTX", "TAC", "CTX"]:
        expanded = sorted(set(drug_targets[drug]) & all_nodes)
        seed_file = RESULTS / f"{drug}_targets.csv"
        seeds = sorted(set(pd.read_csv(seed_file)["Gene"].dropna()) & all_nodes)
        variants = [("Expanded targets", expanded), ("Seed targets only", seeds)]
        if drug == "RTX":
            cd20 = [g for g in ["MS4A1"] if g in all_nodes]
            core = [g for g in ["MS4A1", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "C1QA", "C1QB", "C1QC", "C1R", "C1S"] if g in all_nodes]
            variants.extend([("CD20 only", cd20), ("CD20+Fc/complement", core)])

        for disease_label, disease_all in disease_sets.items():
            for label, targets in variants:
                if len(targets) < 2:
                    rows.append(
                        {
                            "Drug": drug,
                            "Disease set": disease_label,
                            "Analysis set": label,
                            "Targets in network": len(targets),
                            "Disease proteins in network": len(disease_all),
                            "Observed distance": "NA",
                            "Degree-matched permutation mean": "NA",
                            "Degree-matched z-score": "NA",
                            "Interpretation": "Not computed: fewer than two targets in network",
                        }
                    )
                    continue
                obs = avg_min_distance(G, targets, disease_all)
                perm = []
                for _ in range(n_perm):
                    sample = sample_degree_matched(G, targets, bins)
                    perm.append(avg_min_distance(G, sample, disease_all))
                std = np.nanstd(perm)
                z = (obs - np.nanmean(perm)) / std if std > 0 else np.nan
                rows.append(
                    {
                        "Drug": drug,
                        "Disease set": disease_label,
                        "Analysis set": label,
                        "Targets in network": len(targets),
                        "Disease proteins in network": len(disease_all),
                        "Observed distance": round(obs, 4),
                        "Degree-matched permutation mean": round(float(np.nanmean(perm)), 4),
                        "Degree-matched z-score": round(float(z), 4) if not np.isnan(z) else "NA",
                        "Interpretation": "Proximal" if not np.isnan(z) and z < -1.5 else "Not proximal by z < -1.5",
                    }
                )

    pd.DataFrame(rows).to_csv(TABLES / "Supplementary_Table_S3_Robustness_checks.csv", index=False)


def load_font(size=42):
    for path in [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf",
    ]:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def panel_label(img, label):
    img = img.convert("RGB")
    draw = ImageDraw.Draw(img)
    font = load_font(max(36, img.width // 35))
    pad = max(12, img.width // 120)
    box = draw.textbbox((0, 0), label, font=font)
    draw.rectangle((pad, pad, pad + box[2] + 20, pad + box[3] + 16), fill="white")
    draw.text((pad + 10, pad + 8), label, fill="black", font=font)
    return img


def resize_to_width(img, width):
    h = int(img.height * width / img.width)
    return img.resize((width, h), Image.LANCZOS)


def merge_grid(images, out_path, width=1800, gutter=60, bg="white"):
    imgs = [resize_to_width(img, width) for img in images]
    rows = []
    for i in range(0, len(imgs), 2):
        rows.append(imgs[i : i + 2])
    row_widths = [sum(im.width for im in row) + gutter * (len(row) - 1) for row in rows]
    row_heights = [max(im.height for im in row) for row in rows]
    canvas = Image.new("RGB", (max(row_widths), sum(row_heights) + gutter * (len(rows) - 1)), bg)
    y = 0
    for row, rh, rw in zip(rows, row_heights, row_widths):
        x = (canvas.width - rw) // 2
        for im in row:
            canvas.paste(im, (x, y))
            x += im.width + gutter
        y += rh + gutter
    canvas.save(out_path, dpi=(300, 300))


def merge_figures():
    fig1 = panel_label(Image.open(FIGURES / "Figure1_Workflow.png"), "A")
    fig2 = panel_label(Image.open(FIGURES / "Figure5_PPI_Network.png"), "B")
    merge_grid([fig1, fig2], SUBFIG / "Figure1_Workflow_DiseaseModule.png", width=1800)

    f2a = panel_label(Image.open(FIGURES / "Figure2_Drug_Disease_Proximity.png"), "A")
    f2b = panel_label(Image.open(FIGURES / "Figure3_Drug_Drug_Separation.png"), "B")
    f2c = panel_label(Image.open(FIGURES / "Figure4_Complementary_Exposure.png"), "C")
    merge_grid([f2a, f2b, f2c], SUBFIG / "Figure2_NetworkMetrics.png", width=1700)

    f3a = panel_label(Image.open(FIGURES / "Supplementary_Figure_CellType_Mapping.png"), "A")
    f3b = panel_label(Image.open(FIGURES / "Figure6_Summary.png"), "B")
    merge_grid([f3a, f3b], SUBFIG / "Figure3_CellType_Mechanism.png", width=1800)

    f4 = panel_label(Image.open(FIGURES / "Figure7_Molecular_Docking.png"), "S1")
    f4.save(SUBFIG / "Supplementary_Figure_S1_Docking.png", dpi=(300, 300))


def main():
    write_tables()
    robustness_checks()
    merge_figures()
    print("Generated reviewer revision assets")


if __name__ == "__main__":
    main()
