#!/usr/bin/env python3
"""Methodology strengthening analyses for the MN network medicine manuscript.

Generates reviewer-facing supplementary analyses:
- corrected primary proximity/separation after fixing PLCG2 spelling
- RTX single-node/layer ablation
- STRING edge-threshold sensitivity
- comparator/negative-control target sets
- target curation evidence table
"""
from __future__ import annotations

import math
import pickle
import random
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data"
RESULTS = BASE / "results"
TABLES = BASE / "tables"
TABLES.mkdir(exist_ok=True)
SEED = 20260519
GWAS = {"PLA2R1", "NFKB1", "IRF4", "HLA-DQA1", "HLA-DRB1"}

RTX_LAYERS = {
    "CD20 only": {"MS4A1"},
    "CD20 + Fc/complement": {"MS4A1", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B", "C1QA", "C1QB", "C1QC", "C1R", "C1S"},
    "BCR/signaling layer only": {"SYK", "LYN", "BTK", "PLCG2", "PIK3CA", "PIK3CB", "PIK3CD", "AKT1", "AKT2", "MTOR", "MAPK1", "MAPK3", "MAPK8", "NFKB1", "RELA", "IKBKB", "IKBKG", "CHUK"},
    "Immune regulation layer only": {"TGFB1", "TGFBR1", "TGFBR2", "FOXP3", "IL2RA", "CD4", "PRF1", "GZMB", "IFNG", "TNFSF10", "CD19", "CD22", "CD79A", "CD79B", "BLK", "BANK1", "IL10", "IL6", "TNF", "IL2"},
}
RTX_SEEDS = set().union(*RTX_LAYERS.values())
TAC_SEEDS = {"FKBP1A", "FKBP1B", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", "NFATC1", "NFATC2", "NFATC3", "NFATC4", "IL2", "IL2RA", "IL2RB", "IL2RG", "JAK1", "JAK3", "STAT5A", "STAT5B", "MYC", "CCND1", "CDK4"}
CTX_SEEDS = {"ALDH1A1", "CYP2B6", "CYP3A4", "CYP2C9", "TP53", "CDKN1A", "BAX", "BCL2", "CASP3", "CASP8", "CASP9", "ATM", "ATR", "CHEK1", "CHEK2", "FAS", "FASLG", "GSTA1", "GSTP1", "MGMT"}
SEEDS = {"RTX": RTX_SEEDS, "TAC": TAC_SEEDS, "CTX": CTX_SEEDS}
COMPARATORS = {
    "Cyclosporine calcineurin control": {"PPIA", "PPIB", "PPIF", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", "NFATC1", "NFATC2", "NFATC3", "NFATC4", "IL2"},
    "EGFR/ERBB oncology control": {"EGFR", "ERBB2", "ERBB3", "ERBB4", "GRB2", "SOS1", "KRAS", "NRAS", "HRAS", "RAF1", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "PIK3CA", "AKT1"},
    "BCR-ABL oncology control": {"ABL1", "BCR", "KIT", "PDGFRA", "PDGFRB", "SRC", "LYN", "STAT5A", "STAT5B", "CRKL"},
    "Metabolic/lipid negative control": {"ACACA", "FASN", "CPT1A", "SCD", "SREBF1", "PPARA", "PPARG", "APOE", "APOB", "LDLR"},
    "Ion-channel negative control": {"SCN5A", "KCNH2", "KCNQ1", "CACNA1C", "RYR2", "ATP1A1", "ATP2A2", "KCNA5", "KCNJ2", "CACNB2"},
}


def graph_lcc(path=DATA / "mn_ppi_network.graphml"):
    g = nx.read_graphml(path).to_undirected()
    for _, _, d in g.edges(data=True):
        d["weight"] = float(d.get("weight", 1.0))
    lcc = max(nx.connected_components(g), key=len)
    return g.subgraph(lcc).copy()


def disease_sets(g):
    nodes = set(g.nodes())
    de = pd.read_csv(RESULTS / "DE_MGN_vs_Normal.csv")
    deg = set(de.loc[de["adj.P.Val"] < 0.05, "Gene"].dropna()) & nodes
    gwas = GWAS & nodes
    return {"DEGs": deg, "GWAS": gwas, "Integrated": deg | gwas}


def expand(g, genes, threshold=0.4, max_neighbors=10):
    out = set(genes)
    for x in genes:
        if x not in g:
            continue
        nbrs = sorted(((n, float(g[x][n].get("weight", 0))) for n in g.neighbors(x)), key=lambda z: -z[1])
        out.update(n for n, w in nbrs[:max_neighbors] if w >= threshold)
    return out


def avg_min(g, src, dst):
    src = [x for x in set(src) if x in g]
    dst = set(x for x in dst if x in g)
    vals = []
    for x in src:
        lengths = nx.single_source_shortest_path_length(g, x)
        ds = [d for n, d in lengths.items() if n in dst]
        if ds:
            vals.append(min(ds))
    return float(np.mean(vals)) if vals else math.nan


def nearest_other(g, src):
    src = [x for x in set(src) if x in g]
    vals = []
    for x in src:
        lengths = nx.single_source_shortest_path_length(g, x)
        ds = [d for n, d in lengths.items() if n in src and n != x]
        if ds:
            vals.append(min(ds))
    return float(np.mean(vals)) if vals else math.nan


def prox(g, src, dst, n_perm=300, seed=SEED):
    src = sorted(set(src) & set(g.nodes()))
    dst = sorted(set(dst) & set(g.nodes()))
    if len(src) < 2 or len(dst) < 2:
        return None
    rng = random.Random(seed)
    nodes = list(g.nodes())
    # Precompute distance from every node to the disease set once. This is much
    # faster than running single-source shortest paths for every permutation.
    disease_distance = nx.multi_source_dijkstra_path_length(g, dst, weight=None)
    def mean_distance(source_nodes):
        vals = [disease_distance[x] for x in source_nodes if x in disease_distance]
        return float(np.mean(vals)) if vals else math.nan
    obs = mean_distance(src)
    perm = [mean_distance(rng.sample(nodes, len(src))) for _ in range(n_perm)]
    sd = float(np.nanstd(perm))
    z = float((obs - np.nanmean(perm)) / sd) if sd > 0 else math.nan
    return len(src), len(dst), obs, float(np.nanmean(perm)), z


def sep(g, a, b):
    a = sorted(set(a) & set(g.nodes()))
    b = sorted(set(b) & set(g.nodes()))
    if len(a) < 2 or len(b) < 2:
        return None
    daa = nearest_other(g, a)
    dbb = nearest_other(g, b)
    dab = avg_min(g, a, b)
    return dab - (daa + dbb) / 2, dab, daa, dbb


def filter_graph(g, threshold):
    h = nx.Graph()
    for u, v, d in g.edges(data=True):
        if float(d.get("weight", 0)) >= threshold:
            h.add_edge(u, v, **d)
    if h.number_of_nodes() == 0:
        return h
    lcc = max(nx.connected_components(h), key=len)
    return h.subgraph(lcc).copy()


def result_row(label, disease, res, seed_count=None, note_prefix=""):
    if res is None:
        return {**label, "Seed targets": seed_count, "Targets in network": 0, "Disease proteins in network": len(disease), "Observed distance": "NA", "Permutation mean distance": "NA", "Proximity z-score": "NA", "Interpretation": note_prefix + "Not computed: fewer than two targets or disease proteins in network"}
    nsrc, ndst, obs, pm, z = res
    interp = "Proximal" if z < -1.5 else "Not proximal by z < -1.5"
    return {**label, "Seed targets": seed_count, "Targets in network": nsrc, "Disease proteins in network": ndst, "Observed distance": round(obs, 4), "Permutation mean distance": round(pm, 4), "Proximity z-score": round(z, 4), "Interpretation": note_prefix + interp}


def write_target_evidence(g, expanded):
    rows = []
    for drug, genes in SEEDS.items():
        for gene in sorted(genes):
            if drug == "RTX":
                layer = next(k for k, v in RTX_LAYERS.items() if gene in v)
                source = "DrugBank/literature" if layer == "CD20 only" else "Literature/Harmonizome"
            elif drug == "TAC":
                layer, source = "FKBP/calcineurin/NFAT or T-cell signaling", "DrugBank/STITCH/literature"
            else:
                layer, source = "Metabolism, DNA damage or apoptosis", "DrugBank/STITCH/literature"
            rows.append({"Drug": drug, "Target": gene, "Layer": layer, "Evidence source category": source, "In final PPI LCC": gene in g, "In expanded analysis set": gene in expanded[drug], "Comment": "Curated target/pathway component; not all targets are present in the STRING LCC"})
    pd.DataFrame(rows).to_csv(TABLES / "Supplementary_Table_S7_Target_curation_evidence.csv", index=False)


def main():
    g = graph_lcc()
    ds = disease_sets(g)
    expanded = {d: expand(g, s) for d, s in SEEDS.items()}

    # Persist corrected targets and expanded target sets.
    for d, s in SEEDS.items():
        pd.DataFrame({"Gene": sorted(s)}).to_csv(RESULTS / f"{d}_targets.csv", index=False)
    with open(DATA / "drug_targets.pkl", "wb") as f:
        pickle.dump(expanded, f)

    # Corrected primary result files.
    prox_rows = []
    for d, s in expanded.items():
        for name, disease in [("DEGs", ds["DEGs"]), ("GWAS", ds["GWAS"]), ("ALL", ds["Integrated"] )]:
            r = prox(g, s, disease, n_perm=1000, seed=SEED + len(prox_rows))
            if r:
                _, _, obs, pm, z = r
                prox_rows.append({"Drug": d, "TargetSet": name, "Z_score": z, "Observed": obs, "PermMean": pm})
    pd.DataFrame(prox_rows).to_csv(RESULTS / "drug_disease_proximity.csv", index=False)
    sep_rows = []
    for a, b in [("RTX", "TAC"), ("RTX", "CTX"), ("TAC", "CTX")]:
        r = sep(g, expanded[a], expanded[b])
        if r:
            s, dab, daa, dbb = r
            sep_rows.append({"Pair": f"{a}_{b}", "Drug1": a, "Drug2": b, "Separation": s, "d_ab": dab, "d_aa": daa, "d_bb": dbb})
    pd.DataFrame(sep_rows).to_csv(RESULTS / "drug_drug_separation.csv", index=False)

    # S4: RTX ablation.
    rows = []
    ablations = dict(RTX_LAYERS)
    ablations["Full seed set"] = RTX_SEEDS
    ablations["Full expanded set"] = expanded["RTX"]
    for rep, genes in ablations.items():
        for name, disease in ds.items():
            r = prox(g, genes, disease, n_perm=300, seed=SEED + len(rows))
            rows.append(result_row({"RTX representation": rep, "Disease set": name}, disease, r, seed_count=len(genes)))
    pd.DataFrame(rows).to_csv(TABLES / "Supplementary_Table_S4_RTX_layer_ablation.csv", index=False)

    # S5: threshold sensitivity.
    rows = []
    for t in [0.4, 0.7, 0.9]:
        h = filter_graph(g, t)
        hds = disease_sets(h) if h.number_of_nodes() else {}
        for d, seeds in SEEDS.items():
            hs = expand(h, seeds, threshold=t) if h.number_of_nodes() else set()
            for name, disease in hds.items():
                r = prox(h, hs, disease, n_perm=300, seed=SEED + len(rows))
                label = {"STRING edge threshold": t, "Drug": d, "Disease set": name, "Network nodes": h.number_of_nodes(), "Network edges": h.number_of_edges()}
                rows.append(result_row(label, disease, r, seed_count=len(seeds)))
    pd.DataFrame(rows).to_csv(TABLES / "Supplementary_Table_S5_STRING_threshold_sensitivity.csv", index=False)

    # S6: comparators and negative controls.
    rows = []
    for name, genes in COMPARATORS.items():
        r = prox(g, genes, ds["Integrated"], n_perm=1000, seed=SEED + len(rows))
        kind = "mechanistic comparator" if "control" in name and "negative" not in name else "negative control"
        rows.append(result_row({"Comparator set": name, "Control type": kind, "Disease set": "Integrated"}, ds["Integrated"], r, seed_count=len(genes)))
    pd.DataFrame(rows).to_csv(TABLES / "Supplementary_Table_S6_Comparator_controls.csv", index=False)

    write_target_evidence(g, expanded)
    print("Generated strengthened methodology tables and corrected primary results")


if __name__ == "__main__":
    main()
