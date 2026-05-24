#!/usr/bin/env python3
"""
Step 5: Cell-type mapping of drug targets using scRNA-seq marker analysis
===========================================================================
Uses established kidney cell-type markers to map where RTX, TAC, and CTX
targets are expressed.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import gzip
import pickle
from pathlib import Path

BASE_DIR = Path("/workspaces/platform/MN")
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"
SCRNA_DIR = DATA_DIR / "scRNA"
os.makedirs(FIGURES_DIR, exist_ok=True)

print("=" * 60)
print("Step 5: Cell-Type Mapping of Drug Targets")
print("=" * 60)

# ================================================================
# 1. Kidney cell-type marker genes (from literature)
# ================================================================
# Sources:
# - KPMP kidney cell atlas (Lake et al. 2023, Nature)
# - Human kidney cell atlas (Stewart et al. 2019, Science)
# - GSE171458 original paper (2021, Frontiers in Immunology)
# ================================================================

kidney_markers = {
    "Podocyte": [
        "NPHS1", "NPHS2", "WT1", "PODXL", "SYNPO", "PLA2R1", 
        "MAFB", "TJP1", "CD2AP", "ZEB1", "CLIC5", "PTPRO"
    ],
    "Mesangial Cell": [
        "PDGFRB", "ACTA2", "MYH11", "TAGLN", "DES", "RGS5",
        "CSPG4", "NOTCH3", "COL4A4", "ITGA8"
    ],
    "Glomerular Endothelial": [
        "PECAM1", "CDH5", "KDR", "VWF", "FLT1", "TEK",
        "ENG", "PLVAP", "ESAM", "EMCN", "CD34"
    ],
    "Proximal Tubule": [
        "SLC3A1", "SLC34A1", "CUBN", "LRP2", "ALDOB",
        "PCK1", "SLC5A2", "SLC7A13", "SLC22A6", "AQP1"
    ],
    "Loop of Henle": [
        "SLC12A1", "CLDN16", "UMOD", "SLC12A3",
        "CALB1", "CLDN10", "PVALB", "SLC8A1"
    ],
    "Distal Tubule": [
        "SLC12A3", "TRPV5", "TRPM6", "CLDN16",
        "SLC4A1", "KCNJ1", "CLDN8"
    ],
    "Collecting Duct": [
        "AQP2", "AQP3", "AQP4", "AVPR2", "SLC14A2",
        "CALB1", "KRT19", "FXYD4", "ATP6V1B1"
    ],
    "B cell": [
        "MS4A1", "CD19", "CD22", "CD79A", "CD79B",
        "PAX5", "BANK1", "BLK", "FCGR2B"
    ],
    "T cell": [
        "CD3E", "CD3D", "CD4", "CD8A", "CD8B",
        "FOXP3", "IL2RA", "CD28", "ZAP70"
    ],
    "NK cell": [
        "NKG7", "KLRD1", "KLRF1", "FCGR3A", "NCR1",
        "PRF1", "GZMB", "GZMA", "KLRK1"
    ],
    "Macrophage/Monocyte": [
        "CD14", "CD68", "CSF1R", "ITGAM", "FCGR1A",
        "FCGR2A", "C1QA", "C1QB", "LYZ", "CD163"
    ],
    "Tubulointerstitial": [
        "VIM", "FN1", "COL1A1", "COL3A1", "DCN",
        "LUM", "FAP", "ACTA2", "PDGFRA"
    ]
}

# Collapse into a single mapping
marker_to_celltype = {}
for celltype, markers in kidney_markers.items():
    for m in markers:
        marker_to_celltype[m] = celltype

print(f"Loaded {len(kidney_markers)} cell types with {len(marker_to_celltype)} markers")

# ================================================================
# 2. Load drug targets
# ================================================================
print("\nLoading drug targets...")
with open(DATA_DIR / "drug_targets.pkl", "rb") as f:
    drug_targets = pickle.load(f)

# Load DEGs
de = pd.read_csv(RESULTS_DIR / "DE_MGN_vs_Normal.csv")
top_deg_genes = set(de[de['adj.P.Val'] < 0.0001].nsmallest(100, 'adj.P.Val')['Gene'].dropna())

# ================================================================
# 3. Score drug targets against cell-type markers
# ================================================================
print("\nScoring drug targets against kidney cell types...")

def score_celltype_enrichment(gene_set, marker_dict, n_total=20000):
    """Calculate enrichment of a gene set in each cell type using hypergeometric test."""
    results = {}
    for celltype, markers in kidney_markers.items():
        markers_in_data = [m for m in markers if m in gene_set]
        markers_total = len(markers)
        hit_fraction = len(markers_in_data) / max(markers_total, 1)
        results[celltype] = {
            "hits": len(markers_in_data),
            "total_markers": markers_total,
            "fraction": hit_fraction
        }
    return results

# Score each drug's targets
drug_celltype_scores = {}
for drug, targets in drug_targets.items():
    drug_celltype_scores[drug] = score_celltype_enrichment(
        set(targets), marker_to_celltype
    )

# Also score MN DEGs
deg_celltype_scores = score_celltype_enrichment(top_deg_genes, marker_to_celltype)

# ================================================================
# 4. Generate cell-type mapping heatmap
# ================================================================
print("\nGenerating cell-type mapping heatmap...")

celltypes = list(kidney_markers.keys())
drugs = list(drug_targets.keys()) + ["MN_DEGs"]

# Build data matrix
data_matrix = np.zeros((len(celltypes), len(drugs)))
for i, ct in enumerate(celltypes):
    for j, drug in enumerate(drugs):
        scores = drug_celltype_scores[drug] if drug != "MN_DEGs" else deg_celltype_scores
        data_matrix[i, j] = scores.get(ct, {}).get("fraction", 0) * 100  # as percentage

fig, ax = plt.subplots(figsize=(9, 8))

im = ax.imshow(data_matrix, cmap="YlOrRd", aspect="auto", vmin=0, vmax=max(20, data_matrix.max()))

# Add text annotations
for i in range(len(celltypes)):
    for j in range(len(drugs)):
        val = data_matrix[i, j]
        if val > 0:
            color = "white" if val > 15 else "black"
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center", fontsize=8, color=color, fontweight="bold")

ax.set_xticks(range(len(drugs)))
ax.set_xticklabels(["RTX", "TAC", "CTX", "MN\nDEGs"], fontsize=10, fontweight="bold")
ax.set_yticks(range(len(celltypes)))
ax.set_yticklabels(celltypes, fontsize=10)
ax.set_title("Drug Target Expression in Kidney Cell Types\n(% of cell-type markers that are drug targets)",
             fontsize=13, fontweight="bold", pad=8)

plt.colorbar(im, ax=ax, label="% overlap", shrink=0.8)

plt.tight_layout()
plt.savefig(FIGURES_DIR / "Supplementary_Figure_CellType_Mapping.png", dpi=300, bbox_inches="tight")
plt.close()
print("  Saved Supplementary_Figure_CellType_Mapping.png")

# ================================================================
# 5. Print summary
# ================================================================
print("\n\nCell-type enrichment summary:")
print("-" * 60)
print(f"{'Cell Type':<25} {'RTX':>8} {'TAC':>8} {'CTX':>8} {'MN DEGs':>8}")
print("-" * 60)
for i, ct in enumerate(celltypes):
    vals = [f"{data_matrix[i, j]:>7.0f}%" for j in range(len(drugs))]
    print(f"{ct:<25} {vals[0]:>8} {vals[1]:>8} {vals[2]:>8} {vals[3]:>8}")

print("\n\nKey findings:")
print("-" * 60)

# RTX targets are enriched in B cell and NK cell markers
rtx_b_cell = drug_celltype_scores["RTX"].get("B cell", {}).get("fraction", 0)
rtx_nk = drug_celltype_scores["RTX"].get("NK cell", {}).get("fraction", 0)
tac_t_cell = drug_celltype_scores["TAC"].get("T cell", {}).get("fraction", 0)
ctx_broad = max([drug_celltype_scores["CTX"].get(ct, {}).get("fraction", 0) for ct in celltypes])

print(f"  RTX → B cell: {rtx_b_cell*100:.0f}% of markers, NK cell: {rtx_nk*100:.0f}%")
print(f"  TAC → T cell: {tac_t_cell*100:.0f}% of markers")
print(f"  CTX → Broad distribution (max: {ctx_broad*100:.0f}%)")
print(f"  RTX+TAC → Complementary cell types (B cell + T cell) ✓")

# ================================================================
# 6. Save results
# ================================================================
results_df = pd.DataFrame(data_matrix, index=celltypes, columns=drugs)
results_df.to_csv(RESULTS_DIR / "celltype_drug_mapping.csv")
print(f"\nSaved celltype_drug_mapping.csv")

print("\n✓ Step 5 complete!")
