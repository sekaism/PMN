#!/usr/bin/env python3
"""
Step 4: Generate paper figures
================================
Creates all figures for the manuscript.
"""

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import networkx as nx
import numpy as np
import pandas as pd
import pickle
import os
from pathlib import Path
import textwrap

BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "results"
DATA_DIR = BASE_DIR / "data"
FIGURES_DIR = BASE_DIR / "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# Set publication-quality style
plt.rcParams.update({
    'figure.dpi': 300,
    'font.family': 'sans-serif',
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'legend.fontsize': 10,
    'figure.figsize': (8, 6)
})

# ================================================================
# Figure 1: Study Overview / Workflow
# ================================================================
print("Generating Figure 1: Workflow...")
fig, ax = plt.subplots(figsize=(10, 7))
ax.axis('off')

# Create a workflow diagram using text boxes
steps = [
    ("Data Collection", 
     "• GSE108113: 87 MGN + 11 Normal\n• MN GWAS (Xie et al. 2020)\n• Drug targets (DrugBank)\n• STRING PPI network"),
    ("Transcriptomic\nDE Analysis",
     "• 3,063 DEGs (adj.P<0.05)\n• 1,601 up in MN\n• 1,462 down in MN"),
    ("Disease Module\nConstruction",
     "• DEGs + GWAS genes\n• STRING PPI expansion\n• 2,176 nodes, 3,098 edges"),
    ("Network Medicine\nAnalysis",
     "• Drug-Disease PROXIMITY\n• Drug-Drug SEPARATION\n• Complementary Exposure"),
    ("Validation &\nInterpretation",
     "• scRNA-seq cell-type mapping\n• Enrichment analysis\n• Clinical trial comparison")
]

# Position boxes
y_positions = [0.85, 0.67, 0.49, 0.31, 0.13]
colors = ['#4472C4', '#ED7D31', '#70AD47', '#FFC000', '#5B9BD5']

for i, (title, content) in enumerate(steps):
    y = y_positions[i]
    # Box
    ax.add_patch(plt.Rectangle((0.05, y - 0.09), 0.25, 0.18, 
                                facecolor=colors[i], alpha=0.8, edgecolor='black', linewidth=1))
    ax.text(0.175, y + 0.03, title, ha='center', va='center', fontsize=11, fontweight='bold', color='white')
    
    # Content
    ax.text(0.35, y, content, ha='left', va='center', fontsize=9, fontfamily='monospace')
    
    # Arrow
    if i < len(steps) - 1:
        mid_y = (y_positions[i] + y_positions[i+1]) / 2
        ax.annotate('', xy=(0.175, y_positions[i+1] + 0.09), xytext=(0.175, y - 0.09),
                   arrowprops=dict(arrowstyle='->', lw=2, color='gray'))

fig.suptitle('Study Workflow: Multi-Layer Network Medicine Framework for\n'
             'Rituximab in Membranous Nephropathy', 
             fontsize=15, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig(FIGURES_DIR / "Figure1_Workflow.png", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved Figure1_Workflow.png")

# ================================================================
# Figure 2: Drug-Disease Proximity (bar chart)
# ================================================================
print("Generating Figure 2: Drug-Disease Proximity...")

prox = pd.read_csv(RESULTS_DIR / "drug_disease_proximity.csv")
prox_all = prox[prox['TargetSet'] == 'ALL'].copy()

fig, ax = plt.subplots(figsize=(8, 5))

drugs = prox_all['Drug'].tolist()
z_scores = prox_all['Z_score'].tolist()
colors2 = ['#4472C4' if z < -1.5 else '#FFC000' if z < 0 else '#C00000' 
           for z in z_scores]

bars = ax.bar(drugs, z_scores, color=colors2, edgecolor='black', linewidth=1.2, width=0.5)

# Add value labels
for bar, z in zip(bars, z_scores):
    ax.text(bar.get_x() + bar.get_width()/2, 
            bar.get_height() - 0.3 if z < 0 else bar.get_height() + 0.3,
            f'z = {z:.2f}', ha='center', va='bottom' if z < 0 else 'bottom',
            fontweight='bold', fontsize=12)

# Threshold line
ax.axhline(y=-1.5, color='red', linestyle='--', alpha=0.5, label='Significance threshold (z=-1.5)')
ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)

ax.set_ylabel('Proximity Z-score (more negative = closer to disease)', fontsize=12)
ax.set_xlabel('Drug', fontsize=12)
ax.set_title('Drug-Disease Network Proximity to MN Disease Module', fontsize=14, fontweight='bold')
ax.legend(loc='lower right')
ax.set_ylim(min(z_scores) - 1, 1)

# Add annotation
ax.text(0.5, -0.30, 'More proximal to MN disease module →', 
        transform=ax.transAxes, ha='center', fontsize=9, color='#4472C4', style='italic')

plt.tight_layout()
plt.savefig(FIGURES_DIR / "Figure2_Drug_Disease_Proximity.png", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved Figure2_Drug_Disease_Proximity.png")

# ================================================================
# Figure 3: Drug-Drug Separation (bar chart)
# ================================================================
print("Generating Figure 3: Drug-Drug Separation...")

sep = pd.read_csv(RESULTS_DIR / "drug_drug_separation.csv")

fig, ax = plt.subplots(figsize=(8, 5))

pairs = sep['Pair'].tolist()
scores = sep['Separation'].tolist()
# Color based on separation interpretation
colors3 = ['#70AD47' if s > 0.5 else '#FFC000' for s in scores]

bars = ax.bar(pairs, scores, color=colors3, edgecolor='black', linewidth=1.2, width=0.5)

for bar, s in zip(bars, scores):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
            f's = {s:.2f}', ha='center', fontweight='bold', fontsize=12)

# Threshold
ax.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
ax.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5, label='Separation threshold (s=0.5)')

ax.set_ylabel('Separation Score (positive = separated targets)', fontsize=12)
ax.set_xlabel('Drug Pair', fontsize=12)
ax.set_title('Drug-Drug Target Network Separation', fontsize=14, fontweight='bold')
ax.legend(loc='lower right')
ax.set_ylim(min(scores) - 0.3, max(scores) + 0.5)

# Add interpretation note
ax.text(0.5, -0.35, 'Positive separation → targets in different network neighborhoods → Complementary Exposure',
        transform=ax.transAxes, ha='center', fontsize=9, style='italic')

plt.tight_layout()
plt.savefig(FIGURES_DIR / "Figure3_Drug_Drug_Separation.png", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved Figure3_Drug_Drug_Separation.png")

# ================================================================
# Figure 4: Complementary Exposure Network Diagram
# ================================================================
print("Generating Figure 4: Complementary Exposure Diagram...")

fig, ax = plt.subplots(figsize=(10, 8))
ax.axis('off')

# Draw a visual representation of the Complementary Exposure concept

# Disease module (large circle)
circle_disease = plt.Circle((0.3, 0.5), 0.25, color='#4472C4', alpha=0.2, 
                            ec='#4472C4', linewidth=2)
ax.add_patch(circle_disease)
ax.text(0.3, 0.5, 'MN Disease\nModule', ha='center', va='center', 
        fontsize=13, fontweight='bold', color='#4472C4')

# RTX target module
circle_rtx = plt.Circle((0.15, 0.35), 0.12, color='#ED7D31', alpha=0.4, 
                        ec='#ED7D31', linewidth=2)
ax.add_patch(circle_rtx)
ax.text(0.15, 0.35, 'RTX\nTargets', ha='center', va='center', 
        fontsize=11, fontweight='bold', color='white')

# TAC target module  
circle_tac = plt.Circle((0.5, 0.35), 0.12, color='#70AD47', alpha=0.4, 
                        ec='#70AD47', linewidth=2)
ax.add_patch(circle_tac)
ax.text(0.5, 0.35, 'TAC\nTargets', ha='center', va='center', 
        fontsize=11, fontweight='bold', color='white')

# CTX target module
circle_ctx = plt.Circle((0.8, 0.35), 0.12, color='#C00000', alpha=0.4,
                        ec='#C00000', linewidth=2)
ax.add_patch(circle_ctx)
ax.text(0.8, 0.35, 'CTX\nTargets', ha='center', va='center',
        fontsize=11, fontweight='bold', color='white')

# Complementary Exposure arrows for RTX+TAC
# Both RTX and TAC overlap with disease module
ax.annotate('', xy=(0.3, 0.5), xytext=(0.15, 0.35),
           arrowprops=dict(arrowstyle='->', color='#ED7D31', lw=2, connectionstyle='arc3,rad=0.3'))
ax.annotate('', xy=(0.3, 0.5), xytext=(0.5, 0.35),
           arrowprops=dict(arrowstyle='->', color='#70AD47', lw=2, connectionstyle='arc3,rad=-0.3'))

# Separation between RTX and TAC
ax.annotate('', xy=(0.27, 0.33), xytext=(0.38, 0.33),
           arrowprops=dict(arrowstyle='<->', color='gray', lw=2, linestyle='--'))
ax.text(0.325, 0.30, 's = 1.14\n(SEPARATED)', ha='center', va='top', 
        fontsize=9, color='gray')

# Labels
ax.text(0.3, 0.85, 'Complementary Exposure Pattern', 
        ha='center', fontsize=15, fontweight='bold')
ax.text(0.3, 0.78, 'Both drugs target the disease module in separate neighborhoods',
        ha='center', fontsize=11, style='italic')

# Legend
ax.text(0.7, 0.85, 'Proximity to MN disease module:', fontsize=10, fontweight='bold')
ax.text(0.7, 0.78, 'RTX: z = -5.51  ✓ Primary-model proximal', fontsize=9, color='#ED7D31')
ax.text(0.7, 0.72, 'TAC: z = -1.85  ✓ Primary-model proximal', fontsize=9, color='#70AD47')
ax.text(0.7, 0.66, 'CTX: z = -2.61  ✓ Primary-model proximal', fontsize=9, color='#C00000')

# Clinical implication box
ax.add_patch(FancyBboxPatch((0.05, 0.02), 0.9, 0.10, 
                             boxstyle="round,pad=0.05",
                             facecolor='#FFF2CC', alpha=0.8))
ax.text(0.5, 0.07, 
        'Interpretation: RTX + TAC shows primary-model Complementary Exposure (s = 1.14)\n'
        'This is a testable topological hypothesis; TAC proximity is not robust under stricter sensitivity analyses.',
        ha='center', va='center', fontsize=10, style='italic')

plt.tight_layout()
plt.savefig(FIGURES_DIR / "Figure4_Complementary_Exposure.png", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved Figure4_Complementary_Exposure.png")

# ================================================================
# Figure 5: PPI Network Visualization (disease module with drug targets)
# ================================================================
print("Generating Figure 5: PPI Network Visualization...")

G = nx.read_graphml(DATA_DIR / "mn_ppi_network.graphml")
G = G.to_undirected()
components = list(nx.connected_components(G))
lcc = max(components, key=len)
G = G.subgraph(lcc).copy()

# Load drug targets
with open(DATA_DIR / "drug_targets.pkl", "rb") as f:
    drug_targets = pickle.load(f)

# Load DEGs
de = pd.read_csv(RESULTS_DIR / "DE_MGN_vs_Normal.csv")
sig_deg = set(de[de['adj.P.Val'] < 0.0001].nsmallest(50, 'adj.P.Val')['Gene'].dropna())

# Also get GWAS genes
gwas_genes = {"PLA2R1", "NFKB1", "IRF4", "HLA-DQA1", "HLA-DRB1"}

# Select a subnetwork for visualization: drug targets + top DEGs + their 1-step neighbors
selected_nodes = set()
for drug in drug_targets:
    selected_nodes.update([n for n in drug_targets[drug] if n in G])
selected_nodes.update([n for n in sig_deg if n in G])
selected_nodes.update([n for n in gwas_genes if n in G])

# Add 1-step neighbors (to show network context)
all_neighbors = set()
for node in list(selected_nodes):
    if node in G:
        all_neighbors.update(list(G.neighbors(node)))
selected_nodes.update(list(all_neighbors)[:50])  # Limit for visualization

# Extract subgraph
subG = G.subgraph(selected_nodes).copy()

# Keep only LCC of subgraph
sub_comp = list(nx.connected_components(subG))
sub_lcc = max(sub_comp, key=len) if sub_comp else set()
subG = subG.subgraph(sub_lcc).copy()

print(f"  Subnetwork: {subG.number_of_nodes()} nodes, {subG.number_of_edges()} edges")

# Node colors
node_colors = []
node_sizes = []
for node in subG.nodes():
    if node in drug_targets.get("RTX", set()):
        node_colors.append('#ED7D31')  # Orange for RTX
        node_sizes.append(120)
    elif node in drug_targets.get("TAC", set()):
        node_colors.append('#70AD47')  # Green for TAC
        node_sizes.append(100)
    elif node in drug_targets.get("CTX", set()):
        node_colors.append('#C00000')  # Red for CTX
        node_sizes.append(100)
    elif node in gwas_genes:
        node_colors.append('#FF0000')  # Bright red for GWAS
        node_sizes.append(150)
    elif node in sig_deg:
        node_colors.append('#4472C4')  # Blue for DEGs
        node_sizes.append(60)
    else:
        node_colors.append('#D3D3D3')  # Gray for others
        node_sizes.append(30)

# Layout
pos = nx.spring_layout(subG, k=0.3, iterations=50, seed=42)

fig, ax = plt.subplots(figsize=(14, 12))

nx.draw_networkx_edges(subG, pos, alpha=0.15, edge_color='gray', ax=ax)
nx.draw_networkx_nodes(subG, pos, node_color=node_colors, node_size=node_sizes, 
                        edgecolors='black', linewidths=0.5, ax=ax)

# Label key nodes
label_nodes = {}
for node in subG.nodes():
    if node in drug_targets.get("RTX", set()) and node in ["MS4A1", "SYK", "LYN", "BTK", "FCGR3A", "AKT1"]:
        label_nodes[node] = node
    elif node in gwas_genes:
        label_nodes[node] = node
    elif node in drug_targets.get("TAC", set()) and node in ["FKBP1A", "PPP3CA", "NFATC1", "IL2"]:
        label_nodes[node] = node
    elif node in ["CD19", "CD4", "FOXP3", "PLCG2", "NFKB1", "TNF", "IL6"]:
        label_nodes[node] = node

nx.draw_networkx_labels(subG, pos, labels=label_nodes, font_size=9, font_weight='bold', ax=ax)

# Legend
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#ED7D31', markersize=12, label='RTX targets'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#70AD47', markersize=10, label='TAC targets'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#C00000', markersize=10, label='CTX targets'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF0000', markersize=14, label='GWAS risk genes'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#4472C4', markersize=8, label='MN DEGs'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10, framealpha=0.9)

ax.set_title('Protein-Protein Interaction Network: MN Disease Module\n'
             'with Drug Target Overlay',
             fontsize=14, fontweight='bold')
ax.axis('off')

plt.tight_layout()
plt.savefig(FIGURES_DIR / "Figure5_PPI_Network.png", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved Figure5_PPI_Network.png")

# ================================================================
# Figure 6: Summary Plot - Volcano + Mechanism
# ================================================================
print("Generating Figure 6: Summary...")

# Volcano plot of DE results
de = pd.read_csv(RESULTS_DIR / "DE_MGN_vs_Normal.csv")
de['significant'] = de['adj.P.Val'] < 0.05
de['-log10p'] = -np.log10(de['adj.P.Val'].clip(lower=1e-300))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left: Volcano plot
ax1.scatter(de['logFC'], de['-log10p'], c=['#C00000' if (r and f > 2) 
             else '#4472C4' if (r and f < -2)
             else '#4472C4' if r else '#D3D3D3' 
             for r, f in zip(de['significant'], de['logFC'])], 
        s=2, alpha=0.5)
ax1.axvline(x=0, color='gray', linestyle='-', alpha=0.3)
ax1.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.3)
ax1.set_xlabel('log2(Fold Change) MN vs Normal', fontsize=12)
ax1.set_ylabel('-log10(adj.P.Val)', fontsize=12)
ax1.set_title('MN vs Normal: Differentially Expressed Genes\n'
              f'({sum(de["significant"])} DEGs at FDR<0.05)', fontsize=13, fontweight='bold')

# Label top genes on volcano
top_genes = de[de['significant']].nsmallest(10, 'adj.P.Val')
for _, row in top_genes.iterrows():
    ax1.annotate(row['Gene'], (row['logFC'], row['-log10p']),
                fontsize=6, alpha=0.8)

# Right: Mechanism summary
ax2.axis('off')
mechanism_text = """
MULTI-LAYER NETWORK MEDICINE FINDINGS

1. MN Disease Module:
   • 3,063 DEGs (87 MGN vs 11 Normal glomeruli)
   • 5 GWAS risk loci (PLA2R1, NFKB1, IRF4, HLA)
   • 2,176 proteins in PPI disease module

2. Drug-Disease Proximity:
   • RTX: z = -5.51  (strongest primary proximity)
   • CTX: z = -2.61  (primary proximity)
   • TAC: z = -1.85  (primary proximity)

3. Methodological Robustness:
   • RTX integrated proximity retained with
     degree-binned and STRING-threshold checks
   • CD20-only model not evaluable in final LCC
   • Full expanded RTX model proximal across
     DEG, GWAS and integrated disease sets

4. Mechanistic Interpretation:
   • RTX targets: CD20 + BCR signaling +
     FcγR effectors + complement
   • T-cell/NK-cell regulatory layer included
   • RTX+TAC remains a testable topology-derived
     combination hypothesis, not proof of synergy
"""

ax2.text(0.1, 0.5, mechanism_text, fontsize=10, fontfamily='monospace',
        va='center', ha='left', linespacing=1.5)
ax2.set_title('Network Medicine Analysis Summary', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(FIGURES_DIR / "Figure6_Summary.png", dpi=300, bbox_inches='tight')
plt.close()
print("  Saved Figure6_Summary.png")

print("\nAll figures generated successfully!")
print(f"Output directory: {FIGURES_DIR}")
for f in sorted(os.listdir(FIGURES_DIR)):
    size = os.path.getsize(FIGURES_DIR / f) / 1024
    print(f"  {f} ({size:.0f} KB)")
