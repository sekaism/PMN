#!/usr/bin/env python3
"""
Step 2: Build MN Disease Module + Collect Drug Targets
======================================================
1. Load MN DEGs from Step 1
2. Add GWAS risk genes (PLA2R1, NFKB1, IRF4, HLA-DQA1, HLA-DRB1)
3. Expand to disease module via STRING PPI (using networkx)
4. Collect drug targets for RTX, TAC, CTX
5. Save all data for Step 3 network medicine analysis
"""

import pandas as pd
import numpy as np
import networkx as nx
import requests
import json
import os
import pickle
import re
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
SCRIPTS_DIR = BASE_DIR / "scripts"

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 60)
print("Loading MN DEGs from Step 1...")

# Load MGN vs Normal DEGs
de_mn = pd.read_csv(RESULTS_DIR / "DE_MGN_vs_Normal.csv")
print(f"DEGs loaded: {len(de_mn)} total probes")

# Get significant DEGs
sig_de = de_mn[de_mn['adj.P.Val'] < 0.05].copy()
print(f"Significant DEGs: {len(sig_de)}")
print(f"  Up in MN: {(sig_de['logFC'] > 0).sum()}")
print(f"  Down in MN: {(sig_de['logFC'] < 0).sum()}")

# Get gene symbols (the 'Gene' column from limma output)
deg_genes = set(sig_de['Gene'].dropna().unique())
print(f"Unique gene symbols: {len(deg_genes)}")

# ================================================================
# MN GWAS risk genes (Xie et al. 2020, Nature Communications)
# ================================================================
print("\n" + "=" * 60)
print("Adding MN GWAS risk genes...")
gwas_genes = {
    "PLA2R1",    # Primary antigen, chr2q24.2, OR=2.25
    "NFKB1",     # NF-κB subunit, chr4q24, OR=1.25
    "IRF4",      # Interferon regulatory factor, chr6p25.3, OR=1.29
    "HLA-DQA1",  # MHC class II, chr6p21, OR=2.41
    "HLA-DRB1",  # MHC class II, chr6p21
}
print(f"GWAS risk genes: {gwas_genes}")

# Combine DEGs with GWAS genes
mn_seed_genes = deg_genes | gwas_genes
print(f"Combined seed genes (DEGs + GWAS): {len(mn_seed_genes)}")

# ================================================================
# STRING PPI Network Construction
# ================================================================
print("\n" + "=" * 60)
print("Querying STRING PPI network for seed genes...")

# Use STRING API to get interactions
# STRING API: https://string-db.org/api/
STRING_URL = "https://string-db.org/api"

def query_string_network(gene_list, species=9606, confidence=0.4):
    """Query STRING API for protein-protein interactions."""
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "required_score": int(confidence * 1000),
        "network_type": "physical",
        "limit": 50  # Number of interactors per protein
    }
    
    # Step 1: Get interaction partners
    url = f"{STRING_URL}/json/interaction_partners"
    try:
        r = requests.post(url, data=params, timeout=30)
        if r.status_code == 200:
            return r.json()
        else:
            print(f"  STRING API error: {r.status_code}")
            return None
    except Exception as e:
        print(f"  STRING API exception: {e}")
        return None

# Query in batches (STRING API has limits)
all_interactions = []
batch_size = 50
gene_list = sorted(mn_seed_genes)

# Only query a subset to start (top DEGs + GWAS)
top_deg = set(sig_de.nsmallest(200, 'adj.P.Val')['Gene'].dropna())
query_genes = list(top_deg | gwas_genes)
print(f"Querying STRING for {len(query_genes)} top genes...")

# Split into batches
for i in range(0, len(query_genes), batch_size):
    batch = query_genes[i:i+batch_size]
    result = query_string_network(batch)
    if result:
        all_interactions.extend(result)
    print(f"  Batch {i//batch_size + 1}/{(len(query_genes)-1)//batch_size + 1}: {len(batch)} genes, got {len(result) if result else 0} interactions")

print(f"Total interactions from STRING: {len(all_interactions)}")

# Build network
G = nx.Graph()
edge_list = []
for item in all_interactions:
    n1 = item.get("preferredName_A", item.get("stringId_A", ""))
    n2 = item.get("preferredName_B", item.get("stringId_B", ""))
    score = item.get("score", 0)
    if n1 and n2 and score >= 0.4:
        G.add_edge(n1, n2, weight=score)
        edge_list.append((n1, n2, score))

print(f"\nPPI network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Extract largest connected component (LCC)
components = list(nx.connected_components(G))
lcc = max(components, key=len) if components else set()
G_lcc = G.subgraph(lcc).copy()
print(f"LCC: {G_lcc.number_of_nodes()} nodes, {G_lcc.number_of_edges()} edges")

# Identify disease module: seed genes that are in the LCC
seed_in_lcc = mn_seed_genes & lcc
print(f"Seed genes in LCC: {len(seed_in_lcc)}")

# ================================================================
# Drug Target Collection
# ================================================================
print("\n" + "=" * 60)
print("Collecting drug targets...")

# --- RTX targets ---
# Direct target from DrugBank: MS4A1 (CD20)
# Fcγ receptors (ADCC effector function): FCGR1A, FCGR2A, FCGR2B, FCGR3A, FCGR3B
# Complement components (CDC): C1QA, C1QB, C1QC, C1R, C1S
# BCR signaling pathway (downstream of CD20)
# T cell modulation pathway (new findings 2024-2025)

rtx_targets = {
    # Direct target
    "MS4A1",  # CD20
    
    # Fcγ receptors (ADCC)
    "FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B",
    
    # Complement (CDC)
    "C1QA", "C1QB", "C1QC", "C1R", "C1S",
    
    # BCR signaling (downstream of CD20 engagement)
    "SYK", "LYN", "BTK", "PLCG2", "PIK3CA", "PIK3CB", "PIK3CD",
    "AKT1", "AKT2", "MTOR", "MAPK1", "MAPK3", "MAPK8", "NFKB1",
    "RELA", "IKBKB", "IKBKG", "CHUK",
    
    # T cell regulation (RTX-induced Treg pathway)
    "TGFB1", "TGFBR1", "TGFBR2", "FOXP3", "IL2RA", "CD4",
    
    # NK cell interaction (ADCC effector)
    "CD16A", "PRF1", "GZMB", "IFNG", "TNFSF10",
    
    # B cell markers and function
    "CD19", "CD22", "CD79A", "CD79B", "BLK", "BANK1",
    
    # Additional immune modulation
    "IL10", "IL6", "TNF", "IL2",
}

# Remove invalid entries and map CD16A to FCGR3A
rtx_targets.discard("CD16A")
rtx_targets.add("FCGR3A")

print(f"RTX targets: {len(rtx_targets)}")
print(f"  Direct targets: MS4A1 (CD20)")
print(f"  Effector targets: FcγRs + Complement + BCR pathway")

# --- TAC (Tacrolimus) targets ---
tac_targets = {
    "FKBP1A",    # FK506-binding protein (primary target)
    "FKBP1B",    # FK506-binding protein
    "PPP3CA",    # Calcineurin catalytic subunit Aα
    "PPP3CB",    # Calcineurin catalytic subunit Aβ
    "PPP3CC",    # Calcineurin catalytic subunit Aγ
    "PPP3R1",    # Calcineurin regulatory subunit B
    "PPP3R2",    # Calcineurin regulatory subunit B
    "NFATC1",    # NF-AT transcription factor
    "NFATC2",    # NF-AT transcription factor
    "NFATC3",    # NF-AT transcription factor
    "NFATC4",    # NF-AT transcription factor
    "IL2",       # IL-2 (downstream)
    "IL2RA",     # IL-2 receptor α
    "IL2RB",     # IL-2 receptor β
    "IL2RG",     # IL-2 receptor γ
    "JAK1", "JAK3", "STAT5A", "STAT5B",  # JAK-STAT downstream
    "MYC", "CCND1", "CDK4",  # Cell cycle targets
}
print(f"TAC (tacrolimus) targets: {len(tac_targets)}")

# --- CTX (Cyclophosphamide) targets ---
ctx_targets = {
    "ALDH1A1",  # Aldehyde dehydrogenase (activates CTX)
    "CYP2B6",   # Cytochrome P450 (activates CTX)
    "CYP3A4",   # Cytochrome P450
    "CYP2C9",   # Cytochrome P450
    "TP53",     # p53 (DNA damage response)
    "CDKN1A",   # p21 (cell cycle arrest)
    "BAX",      # Apoptosis
    "BCL2",     # Anti-apoptosis
    "CASP3",    # Caspase 3
    "CASP8",    # Caspase 8
    "CASP9",    # Caspase 9
    "ATM",      # DNA damage sensor
    "ATR",      # DNA damage sensor
    "CHEK1",    # Checkpoint kinase
    "CHEK2",    # Checkpoint kinase
    "FAS",      # Death receptor
    "FASLG",    # Fas ligand
    "GSTA1",    # Glutathione S-transferase (detoxification)
    "GSTP1",    # Glutathione S-transferase
    "MGMT",     # DNA repair
}
print(f"CTX (cyclophosphamide) targets: {len(ctx_targets)}")

# Save drug targets
drug_targets = {
    "RTX": rtx_targets,
    "TAC": tac_targets,
    "CTX": ctx_targets,
}

for drug, targets in drug_targets.items():
    # Save as CSV
    pd.DataFrame({"Gene": sorted(targets)}).to_csv(
        RESULTS_DIR / f"{drug}_targets.csv", index=False
    )
    print(f"  Saved {drug}_targets.csv ({len(targets)} genes)")

# ================================================================
# Expand drug targets via STRING (1st order neighbors in PPI)
# ================================================================
print("\n" + "=" * 60)
print("Expanding drug targets via STRING PPI...")

def expand_via_string(targets, network, max_neighbors=10):
    """Expand target set by adding directly connected neighbors in the PPI."""
    expanded = set(targets)
    for node in targets:
        if node in network:
            neighbors = list(network.neighbors(node))
            # Add top weighted neighbors
            weighted = [(n, network.get_edge_data(node, n).get('weight', 0)) 
                       for n in neighbors]
            weighted.sort(key=lambda x: -x[1])
            for n, w in weighted[:max_neighbors]:
                if w >= 0.4:
                    expanded.add(n)
    return expanded

# Expand each drug's targets
drug_targets_expanded = {}
for drug, targets in drug_targets.items():
    expanded = expand_via_string(targets, G_lcc)
    drug_targets_expanded[drug] = expanded
    print(f"  {drug}: {len(targets)} → {len(expanded)} (expanded)")

# ================================================================
# Save all data for Step 3
# ================================================================
print("\n" + "=" * 60)
print("Saving data for Step 3...")

# Save the PPI network
nx.write_graphml(G_lcc, DATA_DIR / "mn_ppi_network.graphml")
print(f"  PPI network saved: {DATA_DIR / 'mn_ppi_network.graphml'}")

# Save disease module genes
disease_module = {
    "seed_genes": sorted(mn_seed_genes),
    "deg_genes": sorted(deg_genes),
    "gwas_genes": sorted(gwas_genes),
    "lcc_genes": sorted(lcc),
}
with open(DATA_DIR / "disease_module.pkl", "wb") as f:
    pickle.dump(disease_module, f)
print(f"  Disease module saved: {DATA_DIR / 'disease_module.pkl'}")

# Save drug targets
with open(DATA_DIR / "drug_targets.pkl", "wb") as f:
    pickle.dump(drug_targets_expanded, f)
print(f"  Drug targets saved: {DATA_DIR / 'drug_targets.pkl'}")

# Save DEGs for reference
# Top 1000 DEGs for network analysis
top_1000_genes = set(sig_de.nsmallest(1000, 'adj.P.Val')['Gene'].dropna())
with open(DATA_DIR / "top1000_degs.pkl", "wb") as f:
    pickle.dump(top_1000_genes, f)

# Also save the full significant gene list
sig_genes = set(sig_de['Gene'].dropna())
with open(DATA_DIR / "sig_degs.pkl", "wb") as f:
    pickle.dump(sig_genes, f)

print(f"\nStep 2 complete! Ready for Step 3 (network medicine analysis).")
print(f"\nSummary:")
print(f"  MN DEGs (significant): {len(deg_genes)}")
print(f"  GWAS risk genes: {len(gwas_genes)}")
print(f"  PPI network LCC: {G_lcc.number_of_nodes()} nodes, {G_lcc.number_of_edges()} edges")
print(f"  RTX targets: {len(drug_targets_expanded['RTX'])}")
print(f"  TAC targets: {len(drug_targets_expanded['TAC'])}")
print(f"  CTX targets: {len(drug_targets_expanded['CTX'])}")
