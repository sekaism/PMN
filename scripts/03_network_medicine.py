#!/usr/bin/env python3
"""
Step 3: Network Medicine Analysis
==================================
Manual implementation (more robust than NetMedPy all_pair_distances).
Computes drug-disease proximity and drug-drug separation.
"""

import networkx as nx
import numpy as np
import pandas as pd
import pickle
import os
import random
from pathlib import Path
from scipy import stats

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

def main():
    print("=" * 60)
    print("STEP 3: Network Medicine Analysis")
    print("=" * 60)
    
    # ================================================================
    # 1. Load data
    # ================================================================
    print("\n1. Loading data...")
    G = nx.read_graphml(DATA_DIR / "mn_ppi_network.graphml")
    G = G.to_undirected()
    
    # Ensure LCC
    components = list(nx.connected_components(G))
    lcc_nodes = max(components, key=len)
    G = G.subgraph(lcc_nodes).copy()
    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    with open(DATA_DIR / "disease_module.pkl", "rb") as f:
        disease_module = pickle.load(f)
    
    with open(DATA_DIR / "drug_targets.pkl", "rb") as f:
        drug_targets = pickle.load(f)
    
    # Also load DEGs from CSV
    de = pd.read_csv(RESULTS_DIR / "DE_MGN_vs_Normal.csv")
    sig_de = de[de['adj.P.Val'] < 0.05]
    deg_genes = set(sig_de['Gene'].dropna())
    gwas_genes = {"PLA2R1", "NFKB1", "IRF4", "HLA-DQA1", "HLA-DRB1"}
    
    # Filter to PPI nodes
    all_ppi = set(G.nodes())
    deg_ppi = deg_genes & all_ppi
    gwas_ppi = gwas_genes & all_ppi
    disease_ppi = deg_ppi | gwas_ppi
    
    drug_ppi = {}
    for drug, targets in drug_targets.items():
        drug_ppi[drug] = set(targets) & all_ppi
    
    print(f"  DEGs in PPI: {len(deg_ppi)}")
    print(f"  GWAS in PPI: {len(gwas_ppi)}")
    for d in drug_ppi:
        print(f"  {d} targets in PPI: {len(drug_ppi[d])}")
    
    # ================================================================
    # 2. Compute shortest path distances
    # ================================================================
    print("\n2. Computing shortest path distances...")
    
    # We only need distances between drug targets and disease genes
    # This is much faster than all-pair computation
    nodes_of_interest = set()
    for d in drug_ppi:
        nodes_of_interest.update(drug_ppi[d])
    nodes_of_interest.update(disease_ppi)
    
    print(f"  Nodes of interest: {len(nodes_of_interest)}")
    
    # Compute shortest path lengths from each drug target to all nodes
    distance_cache = {}
    
    def get_distances(source_nodes, target_nodes):
        """Get shortest path distances between two sets of nodes."""
        key = (frozenset(source_nodes), frozenset(target_nodes))
        if key in distance_cache:
            return distance_cache[key]
        
        result = {}
        for s in source_nodes:
            if s not in result:
                try:
                    lengths = nx.single_source_shortest_path_length(G, s)
                    result[s] = lengths
                except:
                    result[s] = {}
        distance_cache[key] = result
        return result
    
    # ================================================================
    # 3. Drug-Disease PROXIMITY
    # ================================================================
    print("\n3. Drug-Disease PROXIMITY Analysis")
    print("-" * 40)
    
    def proximity_zscore(drug_set, disease_set, n_perm=1000):
        """Calculate proximity z-score via permutation test."""
        drug_set = [n for n in drug_set if n in G]
        disease_set = [n for n in disease_set if n in G]
        
        if len(drug_set) < 2 or len(disease_set) < 2:
            return None
        
        # Observed: average minimum distance from each drug target to disease
        drug_distances = get_distances(drug_set, disease_set)
        obs_min_dists = []
        for s in drug_set:
            dists_to_disease = [drug_distances[s].get(t, 100) for t in disease_set 
                              if t in drug_distances[s]]
            if dists_to_disease:
                obs_min_dists.append(min(dists_to_disease))
        observed = np.mean(obs_min_dists) if obs_min_dists else 100
        
        # Permutation: random drug targets
        all_nodes = list(G.nodes())
        perm_means = []
        for _ in range(n_perm):
            perm_drugs = random.sample(all_nodes, min(len(drug_set), len(all_nodes)))
            perm_distances = get_distances(perm_drugs, disease_set)
            perm_mins = []
            for s in perm_drugs:
                dists = [perm_distances[s].get(t, 100) for t in disease_set 
                        if t in perm_distances[s]]
                if dists:
                    perm_mins.append(min(dists))
            if perm_mins:
                perm_means.append(np.mean(perm_mins))
        
        z = (observed - np.mean(perm_means)) / np.std(perm_means) if np.std(perm_means) > 0 else 0
        
        return {
            "z_score": z,
            "observed": observed,
            "perm_mean": np.mean(perm_means),
            "perm_std": np.std(perm_means),
            "n_drug": len(drug_set),
            "n_disease": len(disease_set)
        }
    
    print("\n  Proximity to MN Disease Module:")
    prox_results = {}
    
    for drug in ["RTX", "TAC", "CTX"]:
        targets = list(drug_ppi[drug])
        if len(targets) < 2:
            print(f"  {drug}: SKIP (too few targets)")
            continue
        
        print(f"\n  {drug} ({len(targets)} targets):")
        
        # vs DEGs
        res_deg = proximity_zscore(targets, list(deg_ppi), n_perm=500)
        if res_deg:
            print(f"    vs DEGs: z = {res_deg['z_score']:.2f}, "
                  f"obs={res_deg['observed']:.2f}, perm={res_deg['perm_mean']:.2f}")
        
        # vs GWAS
        res_gwas = proximity_zscore(targets, list(gwas_ppi), n_perm=500)
        if res_gwas:
            print(f"    vs GWAS: z = {res_gwas['z_score']:.2f}, "
                  f"obs={res_gwas['observed']:.2f}, perm={res_gwas['perm_mean']:.2f}")
        
        # vs Disease module
        res_all = proximity_zscore(targets, list(disease_ppi), n_perm=500)
        if res_all:
            print(f"    vs ALL:  z = {res_all['z_score']:.2f}, "
                  f"obs={res_all['observed']:.2f}, perm={res_all['perm_mean']:.2f}")
        
        prox_results[drug] = {"DEGs": res_deg, "GWAS": res_gwas, "ALL": res_all}
    
    # ================================================================
    # 4. Drug-Drug SEPARATION
    # ================================================================
    print("\n" + "=" * 60)
    print("4. Drug-Drug SEPARATION")
    print("-" * 40)
    
    def separation(set_a, set_b):
        """Calculate network separation between two drug target sets."""
        set_a = [n for n in set_a if n in G]
        set_b = [n for n in set_b if n in G]
        
        if len(set_a) < 2 or len(set_b) < 2:
            return None
        
        # Compute distances
        dists_a = get_distances(set_a, set_a)
        dists_b = get_distances(set_b, set_b)
        dists_ab = get_distances(set_a, set_b)
        
        # d_aa: mean shortest distance within set A
        d_aa_vals = []
        for s in set_a:
            vals = [dists_a[s].get(t, 100) for t in set_a if t != s and t in dists_a[s]]
            if vals:
                d_aa_vals.append(min(vals))
        d_aa = np.mean(d_aa_vals) if d_aa_vals else 100
        
        # d_bb: mean shortest distance within set B
        d_bb_vals = []
        for s in set_b:
            vals = [dists_b[s].get(t, 100) for t in set_b if t != s and t in dists_b[s]]
            if vals:
                d_bb_vals.append(min(vals))
        d_bb = np.mean(d_bb_vals) if d_bb_vals else 100
        
        # d_ab: mean shortest distance between A and B
        d_ab_vals = []
        for s in set_a:
            vals = [dists_ab[s].get(t, 100) for t in set_b if t in dists_ab[s]]
            if vals:
                d_ab_vals.append(min(vals))
        d_ab = np.mean(d_ab_vals) if d_ab_vals else 100
        
        s = d_ab - (d_aa + d_bb) / 2
        
        return {"separation": s, "d_ab": d_ab, "d_aa": d_aa, "d_bb": d_bb,
                "n_a": len(set_a), "n_b": len(set_b)}
    
    print("\n  Drug Pair Separation:")
    sep_results = {}
    pairs = [("RTX", "TAC"), ("RTX", "CTX"), ("TAC", "CTX")]
    
    for d1, d2 in pairs:
        res = separation(drug_ppi[d1], drug_ppi[d2])
        if res:
            sep_results[f"{d1}_{d2}"] = res
            print(f"\n  {d1} vs {d2}:")
            print(f"    Separation (s) = {res['separation']:.2f}")
            print(f"    d_ab = {res['d_ab']:.2f}, d_aa = {res['d_aa']:.2f}, d_bb = {res['d_bb']:.2f}")
            
            if res['separation'] < -0.5:
                print(f"    → OVERLAPPING targets")
            elif res['separation'] < 0.5:
                print(f"    → INTERMEDIATE")
            else:
                print(f"    → SEPARATED (Complementary Exposure potential)")
    
    # ================================================================
    # 5. Complementary Exposure Test
    # ================================================================
    print("\n" + "=" * 60)
    print("5. Complementary Exposure Assessment")
    print("-" * 40)
    
    for d1, d2 in pairs:
        key = f"{d1}_{d2}"
        if key not in sep_results:
            continue
        
        p1 = prox_results.get(d1, {}).get("ALL", {})
        p2 = prox_results.get(d2, {}).get("ALL", {})
        s_val = sep_results[key]["separation"]
        
        print(f"\n  {d1}+{d2}:")
        if p1:
            print(f"    {d1}→MN z = {p1.get('z_score', 0):.2f}")
        if p2:
            print(f"    {d2}→MN z = {p2.get('z_score', 0):.2f}")
        print(f"    Separation s = {s_val:.2f}")
        
        # Complementary Exposure: both close to disease (z < -1) but separated (s > 0)
        z1 = p1.get('z_score', 0) if p1 else 0
        z2 = p2.get('z_score', 0) if p2 else 0
        
        if z1 < -1 and z2 < -1 and s_val > 0:
            print(f"    ✓ COMPLEMENTARY EXPOSURE → synergistic potential!")
        elif z1 < -1 and z2 < -1 and s_val < 0:
            print(f"    → OVERLAPPING EXPOSURE → redundant mechanisms")
        else:
            print(f"    → Neither drug shows strong proximity to MN module")
    
    # ================================================================
    # 6. Save results
    # ================================================================
    print("\n6. Saving results...")
    
    # Proximity results
    prox_rows = []
    for drug in prox_results:
        for target_set in prox_results[drug]:
            r = prox_results[drug][target_set]
            if r:
                prox_rows.append({
                    "Drug": drug, "TargetSet": target_set,
                    "Z_score": r["z_score"], "Observed": r["observed"],
                    "PermMean": r["perm_mean"]
                })
    pd.DataFrame(prox_rows).to_csv(RESULTS_DIR / "drug_disease_proximity.csv", index=False)
    
    # Separation results
    sep_rows = []
    for pair in sep_results:
        r = sep_results[pair]
        sep_rows.append({
            "Pair": pair, "Drug1": pair.split("_")[0], "Drug2": pair.split("_")[1],
            "Separation": r["separation"], "d_ab": r["d_ab"],
            "d_aa": r["d_aa"], "d_bb": r["d_bb"]
        })
    pd.DataFrame(sep_rows).to_csv(RESULTS_DIR / "drug_drug_separation.csv", index=False)
    
    # Save combined
    with open(RESULTS_DIR / "network_medicine_results.pkl", "wb") as f:
        pickle.dump({"proximity": prox_results, "separation": sep_results}, f)
    
    print("  Saved: drug_disease_proximity.csv")
    print("  Saved: drug_drug_separation.csv")
    print("\n" + "=" * 60)
    print("STEP 3 COMPLETE!")
    print("=" * 60)

if __name__ == "__main__":
    main()
