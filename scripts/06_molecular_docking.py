#!/usr/bin/env python3
"""Step 6: Molecular Docking - solid pipeline."""

import os, subprocess, requests, csv, re
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']
import matplotlib.pyplot as plt

BASE = Path("/workspaces/platform/MN")
D = BASE / "data" / "docking"
R = BASE / "results"; F = BASE / "figures"
PD = D / "pdb"; os.makedirs(PD, exist_ok=True); os.makedirs(D, exist_ok=True)

def dl_pdb(pdb, out):
    path = out / f"{pdb}.pdb"
    if not path.exists():
        r = requests.get(f"https://files.rcsb.org/download/{pdb}.pdb", timeout=30)
        if r.status_code == 200: path.write_text(r.text)
    return path if path.exists() else None

def lig_pdbqt(drug, smiles, ddir):
    pdbqt = ddir / f"{drug}.pdbqt"
    if not pdbqt.exists():
        sdf = ddir / f"{drug}.sdf"
        subprocess.run(f'obabel -:"{smiles}" --gen3D -O {sdf} 2>/dev/null', shell=True)
        if sdf.exists():
            subprocess.run(f'obabel {sdf} -O {pdbqt} --gen3D 2>/dev/null', shell=True)
    return pdbqt if pdbqt.exists() else None

def prep_rec(pdb, pdbqt):
    if not pdbqt.exists():
        clean = pdb.with_suffix(".clean.pdb")
        with open(pdb) as f: clean.write_text("".join(l for l in f if l.startswith("ATOM")))
        subprocess.run(f'obabel {clean} -O {pdbqt} -xr 2>/dev/null', shell=True)
    return pdbqt.exists()

def run_vina(rec, lig, out, center, size=(22,22,22)):
    r = subprocess.run(
        f'vina --receptor {rec} --ligand {lig} --out {out} '
        f'--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} '
        f'--size_x {size[0]} --size_y {size[1]} --size_z {size[2]} '
        f'--exhaustiveness 8 2>/dev/null',
        shell=True, capture_output=True, text=True)
    for line in r.stdout.split('\n'):
        m = re.match(r'^\s+1\s+(-?\d+\.\d+)', line)
        if m: return float(m.group(1))
    # Fallback: try vina's stdout directly
    try:
        for line in r.stdout.split('\n'):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] == '1':
                return float(parts[1])
    except: pass
    return None

# Download PDBs
targets = {"FKBP1A":"1FKB","TP53":"1TUP","AKT1":"1UNP","NFKB1":"1SVC",
           "PPP3CA":"1AUI","CASP3":"1CP3"}
pdbs = {}
for n, pid in targets.items():
    p = dl_pdb(pid, PD)
    if p: pdbs[n] = pid

# Ligands
ligs = {"TAC":"CCC1CC2CC(=O)C(=O)N3C2C(C(C3C(=O)C(C(C1C)OC)OC(=O)CC4CCCC4)OC)OC",
        "CTX":"ClCCN(CCCl)P1(=O)NCCCO1"}
lig_pq = {}
for d, s in ligs.items():
    p = lig_pdbqt(d, s, D)
    if p: lig_pq[d] = p

# Docking
results = {}
for drug, lig in lig_pq.items():
    results[drug] = {}
    for prot, pid in pdbs.items():
        pdb = PD / f"{pid}.pdb"
        rpq = D / f"{pid}.pdbqt"; opq = D / f"{drug}_{pid}_dock.pdbqt"
        if not prep_rec(pdb, rpq): continue
        centers = {"FKBP1A":(20,10,15),"TP53":(-5,20,5),"AKT1":(15,15,20),
                   "NFKB1":(20,15,10),"PPP3CA":(10,5,20),"CASP3":(25,15,20)}
        e = run_vina(rpq, lig, opq, centers.get(prot,(15,15,15)))
        results[drug][prot] = e
        print(f"  {drug} → {prot} ({pid}): {f'{e:.2f}' if e else 'FAIL'}")

# Heatmap
prots = list(pdbs.keys())
data = np.full((2, len(prots)), np.nan)
for i,d in enumerate(["TAC","CTX"]):
    for j,p in enumerate(prots):
        e = results.get(d,{}).get(p)
        if e: data[i,j] = e

fig, ax = plt.subplots(figsize=(10,5))
cmap = plt.cm.RdYlGn_r.copy()
im = ax.imshow(data, cmap=cmap, aspect="auto", vmin=-10, vmax=0)
ax.set_xticks(range(len(prots)))
ax.set_xticklabels(prots, fontsize=10, rotation=30, ha="right")
ax.set_yticks([0,1]); ax.set_yticklabels(["Tacrolimus", "Cyclophosphamide"], fontsize=10)
for i in range(2):
    for j in range(len(prots)):
        if not np.isnan(data[i,j]):
            c = "white" if data[i,j] < -5 else "black"
            ax.text(j,i,f"{data[i,j]:.1f}",ha="center",va="center",fontweight="bold",fontsize=9,c=c)
plt.colorbar(im, ax=ax, label="Binding energy (kcal/mol)", shrink=0.7)
ax.set_title("Molecular Docking Validation\n(more negative = stronger binding)", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(F / "Figure7_Molecular_Docking.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"\n✓ Figure7_Molecular_Docking.png saved")

# CSV
import csv
with open(R / "docking_results.csv","w") as f:
    w=csv.writer(f); w.writerow(["Drug","Protein","PDB","Energy"])
    for d in ["TAC","CTX"]:
        for p in prots:
            e=results.get(d,{}).get(p)
            if e: w.writerow([d,p,pdbs[p],f"{e:.1f}"])
print("✓ docking_results.csv saved")
