#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import seaborn as sns

# ===============================
# CONFIG
# ===============================
DOCK_DIR = Path("results/docking")
RECEPTOR_DIR = Path("receptors")
OUT_BASE = Path("results/all_complexes")
PLIP_BASE = OUT_BASE / "plip"
FIG_DIR = OUT_BASE / "figures"

N_CPU = max(1, cpu_count() - 2)

OUT_BASE.mkdir(parents=True, exist_ok=True)
PLIP_BASE.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ===============================
# UTIL FUNCTIONS
# ===============================

def extract_best_pose(vina_pdbqt, out_pdbqt):
    with open(vina_pdbqt) as inp, open(out_pdbqt, "w") as out:
        write = False
        for line in inp:
            if line.startswith("MODEL 1"):
                write = True
                continue
            if line.startswith("MODEL"):
                break
            if write and not line.startswith("ENDMDL"):
                out.write(line)

def pdbqt_to_pdb(pdbqt_file, pdb_file):
    with open(pdbqt_file) as inp, open(pdb_file, "w") as out:
        for line in inp:
            if line.startswith(("ATOM", "HETATM")):
                out.write(line[:66] + "\n")

def make_complex(receptor_pdbqt, ligand_pdb, complex_pdb):
    with open(complex_pdb, "w") as out:
        with open(receptor_pdbqt) as r:
            for l in r:
                if l.startswith(("ATOM", "HETATM")):
                    out.write(l[:66] + "\n")
        with open(ligand_pdb) as l:
            for ln in l:
                out.write(ln)

# ===============================
# STEP 1: COMPLEX GENERATION
# ===============================

print("?? Generating all complexes...")

for receptor_dir in DOCK_DIR.iterdir():
    if not receptor_dir.is_dir():
        continue

    receptor = receptor_dir.name
    receptor_pdbqt = RECEPTOR_DIR / f"{receptor}.pdbqt"
    out_rec_dir = OUT_BASE / receptor
    out_rec_dir.mkdir(exist_ok=True)

    for file in receptor_dir.glob("*_out.pdbqt"):
        ligand = file.stem.replace("_out", "")

        best_pose = out_rec_dir / f"{ligand}_best.pdbqt"
        ligand_pdb = out_rec_dir / f"{ligand}.pdb"
        complex_pdb = out_rec_dir / f"{ligand}_complex.pdb"

        extract_best_pose(file, best_pose)
        pdbqt_to_pdb(best_pose, ligand_pdb)
        make_complex(receptor_pdbqt, ligand_pdb, complex_pdb)

print("? Complex generation done")

# ===============================
# STEP 2: PARALLEL PLIP
# ===============================

def run_plip_task(args):
    complex_pdb, outdir = args
    subprocess.run(
        ["plip", "-f", str(complex_pdb), "-o", str(outdir), "--xml", "--silent"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

print(f"?? Running PLIP in parallel using {N_CPU} CPUs...")

tasks = []
for receptor_dir in OUT_BASE.iterdir():
    if not receptor_dir.is_dir():
        continue
    receptor = receptor_dir.name

    for complex_pdb in receptor_dir.glob("*_complex.pdb"):
        ligand = complex_pdb.stem.replace("_complex", "")
        outdir = PLIP_BASE / receptor / ligand
        outdir.mkdir(parents=True, exist_ok=True)
        tasks.append((complex_pdb, outdir))

with Pool(N_CPU) as pool:
    pool.map(run_plip_task, tasks)

print("? Parallel PLIP completed")

# ===============================
# STEP 3: PARSE XML
# ===============================

print("?? Parsing PLIP outputs...")

PLIP_BLOCKS = {
    "hydrogen_bonds": "hbond",
    "hydrophobic_interactions": "hydrophobic",
    "salt_bridges": "salt",
    "pi_stacks": "pi_stack",
    "pi_cation_interactions": "pi_cation",
    "halogen_bonds": "halogen",
    "metal_complexes": "metal"
}

rows = []

for receptor_dir in PLIP_BASE.iterdir():
    receptor = receptor_dir.name

    for ligand_dir in receptor_dir.iterdir():
        ligand = ligand_dir.name
        xmls = list(ligand_dir.glob("*_report.xml"))

        if not xmls:
            continue

        try:
            root = ET.parse(xmls[0]).getroot()
        except:
            continue

        for block, itype in PLIP_BLOCKS.items():
            node = root.find(f".//{block}")
            if node is None:
                continue

            for item in node:
                rows.append([
                    receptor,
                    ligand,
                    itype,
                    item.findtext("restype", ""),
                    item.findtext("resnr", "")
                ])

df = pd.DataFrame(rows, columns=["receptor","ligand","type","residue","resnum"])
df["residue_id"] = df["residue"] + df["resnum"]

df.to_csv(OUT_BASE / "plip_interactions_all.csv", index=False)

print("? Interaction table created")

# ===============================
# STEP 4: HEATMAPS
# ===============================

print("?? Generating heatmaps...")

heat = df.groupby(["ligand","residue_id"]).size().unstack(fill_value=0)

if heat.shape[1] > 50:
    heat = heat[heat.sum().sort_values(ascending=False).head(30).index]

plt.figure(figsize=(12, 8))
sns.heatmap(heat, cmap="viridis")
plt.title("All Ligand Residue Interaction Heatmap")
plt.tight_layout()
plt.savefig(FIG_DIR / "all_residue_heatmap.png", dpi=300)
plt.close()

print("? Heatmap generated")

# ===============================
# STEP 5: INTERACTION RANKING
# ===============================

print("?? Computing interaction-based ranking...")

weights = {
    "hbond": 3,
    "salt": 3,
    "pi_stack": 2,
    "pi_cation": 2,
    "hydrophobic": 1,
    "halogen": 2,
    "metal": 3
}

df["score"] = df["type"].map(weights).fillna(1)

ranking = df.groupby("ligand")["score"].sum().sort_values(ascending=False)
ranking.to_csv(OUT_BASE / "interaction_ranking.csv")

print("? Interaction ranking done")

# ===============================
# STEP 6: ML DATASET
# ===============================

print("?? Creating ML-ready dataset...")

fp = df.groupby(["ligand","type"]).size().unstack(fill_value=0)

# Add residue counts
res_fp = df.groupby(["ligand","residue_id"]).size().unstack(fill_value=0)

ml_df = pd.concat([fp, res_fp], axis=1)
ml_df.to_csv(OUT_BASE / "ml_dataset.csv")

print("? ML dataset ready")

print("?? ADVANCED ALL-COMPLEX PIPELINE COMPLETED")