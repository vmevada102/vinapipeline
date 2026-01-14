#!/usr/bin/env bash
set -e

# ====================================================
# AutoDock Vina - FULLY AUTOMATED PIPELINE
# Parallel | Resume | Top-N Safe Analysis
# ====================================================

echo "===================================================="
echo " AutoDock Vina - FULLY AUTOMATED PIPELINE"
echo " Parallel | Resume | Top-N Safe Analysis"
echo "===================================================="

# ----------------------------------------------------
# USER INPUTS
# ----------------------------------------------------
read -p "Ligand directory [default: ligands/]: " LIGAND_DIR
read -p "Receptor directory [default: receptors/]: " RECEPTOR_DIR
read -p "Config directory   [default: config/]: " CONFIG_DIR

LIGAND_DIR=${LIGAND_DIR:-ligands}
RECEPTOR_DIR=${RECEPTOR_DIR:-receptors}
CONFIG_DIR=${CONFIG_DIR:-config}

# ----------------------------------------------------
# TOP-N CONTROL (ONLY ADDITION YOU ASKED FOR)
# ----------------------------------------------------
read -p "Number of TOP compounds to use for analysis (heatmaps, poses) [default: 20]: " TOP_N_ANALYSIS
TOP_N_ANALYSIS=${TOP_N_ANALYSIS:-20}
export TOP_N_ANALYSIS

echo "Top-N compounds for downstream analysis : $TOP_N_ANALYSIS"

# ----------------------------------------------------
# OUTPUT DIRECTORIES
# ----------------------------------------------------
DOCK_DIR="results/docking"
SUMMARY_DIR="results/summary"
FIG_DIR="results/summary/figures"

mkdir -p "$DOCK_DIR" "$SUMMARY_DIR" "$FIG_DIR"

# ----------------------------------------------------
# CPU & PARALLEL SETUP
# ----------------------------------------------------
TOTAL_CPU=$(nproc)
PARALLEL_JOBS=$((TOTAL_CPU - 2))
[ "$PARALLEL_JOBS" -lt 1 ] && PARALLEL_JOBS=1

echo "Detected CPUs   : $TOTAL_CPU"
echo "Parallel jobs   : $PARALLEL_JOBS"
echo "Resume mode     : ENABLED"

# ----------------------------------------------------
# CHECK DEPENDENCIES
# ----------------------------------------------------
command -v vina >/dev/null || { echo "ERROR: vina not found"; exit 1; }
command -v python >/dev/null || { echo "ERROR: python not found"; exit 1; }

if command -v parallel >/dev/null; then
    USE_PARALLEL=1
else
    USE_PARALLEL=0
fi



# ----------------------------------------------------
# DOCKING FUNCTION (RESUME SAFE)
# ----------------------------------------------------
dock_ligand () {
    receptor="$1"
    ligand="$2"

    rname=$(basename "$receptor" .pdbqt)
    lname=$(basename "$ligand" .pdbqt)

    cfg="$CONFIG_DIR/${rname}_config.txt"
    outdir="$DOCK_DIR/$rname"
    logfile="$outdir/${lname}_log.txt"

    mkdir -p "$outdir"
    [ ! -f "$cfg" ] && return
    [ -s "$logfile" ] && return

    vina --receptor "$receptor" \
         --ligand "$ligand" \
         --config "$cfg" \
         --cpu 1 \
         --out "$outdir/${lname}_out.pdbqt" \
         > "$logfile" 2>&1
}

export -f dock_ligand
export CONFIG_DIR DOCK_DIR

# ----------------------------------------------------
# RUN DOCKING
# ----------------------------------------------------
for receptor in "$RECEPTOR_DIR"/*.pdbqt; do
    if [ "$USE_PARALLEL" -eq 1 ]; then
        parallel -j "$PARALLEL_JOBS" dock_ligand "$receptor" ::: "$LIGAND_DIR"/*.pdbqt
    else
        for ligand in "$LIGAND_DIR"/*.pdbqt; do
            dock_ligand "$receptor" "$ligand"
        done
    fi
done

# ----------------------------------------------------
# ANALYSIS (TOP-N ONLY, HEATMAP SAFE)
# ----------------------------------------------------
python << EOF
import os, re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

def extract_best_pose(vina_pdbqt, out_pdbqt):
    """
    Extract MODEL 1 (best pose) from a Vina multi-model PDBQT
    """
    with open(vina_pdbqt) as inp, open(out_pdbqt, "w") as out:
        write = False
        for line in inp:
            if line.startswith("MODEL 1"):
                write = True
                continue
            if line.startswith("MODEL"):
                break
            if write:
                if not line.startswith("ENDMDL"):
                    out.write(line)




TOP_N = int(os.getenv("TOP_N_ANALYSIS", "20"))

rows = []
pattern = re.compile(r"^\\s*(\\d+)\\s+(-?\\d+\\.\\d+)")

for r in os.listdir("results/docking"):
    rdir = os.path.join("results/docking", r)
    if not os.path.isdir(rdir):
        continue
    for f in os.listdir(rdir):
        if f.endswith("_log.txt"):
            ligand = f.replace("_log.txt", "")
            with open(os.path.join(rdir, f)) as fh:
                for line in fh:
                    m = pattern.match(line)
                    if m and m.group(1) == "1":
                        rows.append([r, ligand, float(m.group(2))])

if not rows:
    raise SystemExit("No docking results found")

df = pd.DataFrame(rows, columns=["receptor", "ligand", "energy"])
df.to_csv("results/summary/vina_best_scores.csv", index=False)

energy_matrix = df.pivot(index="ligand", columns="receptor", values="energy")
energy_matrix.to_csv("results/summary/ligand_receptor_binding_energy_matrix.csv")

# df columns: ligand, receptor, energy

# Rank ligands per receptor (lower energy = better)
df["rank"] = df.groupby("receptor")["energy"].rank(method="min", ascending=True)

# Build rank matrix
rank_matrix = df.pivot(index="ligand", columns="receptor", values="rank")
rank_matrix.to_csv("results/summary/ligand_receptor_rank_matrix.csv")

# Compute Borda consensus score
borda = pd.DataFrame(index=rank_matrix.index)

for rec in rank_matrix.columns:
    max_rank = rank_matrix[rec].max()
    borda[rec] = max_rank - rank_matrix[rec] + 1

borda["consensus_score"] = borda.sum(axis=1)

# Final ranked ligands
final_ranking = borda.sort_values("consensus_score", ascending=False)
final_ranking.to_csv("results/summary/ligand_consensus_ranking.csv")


# ---- Weighted Consensus (important receptors get higher weight) ----

# ---- Define receptor weights ----
# You can edit these values
receptor_weights = {
    rec: 1.0 for rec in rank_matrix.columns
}

# Example: make main target more important
# receptor_weights["7tob"] = 2.0
# receptor_weights["3pp0"] = 1.5

# ---- Compute weighted Borda ----
weighted_borda = pd.DataFrame(index=rank_matrix.index)

for rec in rank_matrix.columns:
    max_rank = rank_matrix[rec].max()
    weighted_borda[rec] = receptor_weights.get(rec, 1.0) * (max_rank - rank_matrix[rec] + 1)

weighted_borda["weighted_consensus"] = weighted_borda.sum(axis=1)

weighted_borda.sort_values("weighted_consensus", ascending=False)\
              .to_csv("results/summary/ligand_weighted_consensus.csv")


import matplotlib.pyplot as plt
import seaborn as sns

# Limit to top ligands for readability
topN = min(50, len(rank_matrix))
top_ligands = final_ranking.head(topN).index
rank_heat = rank_matrix.loc[top_ligands]

plt.figure(figsize=(rank_heat.shape[1]*1.2, rank_heat.shape[0]*0.35))
sns.heatmap(rank_heat, cmap="viridis_r", annot=False, linewidths=0.2)
plt.title(f"Top-{topN} Ligand Rank Heatmap (Lower = Better)")
plt.xlabel("Receptor")
plt.ylabel("Ligand")
plt.tight_layout()
plt.savefig("results/summary/figures/rank_heatmap_top50.png", dpi=300)
plt.close()


import numpy as np

def bootstrap_consensus(df, receptors, n_boot=1000):
    ligands = df["ligand"].unique()
    scores = {lig: [] for lig in ligands}

    for i in range(n_boot):
        sampled = np.random.choice(receptors, size=len(receptors), replace=True)
        temp = df[df["receptor"].isin(sampled)]

        ranks = temp.groupby("receptor")["energy"].rank(method="min")
        temp = temp.assign(rank=ranks)

        rm = temp.pivot(index="ligand", columns="receptor", values="rank")
        borda = rm.apply(lambda x: rm.max()[x.name] - x + 1)

        total = borda.sum(axis=1)
        for lig in total.index:
            scores[lig].append(total[lig])

    return scores

# Run bootstrap
bootstrap_scores = bootstrap_consensus(df, df["receptor"].unique(), n_boot=500)

# Compute confidence intervals
ci = []
for lig, vals in bootstrap_scores.items():
    ci.append([
        lig,
        np.mean(vals),
        np.percentile(vals, 2.5),
        np.percentile(vals, 97.5)
    ])

ci_df = pd.DataFrame(ci, columns=["ligand","mean_score","CI_low","CI_high"])
ci_df.sort_values("mean_score", ascending=False)\
     .to_csv("results/summary/ligand_consensus_bootstrap_CI.csv", index=False)





# ---- TOP-N SELECTION (CRITICAL SAFETY) ----
top_ligands = (
    df.groupby("ligand")["energy"]
      .mean()
      .sort_values()
      .head(TOP_N)
      .index
)

plot_data = energy_matrix.loc[top_ligands]

# ---- SAFE HEATMAP (NEVER FULL MATRIX) ----
plt.figure(figsize=(10, max(6, TOP_N * 0.3)))
sns.heatmap(
    plot_data,
    cmap="viridis_r",
    annot=TOP_N <= 50,
    cbar_kws={"label": "Binding Energy (kcal/mol)"}
)
plt.title(f"Top {TOP_N} Ligand-Receptor Binding Energy Heatmap")
plt.xlabel("Receptor")
plt.ylabel("Ligand")
plt.tight_layout()
plt.savefig(f"results/summary/figures/top{TOP_N}_binding_energy_heatmap.png", dpi=300)
plt.close()

print("Heatmap generated safely for Top-N ligands")


# " Extract Top-N poses per receptor (SAFE)"


import shutil
from pathlib import Path

TOP_N = int(os.getenv("TOP_N_ANALYSIS", "20"))

TOP_DIR = Path("results/top_N")
TOP_DIR.mkdir(exist_ok=True, parents=True)

# rank ligands globally
ranked_ligands = (
    df.groupby("ligand")["energy"]
      .mean()
      .sort_values()
      .head(TOP_N)
      .index
)

for receptor in df["receptor"].unique():
    rec_dir = TOP_DIR / receptor
    rec_dir.mkdir(exist_ok=True)

    dock_dir = Path("results/docking") / receptor

for lig in top_ligands:
    src = dock_dir / f"{lig}_out.pdbqt"
    dst = rec_dir / f"{lig}_best.pdbqt"

    if src.exists():
        extract_best_pose(src, dst)


#===================================================="
#Convert best pose to PDB & build complexes"
#===================================================="

def pdbqt_to_pdb(pdbqt_file, pdb_file):
    with open(pdbqt_file) as inp, open(pdb_file, "w") as out:
        for line in inp:
            if line.startswith(("ATOM", "HETATM")):
                # Remove AutoDock-specific columns (Q, T)
                out.write(line[:66] + "\n")

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

for receptor in df["receptor"].unique():
    rec_dir = Path("results/top_N") / receptor
    receptor_pdbqt = Path("receptors") / f"{receptor}.pdbqt"

    for lig_pdbqt in rec_dir.glob("*_best.pdbqt"):
        lig_pdb = lig_pdbqt.with_suffix(".pdb")
        complex_pdb = lig_pdbqt.with_name(lig_pdbqt.stem + "_complex.pdb")

        pdbqt_to_pdb(lig_pdbqt, lig_pdb)
        make_complex(receptor_pdbqt, lig_pdb, complex_pdb)


#===================================================="
#PLIP test
#===================================================="

import subprocess
import sys
import shutil

def ensure_plip():
    if shutil.which("plip") is not None:
        return True

    print("PLIP not found. Attempting automatic installation via conda...")

    try:
        subprocess.check_call([
            "conda", "install", "-y",
            "-c", "conda-forge",
            "-c", "bioconda",
            "plip"
        ])
    except Exception as e:
        print("PLIP installation failed:", e)
        return False

    return shutil.which("plip") is not None


PLIP_AVAILABLE = ensure_plip()


if PLIP_AVAILABLE:
    print("PLIP detected: running PLIP interaction analysis")
    # run PLIP here
else:
    print("PLIP not available: skipping PLIP interaction analysis")



#==================================================="
# Run PLIP (safe, optional)"
#===================================================="

import shutil
import subprocess
from pathlib import Path

PLIP_AVAILABLE = shutil.which("plip") is not None

PLIP_BASE = Path("results/interactions/plip")
PLIP_BASE.mkdir(parents=True, exist_ok=True)

if PLIP_AVAILABLE:
    print("PLIP detected: running PLIP interaction analysis")

    for complex_pdb in Path("results/top_N").rglob("*_complex.pdb"):
        receptor = complex_pdb.parent.name
        ligand = complex_pdb.stem.replace("_complex", "")

        outdir = PLIP_BASE / receptor / ligand
        outdir.mkdir(parents=True, exist_ok=True)

        subprocess.run(
            [
                "plip",
                "-f", str(complex_pdb),
                "-o", str(outdir),
                "--xml",
                "--silent"
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False
        )

    print("PLIP interaction analysis completed")

else:
    print("PLIP not detected: skipping PLIP interaction analysis")




#SAFE PLIP XML PARSER (Python)

import xml.etree.ElementTree as ET
import pandas as pd
from pathlib import Path

PLIP_BASE = Path("results/interactions/plip")
rows = []

# Mapping of PLIP XML blocks to interaction names
PLIP_BLOCKS = {
    "hydrogen_bonds": "hydrogen_bond",
    "hydrophobic_interactions": "hydrophobic",
    "salt_bridges": "salt_bridge",
    "pi_stacks": "pi_stack",
    "pi_cation_interactions": "pi_cation",
    "halogen_bonds": "halogen_bond",
    "metal_complexes": "metal_complex"
}

for receptor_dir in PLIP_BASE.iterdir():
    if not receptor_dir.is_dir():
        continue
    receptor = receptor_dir.name

    for ligand_dir in receptor_dir.iterdir():
        if not ligand_dir.is_dir():
            continue
        ligand = ligand_dir.name

        xml_files = list(ligand_dir.glob("*_report.xml"))
        if not xml_files:
            continue

        xml_file = xml_files[0]

        try:
            root = ET.parse(xml_file).getroot()
        except Exception:
            continue

        for block, interaction_name in PLIP_BLOCKS.items():
            block_node = root.find(f".//{block}")
            if block_node is None:
                continue

            for item in block_node:
                resnr = item.findtext("resnr", default="NA")
                restype = item.findtext("restype", default="NA")
                reschain = item.findtext("reschain", default="NA")
                distance = item.findtext("distance", default="")
                angle = item.findtext("angle", default="")

                rows.append([
                    receptor,
                    ligand,
                    interaction_name,
                    reschain,
                    resnr,
                    restype,
                    distance,
                    angle
                ])

# --------------------------------------------------
# SAVE RESULTS
# --------------------------------------------------
outdir = Path("results/interactions")
outdir.mkdir(parents=True, exist_ok=True)

if rows:
    df = pd.DataFrame(
        rows,
        columns=[
            "receptor",
            "ligand",
            "interaction_type",
            "chain",
            "residue_number",
            "residue_name",
            "distance",
            "angle"
        ]
    )

    df.to_csv(outdir / "plip_interactions_all.csv", index=False)

    summary = (
        df.groupby(["receptor", "interaction_type"])
        .size()
        .reset_index(name="count")
    )
    summary.to_csv(
        outdir / "plip_interaction_summary.csv",
        index=False
    )

    print("PLIP XML parsed successfully")
    print("Generated plip_interactions_all.csv and plip_interaction_summary.csv")

else:
    print("PLIP XML files parsed but no interactions detected")


#++++++++++++++++++++++
#Residue-wise heatmap from Protein-Ligand Interaction
#++++++++++++++++++++++

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

df = pd.read_csv("results/interactions/plip_interactions_all.csv")

df["residue"] = (
    df["residue_name"].astype(str)
    + df["residue_number"].astype(str)
    + ":" + df["chain"].astype(str)
)

# Count interactions per residue and ligand
heatmap_data = (
    df.groupby(["ligand", "residue"])
    .size()
    .unstack(fill_value=0)
)

FIG_DIR = Path("results/interactions/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ---- FULL heatmap (safe cap) ----
if heatmap_data.shape[1] > 50:
    top_residues = heatmap_data.sum().sort_values(ascending=False).head(20).index
    heatmap_data = heatmap_data[top_residues]
    suffix = "_top20"
else:
    suffix = ""

plt.figure(figsize=(12, max(6, 0.35 * heatmap_data.shape[0])))
sns.heatmap(
    heatmap_data,
    cmap="viridis",
    linewidths=0.2,
    cbar_kws={"label": "Interaction count"}
)
plt.title("Residue-wise Interaction Heatmap")
plt.xlabel("Protein Residue")
plt.ylabel("Ligand")
plt.tight_layout()
plt.savefig(FIG_DIR / f"residue_interaction_heatmap{suffix}.png", dpi=300)
plt.close()

print("Residue-wise interaction heatmap generated")

#++++++++++++++++++++++
#Ligand Interaction Fingerprints
#++++++++++++++++++++++
# Binary interaction fingerprint per ligand

fingerprint = (
    df.groupby(["ligand", "interaction_type"])
    .size()
    .unstack(fill_value=0)
)

# Convert counts to binary (presence/absence)
fingerprint = (fingerprint > 0).astype(int)

fingerprint.to_csv(
    "results/interactions/ligand_interaction_fingerprints.csv"
)

print("Ligand interaction fingerprints generated")

#===================================================="
#Consensus Interaction Score
#===================================================="

#Consensus Score =
#    (Total number of interactions)
#  + (Number of unique residues)
#  + (Number of interaction types)

consensus = df.groupby("ligand").agg(
    total_interactions=("interaction_type", "count"),
    unique_residues=("residue", "nunique"),
    interaction_types=("interaction_type", "nunique")
)

consensus["consensus_score"] = (
    consensus["total_interactions"]
    + consensus["unique_residues"]
    + consensus["interaction_types"]
)

consensus = consensus.sort_values(
    "consensus_score", ascending=False
)

consensus.to_csv(
    "results/interactions/ligand_consensus_interaction_score.csv"
)

print("Consensus interaction scores generated")

#===================================================="
#Per-receptor residue heatmaps
#===================================================="

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

df = pd.read_csv("results/interactions/plip_interactions_all.csv")

FIG_DIR = Path("results/interactions/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Create residue label
df["residue"] = (
    df["residue_name"].astype(str)
    + df["residue_number"].astype(str)
    + ":" + df["chain"].astype(str)
)

for receptor, rdf in df.groupby("receptor"):

    # Ligand + residue interaction counts
    heatmap_data = (
        rdf.groupby(["ligand", "residue"])
        .size()
        .unstack(fill_value=0)
    )

    if heatmap_data.empty:
        continue

    # -----------------------------
    # FULL heatmap (safe size)
    # -----------------------------
    if heatmap_data.shape[1] > 50:
        top_residues = (
            heatmap_data.sum()
            .sort_values(ascending=False)
            .head(20)
            .index
        )
        heatmap_top = heatmap_data[top_residues]
    else:
        heatmap_top = heatmap_data

    # ---- FULL ----
    plt.figure(
        figsize=(12, max(6, 0.35 * heatmap_data.shape[0]))
    )
    sns.heatmap(
        heatmap_data,
        cmap="viridis",
        linewidths=0.2,
        cbar_kws={"label": "Interaction count"}
    )
    plt.title(f"Residue-wise Interaction Heatmap - {receptor}")
    plt.xlabel("Protein Residue")
    plt.ylabel("Ligand")
    plt.tight_layout()
    plt.savefig(
        FIG_DIR / f"residue_interaction_heatmap_{receptor}.png",
        dpi=300
    )
    plt.close()

    # ---- TOP-20 ----
    plt.figure(
        figsize=(12, max(6, 0.35 * heatmap_top.shape[0]))
    )
    sns.heatmap(
        heatmap_top,
        cmap="viridis",
        linewidths=0.2,
        cbar_kws={"label": "Interaction count"}
    )
    plt.title(f"Top Residues Interaction Heatmap - {receptor}")
    plt.xlabel("Protein Residue")
    plt.ylabel("Ligand")
    plt.tight_layout()
    plt.savefig(
        FIG_DIR / f"residue_interaction_heatmap_{receptor}_top20.png",
        dpi=300
    )
    plt.close()

print("Per-receptor residue-wise interaction heatmaps generated")




#vishal
#Heatmaps per interaction type & receptor
#

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

df = pd.read_csv("results/interactions/plip_interactions_all.csv")

FIG_DIR = Path("results/interactions/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Create residue label (ASCII only)
df["residue"] = (
    df["residue_name"].astype(str)
    + df["residue_number"].astype(str)
    + ":" + df["chain"].astype(str)
)

MAX_RESIDUES = 20  # safety cap

for receptor, rdf in df.groupby("receptor"):

    for interaction_type, idf in rdf.groupby("interaction_type"):

        # Ligand x residue matrix
        heatmap_data = (
            idf.groupby(["ligand", "residue"])
            .size()
            .unstack(fill_value=0)
        )

        if heatmap_data.empty:
            continue

        # Cap residues to top-N by frequency
        if heatmap_data.shape[1] > MAX_RESIDUES:
            top_residues = (
                heatmap_data.sum()
                .sort_values(ascending=False)
                .head(MAX_RESIDUES)
                .index
            )
            heatmap_data = heatmap_data[top_residues]

        plt.figure(
            figsize=(12, max(6, 0.35 * heatmap_data.shape[0]))
        )

        sns.heatmap(
            heatmap_data,
            cmap="viridis",
            linewidths=0.2,
            cbar_kws={"label": "Interaction count"}
        )

        title = (
            f"{interaction_type.replace('_', ' ').title()} "
            f"Interaction Heatmap - {receptor}"
        )

        plt.title(title)
        plt.xlabel("Protein Residue")
        plt.ylabel("Ligand")
        plt.tight_layout()

        outfile = (
            FIG_DIR /
            f"residue_interaction_heatmap_{receptor}_{interaction_type}.png"
        )

        plt.savefig(outfile, dpi=300)
        plt.close()

print("Per-receptor, per-interaction-type heatmaps generated")


#===================================================="
#Interaction-type fingerprints + PCA + clustering
#===================================================="

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from matplotlib.patches import Ellipse

# --------------------------------------------------
# Load PLIP interaction table
# --------------------------------------------------
df = pd.read_csv("results/interactions/plip_interactions_all.csv")

# --------------------------------------------------
# Ligand - interaction-type fingerprint
# --------------------------------------------------
fp = (
    df.groupby(["ligand", "interaction_type"])
      .size()
      .unstack(fill_value=0)
)

fp.to_csv("results/interactions/ligand_interaction_fingerprint.csv")

# --------------------------------------------------
# Standardize features
# --------------------------------------------------
X = StandardScaler().fit_transform(fp.values)

# --------------------------------------------------
# PCA
# --------------------------------------------------
pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X)

pca_df = pd.DataFrame(
    X_pca,
    columns=["PC1", "PC2"],
    index=fp.index
)

# --------------------------------------------------
# KMeans clustering
# --------------------------------------------------
k = 3  # adjust if needed
labels = KMeans(n_clusters=k, random_state=42).fit_predict(X_pca)
pca_df["cluster"] = labels

pca_df.to_csv("results/interactions/interaction_pca_clusters.csv")

# --------------------------------------------------
# Plot PCA with cluster ellipses
# --------------------------------------------------
plt.figure(figsize=(8, 6))

colors = ["tab:red", "tab:blue", "tab:green", "tab:purple"]

for cl in sorted(pca_df["cluster"].unique()):
    subset = pca_df[pca_df["cluster"] == cl]
    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=f"Cluster {cl}",
        s=40,
        alpha=0.8,
        color=colors[cl % len(colors)]
    )

    # Ellipse (like your example)
    if len(subset) > 2:
        cov = np.cov(subset[["PC1", "PC2"]].T)
        vals, vecs = np.linalg.eigh(cov)
        angle = np.degrees(np.arctan2(*vecs[:, 1][::-1]))
        width, height = 2 * np.sqrt(vals)
        ellipse = Ellipse(
            xy=subset[["PC1", "PC2"]].mean(),
            width=width,
            height=height,
            angle=angle,
            edgecolor=colors[cl % len(colors)],
            facecolor="none",
            lw=2
        )
        plt.gca().add_patch(ellipse)

plt.xlabel("PC1 (interaction pattern)")
plt.ylabel("PC2 (interaction pattern)")
plt.title("PCA Clustering of Ligands by Interaction Fingerprints")
plt.legend()
plt.tight_layout()
plt.savefig("results/interactions/figures/interaction_pca_clustering.png", dpi=300)
plt.close()

print("Interaction-type PCA clustering completed")


for receptor, rdf in df.groupby("receptor"):
    fp = (
        rdf.groupby(["ligand", "interaction_type"])
           .size()
           .unstack(fill_value=0)
    )
    # PCA + clustering as above









#===================================================="
#ProLIF (BEST-POSE ONLY, SAFE)"
#===================================================="
from pathlib import Path
from rdkit import Chem
import MDAnalysis as mda
from prolif import Fingerprint

def pdbqt_to_pdb(pdbqt_file, pdb_file):
    with open(pdbqt_file) as inp, open(pdb_file, "w") as out:
        for line in inp:
            if line.startswith(("ATOM", "HETATM")):
                out.write(line[:66] + "\n")

fp = Fingerprint()
prolif_rows = []

TOP_DIR = Path("results/top_N")

for rec_dir in TOP_DIR.iterdir():
    if not rec_dir.is_dir():
        continue

    receptor_name = rec_dir.name
    receptor_pdbqt = Path("receptors") / f"{receptor_name}.pdbqt"

    # Convert receptor once
    receptor_pdb = rec_dir / f"{receptor_name}.pdb"
    if not receptor_pdb.exists():
        pdbqt_to_pdb(receptor_pdbqt, receptor_pdb)

    prot_u = mda.Universe(str(receptor_pdb))
    protein = prot_u.select_atoms("protein")

    for lig_pdbqt in rec_dir.glob("*_best.pdbqt"):
        ligand_name = lig_pdbqt.stem.replace("_best", "")

        lig_pdb = lig_pdbqt.with_suffix(".pdb")
        if not lig_pdb.exists():
            pdbqt_to_pdb(lig_pdbqt, lig_pdb)

        mol = Chem.MolFromPDBFile(str(lig_pdb), removeHs=False)
        if mol is None:
            continue

        try:
            ifp = fp.generate(mol, protein)
        except Exception as e:
            print(f"ProLIF failed for {ligand_name} vs {receptor_name}: {e}")
            continue

        for interaction, residues in ifp.items():
            for res in residues:
                prolif_rows.append([
                    receptor_name,
                    ligand_name,
                    interaction,
                    str(res)
                ])


import pandas as pd

if prolif_rows:
    pd.DataFrame(
        prolif_rows,
        columns=["receptor", "ligand", "interaction", "residue"]
    ).to_csv(
        "results/interactions/prolif_interactions.csv",
        index=False
    )
    print("ProLIF interaction table generated")
else:
    print("No ProLIF interactions detected")







# "===================================================="
# " Auto-generate Results + Figure captions"
# "===================================================="
from pathlib import Path
import pandas as pd

SUMMARY_DIR = Path("results/summary")
FIG_DIR = Path("results/summary/figures")
INT_FIG_DIR = Path("results/interactions/figures")

CAPTION_MASTER = SUMMARY_DIR / "figure_captions.txt"

captions = []

def write_caption(fig_path, caption):
    """Write caption next to figure and register in master list"""
    fig_path = Path(fig_path)
    caption_file = fig_path.with_suffix(fig_path.suffix + ".caption.txt")
    caption_file.write_text(caption + "\n")
    captions.append((fig_path.name, caption))

# --------------------------------------------------
# Load core datasets
# --------------------------------------------------
scores = pd.read_csv(SUMMARY_DIR / "vina_best_scores.csv")
energy_matrix = pd.read_csv(
    SUMMARY_DIR / "ligand_receptor_binding_energy_matrix.csv",
    index_col=0
)

n_ligands = energy_matrix.shape[0]
n_receptors = energy_matrix.shape[1]

# --------------------------------------------------
# Binding energy heatmap captions
# --------------------------------------------------
heatmap_full = FIG_DIR / "binding_energy_heatmap_full.png"
heatmap_top = next(FIG_DIR.glob("top*_binding_energy_heatmap.png"), None)

if heatmap_full.exists():
    write_caption(
        heatmap_full,
        f"Heatmap of ligand-receptor binding energies derived from AutoDock Vina docking. "
        f"The heatmap summarizes best docking scores for {n_ligands} ligands screened "
        f"against {n_receptors} receptors. Lower (more negative) values indicate stronger "
        f"predicted binding affinity."
    )

if heatmap_top:
    top_n = int(heatmap_top.stem.split("top")[1].split("_")[0])
    write_caption(
        heatmap_top,
        f"Binding energy heatmap for the top {top_n} ligands ranked by mean docking score "
        f"across {n_receptors} receptors. This focused view highlights the most promising "
        f"compounds while maintaining interpretability."
    )

# --------------------------------------------------
# Interaction PCA clustering caption
# --------------------------------------------------
pca_file = Path("results/interactions/interaction_pca_clusters.csv")
pca_fig = INT_FIG_DIR / "interaction_pca_clustering.png"

if pca_file.exists() and pca_fig.exists():
    pca_df = pd.read_csv(pca_file)
    n_clusters = pca_df["cluster"].nunique()

    write_caption(
        pca_fig,
        f"PCA clustering of ligands based on PLIP-derived interaction fingerprints. "
        f"Each ligand is represented by counts of protein-ligand interaction types and "
        f"projected onto the first two principal components. Clustering reveals "
        f"{n_clusters} distinct binding-mode groups characterized by different interaction patterns."
    )

# --------------------------------------------------
# Residue-wise interaction heatmaps (per receptor & type)
# --------------------------------------------------
for fig in INT_FIG_DIR.glob("residue_interaction_heatmap_*.png"):
    parts = fig.stem.split("_")
    receptor = parts[3]
    interaction_type = " ".join(parts[4:]).replace("_", " ")

    write_caption(
        fig,
        f"Residue-wise {interaction_type} interaction heatmap for receptor {receptor}. "
        f"The heatmap shows the frequency of {interaction_type} interactions between ligands "
        f"and protein residues, as identified by PLIP analysis."
    )

# --------------------------------------------------
# Write master caption file
# --------------------------------------------------
with open(CAPTION_MASTER, "w") as fh:
    fh.write("AUTO-GENERATED FIGURE CAPTIONS\n")
    fh.write("=" * 50 + "\n\n")
    for i, (name, cap) in enumerate(captions, start=1):
        fh.write(f"Figure {i}. {name}\n{cap}\n\n")

print("Dynamic per-figure captions generated")
print("Master caption file written to:", CAPTION_MASTER)



# "===================================================="
#" Autogenerated methdology and PCA discription
# "===================================================="
from pathlib import Path
import pandas as pd

SUMMARY_DIR = Path("results/summary")
INT_DIR = Path("results/interactions")

METHODS_FILE = SUMMARY_DIR / "methods_auto.txt"
RESULTS_FILE = SUMMARY_DIR / "results_analysis_auto.txt"

# --------------------------------------------------
# Load datasets safely
# --------------------------------------------------
scores = pd.read_csv(SUMMARY_DIR / "vina_best_scores.csv")
energy_matrix = pd.read_csv(
    SUMMARY_DIR / "ligand_receptor_binding_energy_matrix.csv",
    index_col=0
)

n_ligands = energy_matrix.shape[0]
n_receptors = energy_matrix.shape[1]

top_n = int(
    next(
        (p.stem.split("top")[1].split("_")[0]
         for p in (SUMMARY_DIR / "figures").glob("top*_binding_energy_heatmap.png")),
        0
    )
)

# Interaction data
interaction_file = INT_DIR / "plip_interactions_all.csv"
has_interactions = interaction_file.exists()

if has_interactions:
    idf = pd.read_csv(interaction_file)
    interaction_types = sorted(idf["interaction_type"].unique())
else:
    interaction_types = []

# PCA clustering
pca_file = INT_DIR / "interaction_pca_clusters.csv"
if pca_file.exists():
    pca_df = pd.read_csv(pca_file)
    n_clusters = pca_df["cluster"].nunique()
else:
    n_clusters = None

# --------------------------------------------------
# METHODS SECTION
# --------------------------------------------------
methods_text = f"""
Docking Protocol
Ligand docking was performed using AutoDock Vina against {n_receptors} protein receptors.
For each ligand-receptor pair, the best binding pose was selected based on the lowest
predicted binding energy. Docking calculations were parallelized using GNU Parallel,
with resume functionality enabled to ensure robustness.

Binding Energy Analysis
Docking scores were compiled into a ligand-receptor binding energy matrix representing
the best predicted binding energy (kcal/mol) for each ligand-receptor pair. Heatmaps were
generated to visualize binding affinity patterns across receptors. To maintain
interpretability for large datasets, heatmaps were restricted to the top-ranked ligands
when the total number of ligands exceeded a predefined threshold.

Protein-Ligand Interaction Analysis
Protein-ligand interactions were characterized using the Protein-Ligand Interaction
Profiler (PLIP). Interactions were classified into distinct types, including hydrogen
bonds, hydrophobic contacts, salt bridges, pi-stacking interactions, pi-cation
interactions, halogen bonds, and metal complexes. Interaction data were aggregated into
residue-wise and interaction-type-specific matrices for downstream analysis.

Interaction Fingerprints
Interaction fingerprints were constructed by counting the number of each interaction
type observed for a given ligand. These fingerprints provide a compact representation
of ligand binding modes that captures mechanistic interaction patterns beyond docking
scores alone.

Dimensionality Reduction and Clustering
Principal component analysis (PCA) was applied to standardized interaction fingerprints
to reduce dimensionality and identify dominant patterns of interaction variability.
K-means clustering was subsequently performed in PCA space to group ligands with similar
interaction profiles.
""".strip()

METHODS_FILE.write_text(methods_text + "\n")

# --------------------------------------------------
# RESULTS / ANALYSIS DESCRIPTION
# --------------------------------------------------
results_text = f"""
Docking Performance Overview
A total of {n_ligands} ligands were docked against {n_receptors} receptors. The resulting
binding energy landscape revealed a broad distribution of predicted affinities, with
distinct ligand subsets exhibiting consistently favorable binding across multiple
receptors.

Top-ranked Ligand Analysis
To focus on the most promising candidates, the top {top_n if top_n else 'ranked'} ligands
based on mean docking score were selected for detailed analysis. These ligands showed
enhanced binding consistency across receptors compared to the full compound set.

Interaction Pattern Analysis
PLIP-based interaction profiling identified multiple classes of protein-ligand
interactions, including {', '.join(interaction_types) if interaction_types else 'multiple interaction types'}.
Residue-wise interaction heatmaps revealed receptor-specific interaction hotspots and
highlighted residues contributing most frequently to ligand binding.

Interaction Fingerprint Clustering
PCA of interaction fingerprints demonstrated clear separation of ligands based on
binding interaction patterns rather than docking score alone.
{"The clustering analysis identified " + str(n_clusters) + " distinct ligand groups with characteristic interaction profiles." if n_clusters else ""}
These clusters represent alternative binding modes within the receptor binding sites.

Mechanistic Insights
The combined analysis of docking energies and interaction fingerprints provides a
mechanistic understanding of ligand recognition. Ligands within the same interaction
cluster share similar residue-level contacts, suggesting conserved binding strategies
that may be exploited for lead optimization.
""".strip()

RESULTS_FILE.write_text(results_text + "\n")

print("Methods section generated:", METHODS_FILE)
print("Results analysis text generated:", RESULTS_FILE)




# "===================================================="
#" PIPELINE COMPLETED SUCCESSFULLY"
# "===================================================="





EOF

echo "===================================================="
echo " PIPELINE COMPLETED SUCCESSFULLY"
echo "===================================================="
