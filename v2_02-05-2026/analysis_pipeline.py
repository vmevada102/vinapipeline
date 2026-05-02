# -*- coding: utf-8 -*-

# ====================================================
# AUTO-INSTALL networkx (SAFE)
# ====================================================
def ensure_networkx():
    try:
        import networkx as nx
        return nx
    except ImportError:
        print("[INFO] networkx not found ? installing...")

        import subprocess
        import sys

        subprocess.check_call([
            sys.executable, "-m", "pip", "install", "networkx"
        ])

        import networkx as nx
        print("[INFO] networkx installed successfully")
        return nx

# ?? use this instead of normal import
nx = ensure_networkx()



import sys
import subprocess
import shutil
from pathlib import Path
import pandas as pd
import networkx as nx
import numpy as np
from multiprocessing import Pool, cpu_count
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import roc_curve, auc

# ====================================================
# CONFIG
# ====================================================
TOP_N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
FORCE = "--force" in sys.argv

BASE_DIR = Path.cwd()

ALL_COMPLEX_DIR = Path("results/all_complexes")
RECEPTOR_DIR = Path("receptors")
SUMMARY_DIR = Path("results/summary")
FINAL_DIR = Path("results/final_ranking_per_receptor")
FIG_DIR = FINAL_DIR / "figures"

INT_FILE = ALL_COMPLEX_DIR / "plip_interactions_all.csv"
GNINA_FILE = FINAL_DIR / "gnina_scores_full.csv"

FINAL_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

# ====================================================
# GNINA DETECT
# ====================================================
def find_gnina():
    for p in [BASE_DIR/"gnina", BASE_DIR/"gnina.exe"]:
        if p.exists():
            return str(p.resolve())
    return shutil.which("gnina")

GNINA_BIN = find_gnina()
print(f"[INFO] GNINA: {GNINA_BIN}")

# ====================================================
# GNINA TASK
# ====================================================
def gnina_task(args):
    receptor, rec_file, lig_file = args
    ligand = lig_file.stem.replace("_best", "")

    try:
        result = subprocess.run(
            [
                GNINA_BIN,
                "-r", str(rec_file),
                "-l", str(lig_file),
                "--score_only",
                "--cnn", "crossdock_default2018",
                "--cnn_scoring", "rescore",
                "--cpu", "1"
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        output = (result.stdout or "") + "\n" + (result.stderr or "")

        affinity = None
        cnn_score = None

        for line in output.splitlines():
            if line.startswith("Affinity:"):
                affinity = float(line.split()[1])
            elif line.startswith("CNNscore:"):
                cnn_score = float(line.split()[1])

        if cnn_score is not None:
            return (receptor, ligand, affinity, cnn_score)

    except Exception as e:
        print(f"[ERROR] {ligand}: {e}")

    return None

# ====================================================
# GNINA PARALLEL
# ====================================================
def generate_gnina_scores():
    tasks = []

    for rec_dir in ALL_COMPLEX_DIR.iterdir():
        receptor = rec_dir.name
        rec_file = RECEPTOR_DIR / f"{receptor}.pdbqt"

        for lig_file in rec_dir.glob("*_best.pdbqt"):
            tasks.append((receptor, rec_file, lig_file))

    results = []

    with Pool(max(1, cpu_count()-1)) as pool:
        for r in pool.map(gnina_task, tasks):
            if r:
                results.append(r)

    df = pd.DataFrame(results, columns=[
        "receptor","ligand","affinity","cnn_score"
    ])

    df.to_csv(GNINA_FILE, index=False)
    return df

# ====================================================
# LOAD DATA
# ====================================================
cnn_df = pd.read_csv(GNINA_FILE) if GNINA_FILE.exists() and not FORCE else generate_gnina_scores()

dock_df = pd.concat(
    [pd.read_csv(f) for f in SUMMARY_DIR.glob("**/vina_best_scores.csv")],
    ignore_index=True
)

int_df = pd.read_csv(INT_FILE) if INT_FILE.exists() else pd.DataFrame()

# ====================================================
# CLEAN
# ====================================================
def clean(x):
    return str(x).replace("_best","").replace(".pdbqt","")

dock_df["ligand"] = dock_df["ligand"].apply(clean)
cnn_df["ligand"] = cnn_df["ligand"].apply(clean)

if not int_df.empty:
    int_df["ligand"] = int_df["ligand"].apply(clean)

# ====================================================
# FUNCTIONS
# ====================================================
def interaction_matrix(receptor):
    rec = int_df[int_df["receptor"] == receptor]
    if rec.empty:
        return

    with pd.ExcelWriter(FINAL_DIR / f"{receptor}_interaction_matrix.xlsx") as writer:
        for itype, sdf in rec.groupby("type"):
            sdf = sdf.copy()
            sdf["residue_id"] = sdf["residue"] + sdf["resnum"].astype(str)
            mat = sdf.groupby(["ligand","residue_id"]).size().unstack(fill_value=0)
            mat = mat.astype(bool).replace({True:"Yes", False:""})
            mat.to_excel(writer, sheet_name=itype[:30])

def residue_analysis(receptor):

    rec = int_df[int_df["receptor"] == receptor].copy()   # ? FIX

    if rec.empty:
        return

    rec["residue_id"] = rec["residue"] + rec["resnum"].astype(str)

    freq = rec["residue_id"].value_counts().head(20)

    plt.figure(figsize=(10,6))
    sns.barplot(x=freq.values, y=freq.index)
    plt.title(f"{receptor} Residue Frequency")
    plt.savefig(FINAL_DIR / f"{receptor}_residue_frequency.png")
    plt.close()

    total = rec["ligand"].nunique()
    hotspot = rec.groupby("residue_id")["ligand"].nunique() / total

    hotspot[hotspot > 0.3].to_csv(
        FINAL_DIR / f"{receptor}_hotspots.csv"
    )

    mat = rec.groupby(["ligand","residue_id"]).size().unstack(fill_value=0)
    mat.astype(bool).astype(int).to_csv(
        FINAL_DIR / f"{receptor}_ifp_matrix.csv"
    )

def pca_analysis(receptor):

    file = FINAL_DIR / f"{receptor}_ifp_matrix.csv"
    if not file.exists():
        return

    df = pd.read_csv(file, index_col=0)

    # ? FIX: skip bad matrices
    if df.shape[1] < 2 or df.var().sum() == 0:
        print(f"[SKIP PCA] {receptor} (no variance)")
        return

    coords = PCA(n_components=2).fit_transform(df)

    plt.figure()
    plt.scatter(coords[:,0], coords[:,1])
    plt.title(f"{receptor} PCA")
    plt.savefig(FINAL_DIR / f"{receptor}_pca.png")
    plt.close()


def clustering_analysis(receptor):
    file = FINAL_DIR / f"{receptor}_ifp_matrix.csv"
    if not file.exists():
        return
    df = pd.read_csv(file, index_col=0)

    if df.shape[0] < 3:
        return

    labels = AgglomerativeClustering(n_clusters=3).fit_predict(df)

    pd.DataFrame({"Ligand": df.index, "Cluster": labels}).to_csv(
        FINAL_DIR / f"{receptor}_clusters.csv", index=False
    )

def roc_analysis(receptor):
    file = FINAL_DIR / f"{receptor}_ranking.csv"
    act_file = FINAL_DIR / "actives.csv"

    if not file.exists() or not act_file.exists():
        return

    df = pd.read_csv(file)
    actives = pd.read_csv(act_file)["ligand"].astype(str).tolist()

    df["label"] = df.index.astype(str).apply(lambda x: 1 if x in actives else 0)

    if df["label"].sum() == 0:
        return

    fpr, tpr, _ = roc_curve(df["label"], df["Final_Score"])
    roc_auc = auc(fpr, tpr)

    plt.plot(fpr, tpr, label=f"AUC={roc_auc:.2f}")
    plt.plot([0,1],[0,1],'--')
    plt.legend()
    plt.savefig(FINAL_DIR / f"{receptor}_roc.png")
    plt.close()

def consensus(receptor):

    dock = dock_df[dock_df["receptor"]==receptor]
    cnn = cnn_df[cnn_df["receptor"]==receptor]
    inter = int_df[int_df["receptor"]==receptor]

    dock_score = dock.groupby("ligand")["energy"].mean()
    stability = dock.groupby("ligand")["energy"].std().fillna(0)
    cnn_score = cnn.groupby("ligand")["cnn_score"].mean()
    interaction = inter.groupby("ligand").size()

    df = pd.DataFrame({
        "Docking": dock_score,
        "CNN": cnn_score,
        "Interaction": interaction,
        "Stability": stability
    }).fillna(0)

    # ? CORRECT ranking
    df["Vina_Rank"] = df["Docking"].rank(ascending=True)
    df["CNN_Rank"] = df["CNN"].rank(ascending=False)
    df["Interaction_Rank"] = df["Interaction"].rank(ascending=False)
    df["Stability_Rank"] = df["Stability"].rank(ascending=True)

    df["Avg_Rank"] = df[
        ["Vina_Rank","CNN_Rank","Interaction_Rank","Stability_Rank"]
    ].mean(axis=1)

    df.sort_values("Avg_Rank").to_csv(
        FINAL_DIR / f"{receptor}_consensus_ranking.csv"
    )

# ====================================================
# LIGAND SIMILARITY NETWORK
# ====================================================
def ligand_similarity_network(receptor, threshold=0.3):

    file = FINAL_DIR / f"{receptor}_ifp_matrix.csv"
    if not file.exists():
        print(f"[SKIP] No IFP file for {receptor}")
        return

    df = pd.read_csv(file, index_col=0)

    if df.shape[0] < 2:
        print(f"[SKIP] Not enough ligands for network ({receptor})")
        return

    import numpy as np

    def jaccard(a, b):
        a = np.array(a)
        b = np.array(b)
        inter = np.sum((a & b))
        union = np.sum((a | b))
        return inter / union if union != 0 else 0

    ligands = df.index.tolist()
    edges = []

    for i in range(len(ligands)):
        for j in range(i + 1, len(ligands)):
            sim = jaccard(df.iloc[i], df.iloc[j])
            if sim >= threshold:
                edges.append((ligands[i], ligands[j], sim))

    # Save edges
    edge_df = pd.DataFrame(edges, columns=["Ligand1", "Ligand2", "Similarity"])
    edge_df.to_csv(FINAL_DIR / f"{receptor}_network_edges.csv", index=False)

    # Build graph
    G = nx.Graph()

    for lig in ligands:
        G.add_node(lig)

    for l1, l2, sim in edges:
        G.add_edge(l1, l2, weight=sim)

    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G, seed=42)

    nx.draw(
        G,
        pos,
        with_labels=True,
        node_size=300,
        font_size=6
    )

    plt.title(f"{receptor} Ligand Similarity Network")
    plt.savefig(FINAL_DIR / f"{receptor}_network.png", dpi=300)
    plt.close()




# ====================================================
# MAIN
# ====================================================

# ? define function OUTSIDE loop
def assign_tier(score):
    if score >= 0.75:
        return "Tier 1 (Good)"
    elif score >= 0.5:
        return "Tier 2 (Moderate)"
    else:
        return "Tier 3 (Poor)"


for receptor, rdf in dock_df.groupby("receptor"):

    print(f"[PROCESS] {receptor}")

    # ---- Docking ----
    dock_score = rdf.groupby("ligand")["energy"].mean()
    std = rdf.groupby("ligand")["energy"].std().fillna(0)

    dock_norm = 1 - (dock_score - dock_score.min())/(dock_score.max()-dock_score.min()+1e-6)
    stab_norm = 1 - (std - std.min())/(std.max()-std.min()+1e-6)

    # ---- CNN ----
    cnn = cnn_df[cnn_df["receptor"]==receptor].groupby("ligand")["cnn_score"].mean()
    cnn_norm = (cnn - cnn.min())/(cnn.max()-cnn.min()+1e-6)

    # ---- Interaction ----
    inter = int_df[int_df["receptor"]==receptor].groupby("ligand").size()
    inter_norm = (inter - inter.min())/(inter.max()-inter.min()+1e-6)

    # ---- Combine ----
    combined = pd.DataFrame({
        "Docking": dock_norm,
        "CNN": cnn_norm,
        "Interaction": inter_norm,
        "Stability": stab_norm
    }).fillna(0)

    combined["Final_Score"] = (
        0.4*combined["Docking"] +
        0.25*combined["CNN"] +
        0.2*combined["Interaction"] +
        0.15*combined["Stability"]
    )

    # ? Tier INSIDE loop
    combined["Tier"] = combined["Final_Score"].apply(assign_tier)

    # ---- Save ----
    combined.to_csv(FINAL_DIR / f"{receptor}_ranking.csv")

    # ---- Analyses ----
    interaction_matrix(receptor)
    residue_analysis(receptor)
    pca_analysis(receptor)
    clustering_analysis(receptor)
    roc_analysis(receptor)
    consensus(receptor)

    # ---- Network ----
    ligand_similarity_network(receptor)





print("\n[SUCCESS] FULL PIPELINE WITH ALL ANALYTICS COMPLETED")

