# -*- coding: utf-8 -*-

# ====================================================
# 1. AUTO INSTALL (SAFE)
# ====================================================
def ensure_pkg(pkg):
    try:
        return __import__(pkg)
    except ImportError:
        import subprocess, sys
        print(f"[INFO] Installing {pkg}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
        return __import__(pkg)

nx = ensure_pkg("networkx")
ensure_pkg("pyvis")
from pyvis.network import Network

# ====================================================
# 2. IMPORTS
# ====================================================
import sys
import subprocess
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import roc_curve, auc

# ====================================================
# 3. CONFIG
# ====================================================
BASE_DIR = Path.cwd()
ALL_COMPLEX_DIR = Path("results/all_complexes")
RECEPTOR_DIR = Path("receptors")
SUMMARY_DIR = Path("results/summary")
FINAL_DIR = Path("results/final_ranking_per_receptor")

INT_FILE = ALL_COMPLEX_DIR / "plip_interactions_all.csv"
GNINA_FILE = FINAL_DIR / "gnina_scores_full.csv"

FINAL_DIR.mkdir(parents=True, exist_ok=True)

# ====================================================
# 4. GNINA MODULE
# ====================================================
def find_gnina():
    for p in [BASE_DIR/"gnina", BASE_DIR/"gnina.exe"]:
        if p.exists():
            return str(p.resolve())
    return shutil.which("gnina")

GNINA_BIN = find_gnina()
print(f"[INFO] GNINA: {GNINA_BIN}")

def gnina_task(args):
    receptor, rec_file, lig_file = args
    ligand = lig_file.stem.replace("_best", "")

    try:
        result = subprocess.run(
            [GNINA_BIN, "-r", str(rec_file), "-l", str(lig_file),
             "--score_only","--cnn","crossdock_default2018",
             "--cnn_scoring","rescore","--cpu","1"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )

        out = result.stdout + "\n" + result.stderr
        affinity = None
        cnn = None

        for line in out.splitlines():
            if line.startswith("Affinity:"):
                affinity = float(line.split()[1])
            elif line.startswith("CNNscore:"):
                cnn = float(line.split()[1])

        if cnn is not None:
            return (receptor, ligand, affinity, cnn)

    except Exception as e:
        print(f"[GNINA ERROR] {ligand}: {e}")

    return None

def generate_gnina():
    print("[STEP] Running GNINA...")
    tasks = []

    for rec_dir in ALL_COMPLEX_DIR.iterdir():
        receptor = rec_dir.name
        rec_file = RECEPTOR_DIR / f"{receptor}.pdbqt"

        for lig in rec_dir.glob("*_best.pdbqt"):
            tasks.append((receptor, rec_file, lig))

    results = []
    with Pool(max(1, cpu_count()-1)) as pool:
        for r in pool.map(gnina_task, tasks):
            if r:
                results.append(r)

    df = pd.DataFrame(results, columns=["receptor","ligand","affinity","cnn_score"])
    df.to_csv(GNINA_FILE, index=False)
    return df

# ====================================================
# 5. LOAD DATA
# ====================================================
cnn_df = pd.read_csv(GNINA_FILE) if GNINA_FILE.exists() else generate_gnina()

dock_df = pd.concat([
    pd.read_csv(f) for f in SUMMARY_DIR.glob("**/vina_best_scores.csv")
])

int_df = pd.read_csv(INT_FILE) if INT_FILE.exists() else pd.DataFrame()

clean = lambda x: str(x).replace("_best","").replace(".pdbqt","")
dock_df["ligand"] = dock_df["ligand"].apply(clean)
cnn_df["ligand"] = cnn_df["ligand"].apply(clean)
if not int_df.empty:
    int_df["ligand"] = int_df["ligand"].apply(clean)

# ====================================================
# 6. ANALYSIS FUNCTIONS
# ====================================================
def interaction_matrix(receptor):
    if int_df.empty:
        print("[SKIP] No PLIP data")
        return

    rec = int_df[int_df["receptor"] == receptor]
    if rec.empty:
        return

    out = FINAL_DIR / f"{receptor}_interaction_matrix.xlsx"

    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for itype, sdf in rec.groupby("type"):
            sdf = sdf.copy()
            sdf["residue_id"] = sdf["residue"] + sdf["resnum"].astype(str)

            mat = sdf.groupby(["ligand","residue_id"]).size().unstack(fill_value=0)
            mat = mat.astype(bool).replace({True:"Yes", False:""})
            mat.to_excel(writer, sheet_name=itype[:30])

def residue_analysis(receptor):
    rec = int_df[int_df["receptor"] == receptor].copy()
    if rec.empty:
        return

    rec["residue_id"] = rec["residue"] + rec["resnum"].astype(str)

    freq = rec["residue_id"].value_counts().head(20)
    sns.barplot(x=freq.values, y=freq.index)
    plt.savefig(FINAL_DIR / f"{receptor}_residue_freq.png")
    plt.close()

    mat = rec.groupby(["ligand","residue_id"]).size().unstack(fill_value=0)
    mat.astype(bool).astype(int).to_csv(FINAL_DIR / f"{receptor}_ifp_matrix.csv")

def clustering_analysis(receptor):
    f = FINAL_DIR / f"{receptor}_ifp_matrix.csv"
    if not f.exists():
        return
    df = pd.read_csv(f, index_col=0)

    if len(df) < 3:
        return

    labels = AgglomerativeClustering(n_clusters=3).fit_predict(df)
    pd.DataFrame({"Ligand": df.index, "Cluster": labels}).to_csv(
        FINAL_DIR / f"{receptor}_clusters.csv", index=False)


def ligand_network(receptor, max_ligands=100, top_edges=300):

    f = FINAL_DIR / f"{receptor}_ifp_matrix.csv"
    if not f.exists():
        return

    df = pd.read_csv(f, index_col=0)

    if len(df) > max_ligands:
        df = df.head(max_ligands)

    lig = df.index.tolist()
    edges = []

    def jaccard(a,b):
        a,b = np.array(a), np.array(b)
        inter = np.sum((a & b))
        union = np.sum((a | b))
        return inter/union if union!=0 else 0

    for i in range(len(lig)):
        for j in range(i+1,len(lig)):
            sim = jaccard(df.iloc[i], df.iloc[j])
            edges.append((lig[i], lig[j], sim))

    # ?? LIMIT EDGES
    edges = sorted(edges, key=lambda x: x[2], reverse=True)[:top_edges]

    print(f"[NETWORK] Using {len(edges)} edges")

    import networkx as nx
    import matplotlib.pyplot as plt

    G = nx.Graph()
    for a,b,s in edges:
        G.add_edge(a,b,weight=s)

    plt.figure(figsize=(8,6))
    pos = nx.spring_layout(G, seed=42, iterations=30)

    nx.draw(G, pos, node_size=200, font_size=6)

    out_png = FINAL_DIR / f"{receptor}_network.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[NETWORK] DONE ? {out_png}")


def ligand_top10_network(receptor):

    print(f"[TOP10] Processing {receptor}")

    f = FINAL_DIR / f"{receptor}_ifp_matrix.csv"
    if not f.exists():
        print("[TOP10] No IFP file")
        return

    df = pd.read_csv(f, index_col=0)

    lig = df.index.tolist()
    edges = []

    def jaccard(a,b):
        a,b = np.array(a), np.array(b)
        inter = np.sum((a & b))
        union = np.sum((a | b))
        return inter/union if union!=0 else 0

    # ?? LIMIT comparisons early
    max_check = min(len(lig), 50)   # <-- BIG SPEED FIX

    for i in range(max_check):
        for j in range(i+1, max_check):
            sim = jaccard(df.iloc[i], df.iloc[j])
            edges.append((lig[i], lig[j], sim))

    if not edges:
        return

    edges = sorted(edges, key=lambda x: x[2], reverse=True)[:10]

    import networkx as nx
    import matplotlib.pyplot as plt

    G = nx.Graph()
    for a,b,s in edges:
        G.add_edge(a,b,weight=s)

    plt.figure(figsize=(6,5))
    pos = nx.spring_layout(G, seed=42, iterations=20)

    nx.draw(G, pos, node_size=300, font_size=7)

    out_png = FINAL_DIR / f"{receptor}_top10_network.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[TOP10] DONE ? {out_png}")


def rank_based_consensus(receptor):

    dock = dock_df[dock_df["receptor"] == receptor]
    cnn = cnn_df[cnn_df["receptor"] == receptor]
    inter = int_df[int_df["receptor"] == receptor]

    dock_score = dock.groupby("ligand")["energy"].mean()
    stability = dock.groupby("ligand")["energy"].std().fillna(0)

    cnn_score = cnn.groupby("ligand")["cnn_score"].mean() if not cnn.empty else pd.Series()
    interaction = inter.groupby("ligand").size() if not inter.empty else pd.Series()

    df = pd.DataFrame({
        "Docking": dock_score,
        "CNN": cnn_score,
        "Interaction": interaction,
        "Stability": stability
    }).fillna(0)

    df["Vina_Rank"] = df["Docking"].rank(method="min", ascending=True)
    df["CNN_Rank"] = df["CNN"].rank(method="min", ascending=False)
    df["Interaction_Rank"] = df["Interaction"].rank(method="min", ascending=False)
    df["Stability_Rank"] = df["Stability"].rank(method="min", ascending=True)

    df["Avg_Rank"] = df[
        ["Vina_Rank","CNN_Rank","Interaction_Rank","Stability_Rank"]
    ].mean(axis=1)

    df = df.sort_values("Avg_Rank")

    out = FINAL_DIR / f"{receptor}_rank_consensus.csv"
    df.to_csv(out)

    print(f"[CONSENSUS] Saved ? {out}")


# ====================================================
# 7. MAIN
# ====================================================
def assign_tier(x):
    if x >= 0.75: return "Tier 1"
    elif x >= 0.5: return "Tier 2"
    else: return "Tier 3"

for receptor, rdf in dock_df.groupby("receptor"):

    print(f"[PROCESS] {receptor}")

    dock_score = rdf.groupby("ligand")["energy"].mean()
    std = rdf.groupby("ligand")["energy"].std().fillna(0)

    dock_norm = 1 - (dock_score-dock_score.min())/(dock_score.max()-dock_score.min()+1e-6)
    stab_norm = 1 - (std-std.min())/(std.max()-std.min()+1e-6)

    cnn = cnn_df[cnn_df["receptor"]==receptor].groupby("ligand")["cnn_score"].mean()
    cnn_norm = (cnn-cnn.min())/(cnn.max()-cnn.min()+1e-6)

    combined = pd.DataFrame({
        "Docking": dock_norm,
        "CNN": cnn_norm,
        "Stability": stab_norm
    }).fillna(0)

    combined["Final_Score"] = (
        0.5*combined["Docking"] +
        0.3*combined["CNN"] +
        0.2*combined["Stability"]
    )

    combined["Tier"] = combined["Final_Score"].apply(assign_tier)
    combined.to_csv(FINAL_DIR / f"{receptor}_ranking.csv")

    # ---- ANALYSIS ----
    interaction_matrix(receptor)
    residue_analysis(receptor)
    clustering_analysis(receptor)

    # ---- CONSENSUS ----
    rank_based_consensus(receptor)

    # ---- NETWORKS ----
    ligand_network(receptor)
    ligand_top10_network(receptor)


print("\n[SUCCESS] PIPELINE COMPLETED (WITHOUT PCA)")

print("\n[SUCCESS] FULL PIPELINE COMPLETED")