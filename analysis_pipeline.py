# -*- coding: utf-8 -*-

# ====================================================
# USER CONFIGURATION
# ====================================================
TOP_N = 10  # default

# CLI override
import sys
if len(sys.argv) > 1:
    try:
        TOP_N = int(sys.argv[1])
        print(f"[INFO] Using TOP_N = {TOP_N}")
    except:
        print("[WARNING] Invalid input, using default 10")

# ====================================================
# PACKAGE CHECK
# ====================================================
import subprocess, importlib, shutil

REQUIRED_PACKAGES = ["pandas", "numpy", "matplotlib", "openpyxl", "reportlab"]

def install(pkg):
    try:
        subprocess.check_call(["conda", "install", "-y", pkg])
    except:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])

for p in REQUIRED_PACKAGES:
    try:
        importlib.import_module(p)
    except:
        print(f"[INSTALL] {p}")
        install(p)

# ====================================================
# IMPORTS
# ====================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

# ====================================================
# PATHS
# ====================================================
SUMMARY_DIR = Path("results/summary")
INT_DIR = Path("results/interactions")
FINAL_DIR = Path("results/final_ranking")
FIG_DIR = FINAL_DIR / "figures"

FINAL_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

print("\n[STEP] Loading data...")

# ====================================================
# LOAD MULTI-RECEPTOR DATA
# ====================================================
files = list(SUMMARY_DIR.glob("**/vina_best_scores.csv"))
df_list = []

for f in files:
    tmp = pd.read_csv(f)
    if "receptor" not in tmp.columns:
        tmp["receptor"] = f.parent.name
    df_list.append(tmp)

dock_df = pd.concat(df_list, ignore_index=True)

# ====================================================
# LOAD INTERACTIONS
# ====================================================
int_file = INT_DIR / "plip_interactions_all.csv"
int_df = pd.read_csv(int_file) if int_file.exists() else pd.DataFrame()

# ====================================================
# CLEAN NAMES
# ====================================================
def clean(x):
    return str(x).replace("_best","").replace(".pdbqt","").strip()

dock_df["ligand"] = dock_df["ligand"].apply(clean)
if not int_df.empty:
    int_df["ligand"] = int_df["ligand"].apply(clean)

# ====================================================
# SCORES
# ====================================================
dock_score = dock_df.groupby("ligand")["energy"].mean()
dock_norm = 1 - (dock_score - dock_score.min())/(dock_score.max()-dock_score.min())

std = dock_df.groupby("ligand")["energy"].std().fillna(0)
stab_norm = 1 - (std - std.min())/(std.max()-std.min()+1e-6)

if not int_df.empty:
    int_df["weight"] = 1
    inter = int_df.groupby("ligand")["weight"].sum()
    inter_norm = (inter-inter.min())/(inter.max()-inter.min()+1e-6)
else:
    inter_norm = pd.Series(0, index=dock_score.index)

combined = pd.DataFrame({
    "Docking": dock_norm,
    "Stability": stab_norm,
    "Interaction": inter_norm
}).fillna(0)

combined["Final_Score"] = (
    0.5*combined["Docking"] +
    0.3*combined["Interaction"] +
    0.2*combined["Stability"]
)

combined = combined.sort_values("Final_Score", ascending=False)

def classify(x):
    if x > 0.75: return "Tier 1 (High Confidence)"
    elif x > 0.5: return "Tier 2 (Moderate)"
    else: return "Tier 3 (Low)"

combined["Category"] = combined["Final_Score"].apply(classify)

combined.to_csv(FINAL_DIR/"final_integrated_ranking.csv")

topN = combined.head(TOP_N)

# ====================================================
# EXCEL
# ====================================================
wb = Workbook()
wb.remove(wb.active)

def add(ws_name, df):
    ws = wb.create_sheet(ws_name[:30])
    for r in dataframe_to_rows(df.reset_index(), index=False, header=True):
        ws.append(r)

add("All", combined)
add(f"Top_{TOP_N}", topN)

wb.save(FINAL_DIR/"final_report.xlsx")

# ====================================================
# CLEAN PUBLICATION-QUALITY PLOTS
# ====================================================

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9
})

tier_colors = {
    "Tier 1 (High Confidence)": "#d62728",
    "Tier 2 (Moderate)": "#1f77b4",
    "Tier 3 (Low)": "#7f7f7f"
}

# ====================================================
# FIGURE 1: SCATTER (FIXED)
# ====================================================
fig, ax = plt.subplots(figsize=(6,5))

for tier, color in tier_colors.items():
    sub = combined[combined["Category"] == tier]
    ax.scatter(sub["Docking"], sub["Interaction"],
               c=color, label=tier, alpha=0.7)

ax.scatter(topN["Docking"], topN["Interaction"],
           marker="*", s=140, edgecolor="black",
           label=f"Top {TOP_N}")

ax.set_xlabel("Normalized Docking Score")
ax.set_ylabel("Normalized Interaction Score")
ax.set_title("Figure 1. Docking vs Interaction")

ax.legend(loc="best", frameon=False)

plt.tight_layout()
plt.savefig(FIG_DIR / "Figure1.png", dpi=300)
plt.close()


# ====================================================
# FIGURE 2: HISTOGRAM (FIXED)
# ====================================================
fig, ax = plt.subplots(figsize=(6,5))

ax.hist(combined["Final_Score"], bins=30)

ax.set_xlabel("Final Score")
ax.set_ylabel("Frequency")
ax.set_title("Figure 2. Final Score Distribution")

plt.tight_layout()
plt.savefig(FIG_DIR / "Figure2.png", dpi=300)
plt.close()


# ====================================================
# FIGURE 3: TIER DISTRIBUTION (FIXED)
# ====================================================
tier_counts = combined["Category"].value_counts()

fig, ax = plt.subplots(figsize=(6,5))

ax.bar(tier_counts.index, tier_counts.values,
       color=[tier_colors[t] for t in tier_counts.index])

ax.set_ylabel("Number of Ligands")
ax.set_title("Figure 3. Ligand Distribution by Tier")

ax.set_xticks(range(len(tier_counts.index)))
ax.set_xticklabels(tier_counts.index, rotation=25, ha="right")

plt.tight_layout()
plt.savefig(FIG_DIR / "Figure3.png", dpi=300)
plt.close()


# ====================================================
# FIGURE 4: BOXPLOT (FIXED)
# ====================================================
fig, ax = plt.subplots(figsize=(6,5))

combined.boxplot(column="Final_Score", by="Category", ax=ax)

ax.set_title("Figure 4. Score Distribution Across Tiers")
ax.set_xlabel("Tier")
ax.set_ylabel("Final Score")

plt.suptitle("")
plt.xticks(rotation=25)

plt.tight_layout()
plt.savefig(FIG_DIR / "Figure4.png", dpi=300)
plt.close()


# ====================================================
# FIGURE 5: TOP N (FIXED)
# ====================================================
fig, ax = plt.subplots(figsize=(7,5))

labels = [str(lig)[:30] for lig in topN.index[::-1]]

ax.barh(range(len(topN)), topN["Final_Score"][::-1])

ax.set_yticks(range(len(topN)))
ax.set_yticklabels(labels)

ax.set_xlabel("Final Score")
ax.set_title(f"Figure 5. Top {TOP_N} Ligands")

plt.tight_layout()
plt.savefig(FIG_DIR / "Figure5.png", dpi=300)
plt.close()


# ====================================================
# FIGURE 6: HEATMAP (FIXED)
# ====================================================
pivot = None

if "receptor" in dock_df.columns:

    pivot = dock_df.pivot_table(
        index="ligand",
        columns="receptor",
        values="energy"
    )

    pivot = pivot.head(50)

    fig, ax = plt.subplots(figsize=(7,6))

    c = ax.imshow(pivot, aspect="auto")

    ax.set_xlabel("Receptor")
    ax.set_ylabel("Top Ligands")

    plt.colorbar(c, ax=ax)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "Figure6.png", dpi=300)
    plt.close()# ====================================================
# ====================================================
# FINAL PDF REPORT (CLEAN + ERROR FREE)
# ====================================================

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet

print("\n[STEP] Generating Final PDF Report...")

pdf_file = FINAL_DIR / "Docking_Report.pdf"
styles = getSampleStyleSheet()

doc = SimpleDocTemplate(str(pdf_file), pagesize=A4)

# ? IMPORTANT
story = []

# --------------------------------------------------
# TITLE
# --------------------------------------------------
story.append(Paragraph("Docking Analysis Report", styles["Title"]))
story.append(Spacer(1, 12))

# --------------------------------------------------
# TOP N TABLE
# --------------------------------------------------
table_data = [["Ligand", "Final Score", "Category"]]

for lig, row in topN.iterrows():
    table_data.append([
        str(lig),
        f"{row['Final_Score']:.3f}",
        row["Category"]
    ])

table = Table(table_data)
table.setStyle(TableStyle([
    ("BACKGROUND", (0,0), (-1,0), colors.grey),
    ("TEXTCOLOR",(0,0),(-1,0),colors.whitesmoke),
    ("GRID", (0,0), (-1,-1), 1, colors.black)
]))

story.append(Paragraph(f"Top {TOP_N} Ligands", styles["Heading2"]))
story.append(Spacer(1, 8))
story.append(table)
story.append(Spacer(1, 20))

# ====================================================
# DATA SUMMARY (NOW CORRECT)
# ====================================================
story.append(Paragraph("Data Summary", styles["Heading2"]))
story.append(Spacer(1, 10))

tier_counts = combined["Category"].value_counts()

for tier, count in tier_counts.items():
    story.append(Paragraph(f"{tier}: {count} ligands", styles["BodyText"]))

story.append(Spacer(1, 8))

story.append(Paragraph(
    f"Score Range: {combined['Final_Score'].min():.3f} - {combined['Final_Score'].max():.3f}",
    styles["BodyText"]
))

story.append(Spacer(1, 20))

# --------------------------------------------------
# CAPTIONS
# --------------------------------------------------
captions = {
    1: "Docking vs interaction scatter plot.",
    2: "Distribution of final scores.",
    3: "Ligand classification into tiers.",
    4: "Score distribution across tiers.",
    5: f"Top {TOP_N} ligands ranking.",
    6: "Binding energy heatmap across receptors."
}

# --------------------------------------------------
# FIGURES
# --------------------------------------------------
for i in range(1, 7):
    fig_path = FIG_DIR / f"Figure{i}.png"

    if fig_path.exists():
        story.append(PageBreak())

        story.append(Paragraph(f"Figure {i}", styles["Heading2"]))
        story.append(Spacer(1, 10))

        story.append(Image(str(fig_path), width=450, height=320))
        story.append(Spacer(1, 8))

        story.append(Paragraph(f"<i>{captions[i]}</i>", styles["BodyText"]))

# --------------------------------------------------
# BUILD PDF
# --------------------------------------------------
doc.build(story)

print(f"[SUCCESS] PDF generated -> {pdf_file}")
