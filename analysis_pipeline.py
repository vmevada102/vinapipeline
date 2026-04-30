# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
from pathlib import Path

print("\n[STEP] Running Analysis Pipeline...")

# ====================================================
# PATHS
# ====================================================
SUMMARY_DIR = Path("results/summary")
INT_DIR = Path("results/interactions")
FINAL_DIR = Path("results/final_ranking")
FINAL_DIR.mkdir(parents=True, exist_ok=True)

# ====================================================
# LOAD DOCKING DATA
# ====================================================
dock_df = pd.read_csv(SUMMARY_DIR / "vina_best_scores.csv")

# ====================================================
# LOAD INTERACTION DATA (SAFE)
# ====================================================
interaction_file = INT_DIR / "plip_interactions_all.csv"

if interaction_file.exists():
    int_df = pd.read_csv(interaction_file)
else:
    print("WARNING: No PLIP data found. Using empty interactions.")
    int_df = pd.DataFrame(columns=["ligand", "interaction_type"])

# ====================================================
# DOCKING SCORE
# ====================================================
dock_score = dock_df.groupby("ligand")["energy"].mean()

dock_norm = (dock_score - dock_score.min()) / (dock_score.max() - dock_score.min())
dock_norm = 1 - dock_norm

# ====================================================
# STABILITY SCORE
# ====================================================
dock_std = dock_df.groupby("ligand")["energy"].std().fillna(0)

stab_norm = (dock_std - dock_std.min()) / (dock_std.max() - dock_std.min() + 1e-6)
stab_norm = 1 - stab_norm

# ====================================================
# INTERACTION SCORE
# ====================================================
if not int_df.empty:
    weights = {
        "hydrogen_bond": 3,
        "salt_bridge": 2,
        "pi_stack": 2,
        "pi_cation": 2,
        "hydrophobic": 1,
        "halogen_bond": 2,
        "metal_complex": 3
    }

    int_df["weight"] = int_df["interaction_type"].map(weights).fillna(1)
    interaction_score = int_df.groupby("ligand")["weight"].sum()

    int_norm = (interaction_score - interaction_score.min()) / (
        interaction_score.max() - interaction_score.min() + 1e-6
    )
else:
    int_norm = pd.Series(0, index=dock_score.index)

# ====================================================
# COMBINE SCORES
# ====================================================
combined = pd.DataFrame({
    "Docking": dock_norm,
    "Stability": stab_norm,
    "Interaction": int_norm
}).fillna(0)

combined["Final_Score"] = (
    0.5 * combined["Docking"] +
    0.3 * combined["Interaction"] +
    0.2 * combined["Stability"]
)

combined = combined.sort_values("Final_Score", ascending=False)
combined["Rank"] = range(1, len(combined) + 1)

# ====================================================
# CLASSIFICATION
# ====================================================
def classify(score):
    if score > 0.75:
        return "Tier 1 (High Confidence)"
    elif score > 0.5:
        return "Tier 2 (Moderate)"
    else:
        return "Tier 3 (Low)"

combined["Category"] = combined["Final_Score"].apply(classify)

# ====================================================
# SAVE FINAL RANKING
# ====================================================
out_file = FINAL_DIR / "final_integrated_ranking.csv"
combined.to_csv(out_file)

print("[SUCCESS] Final ranking saved ->", out_file)

# ====================================================
# TOP 10 SUMMARY
# ====================================================
top10 = combined.head(10)

with open(FINAL_DIR / "top10_summary.txt", "w") as f:
    f.write("Top 10 Ligands (Final Integrated Ranking)\n")
    f.write("========================================\n\n")

    for i, (lig, row) in enumerate(top10.iterrows(), 1):
        f.write(f"{i}. {lig} | Score: {row['Final_Score']:.3f} | {row['Category']}\n")

print("[SUCCESS] Top 10 summary generated")

# ====================================================
# INTERACTION SUMMARY
# ====================================================
if not int_df.empty:
    interaction_summary = (
        int_df.groupby(["ligand", "interaction_type"])
        .size()
        .unstack(fill_value=0)
    )
else:
    interaction_summary = pd.DataFrame()

# ====================================================
# DOCKING SUMMARY
# ====================================================
dock_summary = dock_df.groupby("ligand").agg(
    Mean_Energy=("energy", "mean"),
    Std_Energy=("energy", "std")
).fillna(0)

# ====================================================
# TOP COMPARISON TABLE
# ====================================================
top_compare = combined.head(10)
top_compare = top_compare.merge(dock_summary, left_index=True, right_index=True)

if not interaction_summary.empty:
    top_compare = top_compare.merge(interaction_summary, left_index=True, right_index=True)

top_compare.to_csv(FINAL_DIR / "top_ligand_comparison.csv")

print("[SUCCESS] Top ligand comparison table generated")

# ====================================================
# INTERPRETATION FILE
# ====================================================
with open(FINAL_DIR / "final_score_interpretation.txt", "w") as f:
    f.write("FINAL SCORE INTERPRETATION\n")
    f.write("=========================\n\n")

    f.write("Final score combines docking, interaction, and stability.\n\n")
    f.write("Score ranges:\n")
    f.write("  - >0.75  : Strong binding\n")
    f.write("  - 0.50-0.75 : Moderate binding\n")
    f.write("  - <0.50  : Weak binding\n\n")

    f.write(f"Best score: {combined['Final_Score'].max():.3f}\n")
    f.write(f"Mean score: {combined['Final_Score'].mean():.3f}\n")

print("[SUCCESS] Interpretation report generated")

print("\nPIPELINE ANALYSIS COMPLETED SUCCESSFULLY\n")
