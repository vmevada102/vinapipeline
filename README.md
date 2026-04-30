# 🧬 AutoDock Vina Master Pipeline  
**Parallel Docking · Consensus Ranking · Interaction Fingerprints · Publication-Ready Analytics**

---

## 📌 Overview

This repository provides a **fully automated, end-to-end molecular docking and post-docking analysis pipeline** built around **AutoDock Vina**.

The pipeline is designed for:

- High-throughput virtual screening  
- Drug repurposing studies  
- Multi-target docking analysis  
- Manuscript-ready computational workflows  

It integrates **docking, ranking, statistics, interaction profiling, and visualization** into a single, reproducible Bash-driven workflow.

---

## 🚀 Key Features

### Docking
- AutoDock Vina–based docking
- One ligand per CPU core
- Resume-safe execution (skips completed jobs)
- Multi-receptor support

### Ranking & Scoring
- Best-pose energy extraction
- Z-score normalization across receptors
- Stability index (mean / SD)
- Borda rank aggregation
- Final composite ligand score

### Visualization
- Ligand × receptor binding-energy heatmaps
- Rank-matrix heatmaps
- Row-wise and column-wise normalization
- Hierarchical clustering
- PCA, UMAP, t-SNE
- Interactive Plotly dashboards

### Interaction Analysis
- PLIP XML-based interaction profiling
- Residue-wise interaction frequency heatmaps
- Interaction-type fingerprints
- Consensus interaction scores
- Ligand–residue clustering
- Optional LigPlot-style 2D interaction diagrams

### Statistics
- PERMANOVA (scikit-bio)
- Rank robustness analysis
- Cluster validity metrics
- Bootstrapped confidence intervals
- Energy cutoff statistics

### Automation
- Auto-generated figure captions
- Auto-generated **Methods** section
- Auto-generated **Results** narrative
- Headless / HPC-ready execution

---


## 📁 Project Structure

```text
project/
├── vina_pipeline.sh            # Main pipeline script
├── vina_pipeline.yml           # Conda environment file
│
├── ligands/                    # Ligand PDBQT files
├── receptors/                  # Receptor PDBQT files
├── config/                     # AutoDock Vina config files
│   └── <receptor_name>_config.txt
│
├── samplefiles/                # Example test dataset
│   ├── ligands/
│   ├── receptors/
│   └── config/
│
└── results/
    ├── docking/                # Raw docking outputs
    ├── top_N/<receptor>/       # Top-N extracted complexes
    ├── interactions/
    │   └── plip/<receptor>/<ligand>/
    └── summary/
        ├── figures/
        ├── vina_best_scores.csv
        ├── ligand_composite_ranking.csv
        ├── interaction_fingerprints_plip.csv
        ├── methods_auto.txt
        └── results_analysis_auto.txt
```
---

# **Installation**
The pipeline is designed to run inside a dedicated Conda environment for maximum stability.\
\
0\. Download Repository:\
```text
git clone https://github.com/vmevada102/vinapipeline.git
```

1\. Create the Conda environment:\
```text
conda create --name vina_pipeline python=3.10
conda activate vina_pipeline

conda env update --name vina_pipeline --file environment_v1.yml --prune

conda env update --name vina_pipeline --file vina_pipeline.yml --prune
conda install -c conda-forge python-kaleido
conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme

```
\
2\. Install PLIP at the system level:\
```text
sudo apt install -y plip
or
pip install plip
```
3\. Install GNU Parallel and accept citation notice:\
```text
sudo apt install -y parallel
```
4\. Install Vina\
```text
pip install vina
sudo apt install autodock-vina
```

# **Running the Pipeline**
The entire workflow is executed using a single master script:\

```text
bash vina_master_pipeline_auto.sh
```

During execution, the script interactively prompts for input directories and the number of top-ranked ligands to carry forward for advanced analyses such as heatmaps, interaction profiling, and complex generation.

# **Example Test Data**
A minimal, ready-to-run test dataset is provided in the samplefiles/ directory. This allows users to verify correct installation and functionality without preparing custom data. The test dataset completes within minutes on standard hardware.

# **Analyses Performed**
Docking Analysis:\
\- Best binding pose extraction\
\- Per-receptor and consensus energy ranking\
\
Statistical Analysis:\
\- Z-score normalization\
\- Stability index computation\
\- Borda rank aggregation\
\- PERMANOVA and cluster validation\
\
Interaction Analysis:\
\- PLIP residue-level interaction profiling\
\- ProLIF interaction fingerprints\
\- Consensus interaction scoring\
\
Multivariate Analysis:\
\- PCA, UMAP, and t-SNE visualizations


# **Outputs**
The pipeline produces a comprehensive set of outputs, including:\
\
\- CSV files of docking scores and rankings\
\- Heatmaps and clustering figures (300 DPI)\
\- Interaction frequency and fingerprint matrices\
\- LigPlot-style 2D interaction diagrams (optional)\
\- Auto-generated Methods and Results text files

# **Reproducibility and Best Practices**
The pipeline adheres to modern best practices in computational structural biology:\
\- Fully scripted, no manual intervention\
\- Conda-pinned environment\
\- Resume-safe execution\
\- Headless server and HPC compatibility\
\- Transparent statistical methodology

# **Citation**
If you use this pipeline, please cite:\
\
Trott O, Olson AJ. AutoDock Vina. Journal of Computational Chemistry, 2010.\
Schake et al. PLIP 2025. Nucleic Acids Research.\
Tange O. GNU Parallel 2025. Zenodo.

# **Intended Use**
This pipeline is intended for academic and industrial research applications, including:\
\- Virtual screening campaigns\
\- Drug repurposing studies\
\- Binding mode comparison\
\- Manuscript and report generation

# **Analysis using analysis_pipeline.py**
Contributions, bug reports, and feature requests are welcome. Please use GitHub Issues or Pull Requests for communication.


Deterministic statistics

Headless server & HPC compatible

No GUI dependency required

# **Contact and Contributions**
Key Features
1. 🧮 Multi-Criteria Scoring System
Combines three critical components:
Docking Score (binding affinity)
Interaction Score (PLIP-derived interactions)
Stability Score (variance across receptors)

Final Score = 0.5 × Docking + 0.3 × Interaction + 0.2 × Stability

2. 📊 Normalization & Ranking

All scores are normalized (0–1 scale)
Ligands are ranked based on final composite score
Ensures fair comparison across datasets

3. 🧪 Interaction-Aware Scoring
If PLIP interaction data is available:
Assigns weighted scores to interaction types:
Hydrogen bonds
Salt bridges
π-interactions
Hydrophobic contacts

Enhances biological relevance beyond docking scores

4. 🧾 Tier-Based Classification
Ligands are categorized into:
TierDescriptionTier 1 (>0.75)High-confidence bindersTier 2 (0.50–0.75)Moderate candidatesTier 3 (<0.50)Weak binders

5. 📁 Automated Output Generation
The script generates multiple outputs for analysis and reporting:
📊 Main Results
final_integrated_ranking.csv
→ Complete ranked ligand list

top10_summary.txt
→ Quick overview of best candidates

🧾 Comparative Analysis
top_ligand_comparison.csv
→ Combined docking + interaction metrics

🧠 Interpretation Report
final_score_interpretation.txt
→ Human-readable explanation of scoring and results

📂 Input Requirements
The script expects the following directory structure:
results/├── summary/│   └── vina_best_scores.csv└── interactions/    └── plip_interactions_all.csv (optional)
Input Files
FileDescriptionvina_best_scores.csvDocking results (ligand, receptor, energy)plip_interactions_all.csvInteraction data (optional but recommended)

▶️ Usage
Run the script after docking and interaction analysis:
python analysis_pipeline.py

📈 Scientific Significance
This pipeline improves traditional docking workflows by:


Integrating structure-based and interaction-based metrics


Reducing false positives from docking alone


Providing mechanistic insight into ligand binding


Enabling robust lead prioritization



🎯 Use Cases


Virtual screening campaigns


Drug discovery projects


MSc / PhD thesis work


Research publications



🧠 Notes


If PLIP data is not available, the script will still run using docking and stability scores.


Designed to be modular and easily extendable (e.g., ADMET, ML models).



🚀 Future Extensions


ADMET filtering integration


Machine learning-based scoring


Interactive visualization dashboard


Multi-target selectivity analysis



If you want, I can also:


Write a full README.md for your repo


Add badges, installation guide, and screenshots


Or prepare it for GitHub publication (structured project layout)


