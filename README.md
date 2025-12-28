# 🧬 AutoDock Vina Pipeline  
**Parallel Docking · Consensus Ranking · Interaction Fingerprints · Publication-Ready Analytics**

---

## 📌 Overview

This repository provides a **fully automated, high-throughput molecular docking and post-docking analysis pipeline** built around **AutoDock Vina**, designed for:

- Structure-based virtual screening  
- Drug repurposing  
- Multi-target docking studies  
- Manuscript-ready computational analysis  

The pipeline integrates **docking, ranking, statistics, interaction analysis, and visualization** into a single reproducible workflow.

---

## 🚀 Key Features

### Docking
- AutoDock Vina
- One ligand per CPU core
- Resume-safe execution
- Multi-receptor docking support

### Ranking & Scoring
- Best-pose energy extraction
- Z-score normalization across receptors
- Stability index (mean / SD)
- Borda rank aggregation
- Final composite score

### Visualization
- Ligand × receptor binding-energy heatmaps
- Rank-matrix heatmaps
- Row-wise and column-wise normalization
- Hierarchical clustering
- PCA, UMAP, t-SNE
- Interactive Plotly dashboards

### Interaction Analysis
- PLIP XML-based interaction profiling
- Residue-wise interaction heatmaps
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

## 📂 Repository Structure
project/
├── ligands/ # Ligands (.pdbqt)
├── receptors/ # Receptors (.pdbqt)
├── config/ # Vina config files
│ └── receptor_config.txt
├── samplefiles/ # Example test dataset
│ ├── ligands/
│ ├── receptors/
│ └── config/
├── vina_master_pipeline_auto.sh # Main pipeline script
├── vina_pipeline.yml # Conda environment
└── results/
├── docking/
├── top_N/<receptor>/
├── interactions/
│ └── plip/<receptor>/<ligand>/
└── summary/
├── figures/
├── vina_best_scores.csv
├── ligand_composite_ranking.csv
├── interaction_fingerprints_plip.csv
├── methods_auto.txt
└── results_analysis_auto.txt


---
## ⚙️ Installation

### 1️⃣ Create Conda Environment

Create and activate the dedicated Conda environment:

```bash
conda env create -f vina_pipeline.yml
conda activate vina_pipeline

###2️⃣ Install PLIP (System-Level Requirement)

PLIP (Protein–Ligand Interaction Profiler) must be available at the system level.
The pipeline automatically detects PLIP and enables interaction analysis if available.


Ubuntu / Debian (recommended)

Installation

1\. Create Conda Environment:

conda env create -f vina\_pipeline.yml

conda activate vina\_pipeline

2\. Install PLIP (System-Level):

sudo apt install -y plip

or

pip install plip

3\. GNU Parallel:

sudo apt install -y parallel

parallel --citation

Running the Pipeline:

bash vina\_master\_pipeline\_auto.sh

Sample Data:

Use samplefiles/ directory to test pipeline quickly.

Outputs:

\- Docking results

\- Ranking tables

\- Heatmaps and figures

\- PLIP interaction reports

\- Manuscript-ready text files

Citation:

AutoDock Vina, PLIP, GNU Parallel as referenced in README.

Contact:

Open an issue or pull request for collaboration.

