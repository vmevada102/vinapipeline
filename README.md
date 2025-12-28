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

## 📂 Repository Structure

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



#Installation
##1️⃣ Create the Conda Environment

Create and activate the dedicated environment:
```text
conda env create -f vina_pipeline.yml
conda activate vina_pipeline


This installs all required Python packages for:

Docking analysis

Statistical analysis

Visualization

PLIP XML parsing

ProLIF interaction fingerprints

Multivariate analysis (PCA / UMAP / t-SNE)

2️⃣ Install PLIP (System-Level Requirement)

PLIP (Protein–Ligand Interaction Profiler) must be available at the system level.

Ubuntu / Debian (recommended)
sudo apt update
sudo apt install -y plip

Alternative (pip)
pip install plip


Verify installation:

plip -h


⚠️ If PLIP is not detected, the pipeline automatically skips PLIP analysis safely without crashing.

3️⃣ Install GNU Parallel

GNU Parallel is required for high-performance docking.

sudo apt install -y parallel


Silence the citation notice (recommended):

parallel --citation

▶️ Running the Pipeline

Run the pipeline using:

bash vina_pipeline.sh

🧩 Interactive Prompts

During execution, you will be prompted for:

Ligand directory

Receptor directory

Config directory

Top-N ligands to proceed for:

Heatmap generation

Complex extraction

PLIP / ProLIF interaction analysis

Default values are provided for convenience.

🧪 Example Test Data (Quick Start)

A minimal test dataset is provided under samplefiles/ to verify correct installation.

Sample Data Structure
samplefiles/
├── ligands/
│   ├── sample_ligand_1.pdbqt
│   ├── sample_ligand_2.pdbqt
│   └── sample_ligand_3.pdbqt
├── receptors/
│   └── sample_receptor.pdbqt
└── config/
    └── sample_receptor_config.txt

Run Using Sample Data
bash vina_pipeline.sh


Enter the following when prompted:

Ligand directory [default: ligands/]: samplefiles/ligands
Receptor directory [default: receptors/]: samplefiles/receptors
Config directory   [default: config/]: samplefiles/config
Top-N compounds for downstream analysis [default: 50]: 3


✅ The test run completes in minutes.

📊 Outputs
🔬 Docking Results

Best binding energy per ligand–receptor pair

Resume-safe docking logs

🧮 Ranking & Scoring

Mean binding energy

Z-score normalization

Stability index

Borda rank aggregation

Composite ligand score

All tables are saved under:

results/summary/

📈 Visualization Outputs

Binding-energy heatmaps (Top-N safe mode)

Rank-matrix heatmaps

Residue-interaction heatmaps

Interaction-type heatmaps

PCA / UMAP / t-SNE plots

Cluster validity plots

Figures are saved at 300 DPI in:

results/summary/figures/

🔗 Interaction Analysis

PLIP XML interaction reports

Residue-wise interaction frequencies

Consensus interaction scores

Ligand interaction fingerprints

Optional LigPlot-style 2D interaction diagrams

PLIP outputs are stored in:

results/interactions/plip/<receptor>/<ligand>/

📝 Manuscript-Ready Outputs

The pipeline automatically generates:

methods_auto.txt – ready-to-paste Methods section

results_analysis_auto.txt – structured Results narrative

Per-figure captions saved alongside figures

These files are suitable for direct manuscript integration.

🧾 Citation

If you use this pipeline, please cite:

AutoDock Vina

Trott O, Olson AJ.
AutoDock Vina: improving the speed and accuracy of docking.
Journal of Computational Chemistry, 2010.

PLIP

Schäke et al.
PLIP 2025: protein–ligand interaction profiler.
Nucleic Acids Research, 2025.

GNU Parallel

Tange O.
GNU Parallel 2025.
Zenodo.
https://doi.org/10.5281/zenodo.17692695

🔒 Reproducibility & Best Practices

Conda-pinned environment

Resume-safe execution

Deterministic statistics

Headless server & HPC compatible

No GUI dependency required
