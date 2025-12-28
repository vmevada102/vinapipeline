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

---

# **Installation**
The pipeline is designed to run inside a dedicated Conda environment for maximum stability.\
\
1\. Create the Conda environment:\
`   `conda env create -f vina\_pipeline.yml\
`   `conda activate vina\_pipeline\
\
2\. Install PLIP at the system level:\
`   `sudo apt install -y plip\
`   `or\
`   `pip install plip\
\
3\. Install GNU Parallel and accept citation notice:\
`   `sudo apt install -y parallel\
`   `parallel --citation
# **Running the Pipeline**
The entire workflow is executed using a single master script:\
\
bash vina\_master\_pipeline\_auto.sh\
\
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
# **Contact and Contributions**
Contributions, bug reports, and feature requests are welcome. Please use GitHub Issues or Pull Requests for communication.


Deterministic statistics

Headless server & HPC compatible

No GUI dependency required
