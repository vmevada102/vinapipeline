# AutoDock Vina Ligand Preparation Pipeline

A high-throughput ligand preparation pipeline for preparing docking-ready `.pdbqt` files compatible with AutoDock Vina.

This pipeline supports:
- Multi-molecule SDF files
- Multi-SMILES libraries
- MOL / MOL2 / PDB structures
- Batch ligand processing
- Automatic conformer generation
- MMFF94 minimization
- Parallel processing

Designed for:
- Virtual screening
- Molecular docking
- AMR studies
- Phytochemical screening
- Drug discovery workflows

---

# Features

## Supported Input Formats

| Format | Supported |
|---|---|
| SDF | Yes |
| Multi-SDF | Yes |
| MOL | Yes |
| MOL2 | Yes |
| PDB | Yes |
| SMI | Yes |
| SMILES | Yes |
| TXT (SMILES) | Yes |
| CSV (SMILES) | Yes |

---

# Main Capabilities

- Automatic molecule extraction
- Automatic ligand naming
- Multiple conformer generation
- Lowest-energy conformer selection
- MMFF94 optimization
- UFF fallback optimization
- Meeko-based PDBQT generation
- Parallel CPU processing
- Failed ligand tracking
- Molecular descriptor calculation

---

# Output Files

## Generated Ligands

```text
ligands/
    ligand1.pdbqt
    ligand2.pdbqt


# Error Logs

## Detailed Error Log

### `error_log.txt`

Contains:
- traceback information
- molecule parsing errors
- conformer generation failures
- PDBQT conversion failures
- unexpected processing exceptions

---

## Failed Compound List

### `error_convert.txt`

Contains only the names of compounds that failed during processing or conversion.

Useful for:
- reprocessing failed ligands
- quality control
- debugging problematic structures
- filtering invalid compounds

---

# Molecular Descriptors

### `ligand_properties.csv`

Automatically generated molecular descriptor table containing:

- Molecular Weight
- LogP
- TPSA
- H-Bond Donors
- H-Bond Acceptors
- Rotatable Bonds

---

# Installation

## Recommended Environment

```bash
conda create -n vina_pipeline python=3.10 -y
conda activate vina_pipeline

# Install Dependencies

```bash
conda install -c conda-forge rdkit -y
conda install -c conda-forge rdkit gemmi pandas -y

pip install meeko


# Usage

## 1. Place Input Files

Put all ligand files inside the `input/` directory.

Example:

```text
input/
    phytochemicals.sdf
    drug_library.smi
    compounds.csv
## 2. Run the Pipeline

```bash
python ligand_prepare_advanced.py



Supported Input Formats
Format	Supported
SDF	Yes
Multi-SDF	Yes
MOL	Yes
MOL2	Yes
PDB	Yes
SMI	Yes
SMILES	Yes
TXT	Yes
CSV	Yes
Multi-SMILES File Formats

The pipeline supports multiple SMILES input styles.










