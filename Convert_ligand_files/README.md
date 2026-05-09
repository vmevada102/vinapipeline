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





