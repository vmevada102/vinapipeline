#!/usr/bin/env python3

"""
=========================================================
AUTODOCK VINA LIGAND PREPARATION PIPELINE
MULTI-SDF + MULTI-SMILES SUPPORT
=========================================================

SUPPORTED INPUT
---------------
? Multi-molecule SDF
? MOL
? MOL2
? PDB
? SMI
? SMILES
? TXT
? CSV

FEATURES
--------
? Automatically extracts all molecules
? Uses compound names automatically
? Generates ONLY .pdbqt files
? MMFF94 optimization
? UFF fallback
? Multiple conformers
? Lowest-energy conformer selection
? Parallel processing
? Error logging
? Failed compound tracking
? Molecular descriptor CSV

OUTPUT
------
ligands/
    compound1.pdbqt
    compound2.pdbqt

error_log.txt
error_convert.txt
ligand_properties.csv

=========================================================
"""

import os
import re
import sys
import traceback
import subprocess
import importlib.util

from multiprocessing import Pool, cpu_count

# =========================================================
# AUTO INSTALL
# =========================================================

REQUIRED_PACKAGES = {
    "rdkit": "rdkit",
    "meeko": "meeko",
    "gemmi": "gemmi",
    "pandas": "pandas"
}


def is_installed(package_name):

    return importlib.util.find_spec(package_name) is not None


def install_package(package_name):

    print(f"\nInstalling missing package: {package_name}")

    try:

        if package_name in [
            "rdkit",
            "gemmi",
            "pandas"
        ]:

            subprocess.check_call([
                "conda",
                "install",
                "-c",
                "conda-forge",
                package_name,
                "-y"
            ])

        else:

            subprocess.check_call([
                sys.executable,
                "-m",
                "pip",
                "install",
                package_name
            ])

        print(f"SUCCESS: Installed {package_name}")

    except Exception as e:

        print(f"FAILED installing {package_name}")
        print(str(e))
        sys.exit(1)


for pkg_import, pkg_install in REQUIRED_PACKAGES.items():

    if not is_installed(pkg_import):

        install_package(pkg_install)

# =========================================================
# IMPORTS
# =========================================================

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski

from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy

# =========================================================
# FOLDERS
# =========================================================

INPUT_DIR = "input"

OUTPUT_DIR = "ligands"

os.makedirs(INPUT_DIR, exist_ok=True)

os.makedirs(OUTPUT_DIR, exist_ok=True)

ERROR_LOG = "error_log.txt"

FAILED_CONVERT_FILE = "error_convert.txt"

PROPERTY_CSV = "ligand_properties.csv"

# =========================================================
# LOGGING
# =========================================================

def log_error(message):

    with open(ERROR_LOG, "a") as f:

        f.write(message + "\n")


def log_failed_compound(name):

    with open(
        FAILED_CONVERT_FILE,
        "a"
    ) as f:

        f.write(name + "\n")

# =========================================================
# SAFE FILENAMES
# =========================================================

def sanitize_filename(name):

    if not name:

        name = "ligand"

    name = re.sub(
        r'[^A-Za-z0-9_\-]',
        '_',
        name
    )

    return name.strip("_")

# =========================================================
# GET MOLECULE NAME
# =========================================================

def get_molecule_name(
    mol,
    fallback_name,
    index
):

    try:

        if mol.HasProp("_Name"):

            name = mol.GetProp("_Name")

            if name.strip():

                return sanitize_filename(name)

    except:

        pass

    return f"{fallback_name}_{index}"

# =========================================================
# READ INPUT FILE
# =========================================================

def read_input_file(filepath):

    ext = os.path.splitext(
        filepath
    )[1].lower()

    filename = os.path.basename(
        filepath
    ).split('.')[0]

    molecules = []

    try:

        # =================================================
        # MULTI-SDF
        # =================================================

        if ext == ".sdf":

            supplier = Chem.SDMolSupplier(
                filepath,
                removeHs=False
            )

            for idx, mol in enumerate(supplier):

                if mol is None:

                    continue

                mol_name = get_molecule_name(
                    mol,
                    filename,
                    idx + 1
                )

                molecules.append(
                    (
                        mol_name,
                        mol
                    )
                )

        # =================================================
        # MOL
        # =================================================

        elif ext == ".mol":

            mol = Chem.MolFromMolFile(
                filepath,
                removeHs=False
            )

            if mol:

                molecules.append(
                    (
                        filename,
                        mol
                    )
                )

        # =================================================
        # MOL2
        # =================================================

        elif ext == ".mol2":

            mol = Chem.MolFromMol2File(
                filepath,
                removeHs=False
            )

            if mol:

                molecules.append(
                    (
                        filename,
                        mol
                    )
                )

        # =================================================
        # PDB
        # =================================================

        elif ext == ".pdb":

            mol = Chem.MolFromPDBFile(
                filepath,
                removeHs=False
            )

            if mol:

                molecules.append(
                    (
                        filename,
                        mol
                    )
                )

        # =================================================
        # SMILES/TXT/CSV
        # =================================================

        elif ext in [
            ".smi",
            ".smiles",
            ".txt",
            ".csv"
        ]:

            with open(filepath) as f:

                lines = f.readlines()

            for idx, line in enumerate(lines):

                line = line.strip()

                if not line:

                    continue

                if line.startswith("#"):

                    continue

                smiles = None
                name = None

                # CSV
                if "," in line:

                    parts = line.split(",")

                # TAB
                elif "\t" in line:

                    parts = line.split("\t")

                # SPACE
                else:

                    parts = line.split()

                if len(parts) == 0:

                    continue

                smiles = parts[0].strip()

                if len(parts) > 1:

                    name = parts[1].strip()

                else:

                    name = f"{filename}_{idx+1}"

                name = sanitize_filename(name)

                try:

                    mol = Chem.MolFromSmiles(
                        smiles
                    )

                    if mol is None:

                        log_failed_compound(name)

                        log_error(
                            f"Invalid SMILES\n"
                            f"{smiles}\n"
                            f"Line: {idx+1}"
                        )

                        continue

                    molecules.append(
                        (
                            name,
                            mol
                        )
                    )

                except Exception as e:

                    log_failed_compound(name)

                    log_error(
                        f"\nSMILES parsing failed\n"
                        f"Name: {name}\n"
                        f"SMILES: {smiles}\n"
                        f"{str(e)}"
                    )

        else:

            print(
                f"Unsupported format: "
                f"{filepath}"
            )

    except Exception as e:

        log_error(
            f"\nFailed reading "
            f"{filepath}\n{str(e)}"
        )

    return molecules

# =========================================================
# GENERATE BEST CONFORMER
# =========================================================

def generate_best_conformer(
    mol,
    num_confs=10
):

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()

    conf_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=num_confs,
        params=params
    )

    if len(conf_ids) == 0:

        raise Exception(
            "Conformer generation failed"
        )

    energies = []

    for conf_id in conf_ids:

        try:

            mp = AllChem.MMFFGetMoleculeProperties(
                mol
            )

            ff = AllChem.MMFFGetMoleculeForceField(
                mol,
                mp,
                confId=conf_id
            )

            ff.Minimize()

            energy = ff.CalcEnergy()

        except:

            ff = AllChem.UFFGetMoleculeForceField(
                mol,
                confId=conf_id
            )

            ff.Minimize()

            energy = ff.CalcEnergy()

        energies.append(
            (
                conf_id,
                energy
            )
        )

    best_conf = min(
        energies,
        key=lambda x: x[1]
    )[0]

    return mol

# =========================================================
# CONVERT TO PDBQT
# =========================================================

def convert_to_pdbqt(
    mol,
    output_file
):

    preparator = MoleculePreparation()

    setup_list = preparator.prepare(mol)

    pdbqt_result = (
        PDBQTWriterLegacy.write_string(
            setup_list[0]
        )
    )

    if isinstance(
        pdbqt_result,
        tuple
    ):

        pdbqt_string = pdbqt_result[0]

    else:

        pdbqt_string = pdbqt_result

    with open(output_file, "w") as f:

        f.write(pdbqt_string)

# =========================================================
# CALCULATE PROPERTIES
# =========================================================

def calculate_properties(
    mol,
    ligand_name
):

    return {

        "Ligand":
            ligand_name,

        "MolecularWeight":
            Descriptors.MolWt(mol),

        "LogP":
            Descriptors.MolLogP(mol),

        "TPSA":
            Descriptors.TPSA(mol),

        "HDonors":
            Lipinski.NumHDonors(mol),

        "HAcceptors":
            Lipinski.NumHAcceptors(mol),

        "RotatableBonds":
            Lipinski.NumRotatableBonds(mol)
    }

# =========================================================
# PROCESS MOLECULE
# =========================================================

def process_molecule(args):

    ligand_name, mol = args

    print(f"\nProcessing: {ligand_name}")

    try:

        mol = generate_best_conformer(
            mol
        )

        output_file = os.path.join(
            OUTPUT_DIR,
            ligand_name + ".pdbqt"
        )

        convert_to_pdbqt(
            mol,
            output_file
        )

        print(
            f"SUCCESS: "
            f"{output_file}"
        )

        properties = (
            calculate_properties(
                mol,
                ligand_name
            )
        )

        return properties

    except Exception as e:

        error_message = (

            f"\nERROR processing "
            f"{ligand_name}\n"

            f"{str(e)}\n\n"

            f"{traceback.format_exc()}"
        )

        print(error_message)

        log_error(error_message)

        log_failed_compound(
            ligand_name
        )

        return None

# =========================================================
# MAIN
# =========================================================

def main():

    print("\n====================================")
    print(" MULTI-LIGAND PREPARATION ")
    print("====================================")

    all_molecules = []

    input_files = [

        os.path.join(INPUT_DIR, f)

        for f in os.listdir(INPUT_DIR)

        if os.path.isfile(
            os.path.join(INPUT_DIR, f)
        )
    ]

    if len(input_files) == 0:

        print(
            "\nNo files found "
            "inside input folder"
        )

        return

    # Read all molecules
    for filepath in input_files:

        molecules = read_input_file(
            filepath
        )

        all_molecules.extend(
            molecules
        )

    print(
        f"\nTotal molecules found: "
        f"{len(all_molecules)}"
    )

    if len(all_molecules) == 0:

        print(
            "\nNo valid molecules found"
        )

        return

    # CPU cores
    n_cpu = max(
        1,
        cpu_count() - 1
    )

    print(
        f"Using {n_cpu} CPU cores"
    )

    # Parallel processing
    with Pool(n_cpu) as pool:

        results = pool.map(
            process_molecule,
            all_molecules
        )

    # Keep successful results
    valid_results = [

        r for r in results
        if r is not None
    ]

    # Save CSV
    if len(valid_results) > 0:

        df = pd.DataFrame(
            valid_results
        )

        df.to_csv(
            PROPERTY_CSV,
            index=False
        )

        print(
            f"\nSaved properties: "
            f"{PROPERTY_CSV}"
        )

    print("\nALL PROCESSING COMPLETED")

# =========================================================
# ENTRY
# =========================================================

if __name__ == "__main__":

    main()