import os

import glob

import subprocess

import csv

import sys

import time



# 1. Force Matplotlib to use a non-interactive backend (Prevents display errors)

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt



from concurrent.futures import ThreadPoolExecutor, as_completed



# ================= CONFIGURATION =================

LOCAL_GNINA = "./gnina"

TIMEOUT_SECONDS = 60

VINA_CUTOFF = -7.0   # Threshold for "Good Affinity" (More negative is better)

CNN_CUTOFF = 0.9     # Threshold for "High Confidence" (Higher is better)

# =================================================



def get_gnina_command():

    if os.path.exists(LOCAL_GNINA) and os.access(LOCAL_GNINA, os.X_OK):

        return os.path.abspath(LOCAL_GNINA)

    from shutil import which

    return which('gnina')



def clean_pdbqt(input_path, output_path):

    """Clean PDBQT to fix Multi-Model errors."""

    try:

        with open(input_path, 'r') as infile:

            lines = infile.readlines()

        valid_lines = []

        inside_model = False

        has_models = any(l.startswith("MODEL") for l in lines)

        has_atoms = False

        for line in lines:

            if has_models:

                if line.startswith("MODEL") and "1" in line:

                    inside_model = True; continue 

                if line.startswith("ENDMDL"):

                    if inside_model: break 

                    continue

                if inside_model:

                    if line.startswith(("ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):

                        valid_lines.append(line); has_atoms = True

            else:

                if line.startswith(("ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):

                     valid_lines.append(line); has_atoms = True

        if not has_atoms: return False

        with open(output_path, 'w') as outfile: outfile.writelines(valid_lines)

        return True

    except: return False



def run_gnina_worker(args):

    gnina_cmd, receptor_file, original_ligand, clean_dir = args

    filename = os.path.basename(original_ligand)

    clean_file = os.path.join(clean_dir, f"clean_{filename}")

    

    if not clean_pdbqt(original_ligand, clean_file): return None



    cmd = [gnina_cmd, "-r", receptor_file, "-l", clean_file, "--score_only", "--cpu", "1"]

    

    try:

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_SECONDS)

        data = {"Ligand": filename, "Vina_Affinity": None, "CNN_Score": None, "CNN_Affinity": None}

        

        for line in result.stdout.splitlines():

            line = line.strip()

            if line.startswith("Affinity:"):

                try: data["Vina_Affinity"] = float(line.replace("Affinity:", "").replace("(kcal/mol)", "").strip())

                except: pass

            elif line.startswith("CNNscore:"):

                try: data["CNN_Score"] = float(line.split(":")[1].strip())

                except: pass

            elif line.startswith("CNNaffinity:"):

                try: data["CNN_Affinity"] = float(line.split(":")[1].strip())

                except: pass



        if data["CNN_Score"] is not None:

            return data

        return None

    except: return None



def generate_graphs_pure(csv_file, output_dir, receptor_name):

    """

    Generates Graphs using ONLY standard Python lists (No Numpy dependency).

    """

    print("[-] Generating Graphs (Safe Mode)...")

    try:

        vina_vals = []

        cnn_vals = []

        

        actives_count = 0

        decoys_count = 0

        

        # 1. Read Data into simple lists

        with open(csv_file, 'r') as f:

            reader = csv.DictReader(f)

            for row in reader:

                try:

                    v = float(row['Vina_Affinity'])

                    c = float(row['CNN_Score'])

                    

                    vina_vals.append(v)

                    cnn_vals.append(c)

                    

                    if v < VINA_CUTOFF and c > CNN_CUTOFF:

                        actives_count += 1

                    else:

                        decoys_count += 1

                except:

                    continue



        if not vina_vals:

            print("[!] No valid data found for plotting.")

            return 0, 0



        # 2. Consensus Scatter Plot (Using plot 'o' instead of scatter to avoid array checks)

        plt.figure(figsize=(10, 8))

        

        # Plotting simple lists works even if Numpy is broken

        plt.plot(vina_vals, cnn_vals, 'o', color='#1f77b4', alpha=0.6, markersize=5, label='Compounds')

        

        # Add Lines

        plt.axvline(x=VINA_CUTOFF, color='red', linestyle='--', label=f'Vina Cutoff ({VINA_CUTOFF})')

        plt.axhline(y=CNN_CUTOFF, color='green', linestyle='--', label=f'CNN Cutoff ({CNN_CUTOFF})')

        

        # Formatting

        plt.gca().invert_xaxis()

        plt.title(f'{receptor_name}: Consensus Scoring Analysis')

        plt.xlabel('Vina Affinity (kcal/mol)')

        plt.ylabel('CNN Score')

        plt.legend(loc='upper left')

        plt.grid(True, linestyle='--', alpha=0.5)

        

        plot_path = os.path.join(output_dir, f"{receptor_name}_Consensus_Scatter.png")

        plt.savefig(plot_path, dpi=300)

        plt.close()



        # 3. Classification Bar Chart

        plt.figure(figsize=(8, 6))

        bars = plt.bar(['Likely Actives', 'Decoys/Weak'], [actives_count, decoys_count], color=['#2ca02c', '#d62728'])

        plt.title(f'{receptor_name} Screening Summary')

        plt.ylabel('Count')

        

        for bar in bars:

            yval = bar.get_height()

            plt.text(bar.get_x() + bar.get_width()/2, yval + (yval*0.01), int(yval), ha='center', va='bottom', fontweight='bold')

            

        chart_path = os.path.join(output_dir, f"{receptor_name}_Classification_Chart.png")

        plt.savefig(chart_path, dpi=300)

        plt.close()

        

        print(f"[SUCCESS] Graphs saved to: {output_dir}")

        return actives_count, decoys_count

        

    except Exception as e:

        print(f"[!] Graph generation failed: {e}")

        return 0, 0



def main():

    print("=== FINAL PARALLEL SCREENER STABLE (No-Numpy) ===")

    

    gnina_cmd = get_gnina_command()

    if not gnina_cmd: print("[!] gnina not found"); sys.exit(1)



    receptor_name = input("Enter Receptor Name (e.g., 8a4v): ").strip()

    docking_dir = os.path.join("results", "docking", receptor_name)

    receptor_file = os.path.join("receptors", f"{receptor_name}.pdbqt")

    

    output_dir = os.path.join("results", "ml_results", receptor_name)

    clean_dir = os.path.join(output_dir, "temp_clean_files")

    os.makedirs(clean_dir, exist_ok=True)

    os.makedirs(output_dir, exist_ok=True)

    

    main_csv = os.path.join(output_dir, f"{receptor_name}_All_Results.csv")



    if not os.path.exists(docking_dir): print("[!] Input directory missing"); sys.exit(1)



    # RESUME LOGIC

    processed_ligands = set()

    if os.path.exists(main_csv):

        print("[-] Found existing results. Checking for resume...")

        try:

            with open(main_csv, 'r') as f:

                reader = csv.reader(f)

                header = next(reader, None)

                for row in reader:

                    if row: processed_ligands.add(row[0])

            print(f"[-] Resuming: Skipping {len(processed_ligands)} already processed ligands.")

        except:

            print("[!] CSV corrupt or empty. Starting fresh.")



    all_files = glob.glob(os.path.join(docking_dir, "*.pdbqt"))

    files_to_process = [f for f in all_files if os.path.basename(f) not in processed_ligands]

    

    print(f"[-] Ligands to process: {len(files_to_process)}")

    

    # IF DONE -> GENERATE GRAPHS AND EXIT

    if not files_to_process and len(processed_ligands) > 0:

        print("[!] All ligands already processed. Generating graphs only...")

        actives, decoys = generate_graphs_pure(main_csv, output_dir, receptor_name)

    else:

        # PROCESSING

        max_workers = max(1, os.cpu_count() - 1)

        print(f"[-] Starting threads: {max_workers}")

        

        valid_count = 0

        high_conf_count = 0

        

        with open(main_csv, 'a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames=["Ligand", "Vina_Affinity", "CNN_Score", "CNN_Affinity"])

            if os.stat(main_csv).st_size == 0:

                writer.writeheader()

                

            with ThreadPoolExecutor(max_workers=max_workers) as executor:

                futures = {executor.submit(run_gnina_worker, (gnina_cmd, receptor_file, f, clean_dir)): f for f in files_to_process}

                

                completed = 0

                start_time = time.time()

                

                for future in as_completed(futures):

                    res = future.result()

                    completed += 1

                    

                    if res:

                        valid_count += 1

                        if (res['Vina_Affinity'] is not None and res['CNN_Score'] is not None and 

                            res['Vina_Affinity'] < VINA_CUTOFF and res['CNN_Score'] > CNN_CUTOFF):

                            high_conf_count += 1

                            

                        writer.writerow(res)

                        f.flush()

                    

                    if completed % 20 == 0:

                        elapsed = time.time() - start_time

                        rate = completed / elapsed if elapsed > 0 else 0

                        print(f"    Progress: {completed}/{len(files_to_process)} ({rate:.1f} lig/s) | Valid: {valid_count} | High Conf: {high_conf_count}")

        

        actives, decoys = generate_graphs_pure(main_csv, output_dir, receptor_name)



    # GENERATE HIGH CONFIDENCE CSV

    try:

        high_conf_list = []

        with open(main_csv, 'r') as f:

            reader = csv.DictReader(f)

            for row in reader:

                try:

                    if float(row['Vina_Affinity']) < VINA_CUTOFF and float(row['CNN_Score']) > CNN_CUTOFF:

                        high_conf_list.append(row)

                except: continue

        

        high_conf_list.sort(key=lambda x: float(x['CNN_Score']), reverse=True)

        

        high_conf_file = os.path.join(output_dir, f"{receptor_name}_High_Confidence.csv")

        with open(high_conf_file, 'w', newline='') as f:

            writer = csv.DictWriter(f, fieldnames=["Ligand", "Vina_Affinity", "CNN_Score", "CNN_Affinity"])

            writer.writeheader()

            writer.writerows(high_conf_list)

            

        print(f"[SUCCESS] Top {len(high_conf_list)} Hits saved to {high_conf_file}")

        

        top_ligand = high_conf_list[0]['Ligand'] if high_conf_list else "None"

        top_score = high_conf_list[0]['Vina_Affinity'] if high_conf_list else "N/A"



    except Exception as e:

        print(f"[!] Could not save High Confidence CSV: {e}")

        top_ligand = "None"; top_score = "N/A"



    # MANUSCRIPT REPORT

    txt_file = os.path.join(output_dir, f"{receptor_name}_Manuscript_Report.txt")

    with open(txt_file, 'w') as f:

        f.write(f"Screening Report for {receptor_name}\n")

        f.write("===================================\n")

        f.write(f"Total Processed: {len(processed_ligands) + len(files_to_process)}\n")

        f.write(f"High Confidence Hits: {actives}\n")

        if (actives + decoys) > 0:

            f.write(f"Hit Rate: {round((actives/(actives+decoys))*100, 2)}%\n")

        f.write(f"Top Candidate: {top_ligand} (Vina: {top_score})\n")

        f.write(f"Thresholds used: Vina < {VINA_CUTOFF}, CNN Score > {CNN_CUTOFF}\n")

    

    print(f"[SUCCESS] Manuscript Data -> {txt_file}")

    print(f"[DONE] All outputs saved in: {output_dir}")



if __name__ == "__main__":

    main()
