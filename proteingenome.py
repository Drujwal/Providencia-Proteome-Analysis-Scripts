import os
import csv
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# 1. PATHS - Everything is in 'genomes' now
folder_path = r"D:\New_project\Isoelectric\genomes"
output_file = r"D:\New_project\Isoelectric\Whole_Genome_Results.csv"

all_results = []

print(f"--- Starting Analysis on 4MB Genome Files ---")

# 2. SCAN THE FOLDER
for filename in os.listdir(folder_path):
    if filename.lower().endswith(".fa"):
        filepath = os.path.join(folder_path, filename)
        
        p_count = 0
        mw_sum = 0
        pi_sum = 0
        
        # Open the 4000KB file and read every single gene
        with open(filepath, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                try:
                    # DNA -> Protein (Stops at '*' to prevent pI tool crash)
                    prot = record.seq.translate(to_stop=True)
                    
                    # Clean: Keep only 20 standard Amino Acids
                    clean_p = "".join([c for c in str(prot) if c in "ACDEFGHIKLMNPQRSTVWY"])
                    
                    if len(clean_p) > 20:
                        analysis = ProteinAnalysis(clean_p)
                        mw_sum += analysis.molecular_weight()
                        pi_sum += analysis.isoelectric_point()
                        p_count += 1
                except:
                    continue
        
        # 3. COLLECT DATA
        if p_count > 0:
            all_results.append({
                "Species": filename,
                "Dataset": "Whole_Genome",
                "Total_Proteins": p_count,
                "Avg_MW_kDa": round((mw_sum / p_count) / 1000, 2),
                "Avg_pI": round(pi_sum / p_count, 2)
            })
            print(f"   Done: {filename} -> {p_count} genes found.")

# 4. SAVE TO CSV
keys = ["Species", "Dataset", "Total_Proteins", "Avg_MW_kDa", "Avg_pI"]
with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=keys)
    writer.writeheader()
    writer.writerows(all_results)

print(f"\nSUCCESS! Open {output_file} to see the full genome results.")