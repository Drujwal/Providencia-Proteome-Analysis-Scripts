import os, subprocess, shutil
from Bio import SeqIO
import multiprocessing as mp
from functools import partial

# --- CONFIGURATION ---
# UPDATED: No more space in the folder name
BASE_DIR = r"D:\New_project"
GENOME_DIR = os.path.join(BASE_DIR, "Analysis")
RESULTS_DIR = os.path.join(BASE_DIR, "Core_Genes_Results")
DB_DIR = os.path.join(BASE_DIR, "blast_databases")

REF_FILE = os.path.join(GENOME_DIR, "reference.fa")
FINAL_MATRIX = os.path.join(BASE_DIR, "supermatrix.fasta")
BLAST_PATH = r"C:\Program Files\NCBI\blast-2.17.0+\bin"

THREADS = 4
OCCUPANCY_THRESHOLD = 23  
IDENTITY_THRESHOLD = 70.0 
LENGTH_THRESHOLD = 0.70   

def get_exe(name):
    return os.path.join(BLAST_PATH, f"{name}.exe")

def process_gene(ref_rec, all_genomes):
    safe_id = "".join([x for x in ref_rec.id if x.isalnum()])[:20]
    tmp_q = os.path.join(DB_DIR, f"q_{os.getpid()}_{safe_id}.fa")
    species_found = {}
    
    try:
        with open(tmp_q, "w") as f:
            f.write(f">query\n{str(ref_rec.seq)}")
        
        for g_file in all_genomes:
            db_path = os.path.join(DB_DIR, f"{os.path.splitext(g_file)[0]}_db")
            
            # Simplified command (no more space-fighting needed)
            cmd = [
                get_exe("blastn"),
                "-query", tmp_q,
                "-db", db_path,
                "-outfmt", "6 pident length qlen sseq",
                "-max_target_seqs", "1"
            ]
            
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True)
                out = proc.stdout.strip()
                
                if out:
                    parts = out.split('\t')
                    pident = float(parts[0])
                    aln_len = int(parts[1])
                    q_len = int(parts[2])
                    sseq = parts[3].replace("-", "")
                    
                    if pident >= IDENTITY_THRESHOLD and (aln_len / q_len) >= LENGTH_THRESHOLD:
                        species_found[os.path.splitext(g_file)[0]] = sseq
            except:
                continue
                
        if len(species_found) >= OCCUPANCY_THRESHOLD:
            out_path = os.path.join(RESULTS_DIR, f"Gene_{safe_id}.fasta")
            with open(out_path, "w") as f:
                f.write(f">reference\n{str(ref_rec.seq)}\n")
                for sp, seq in species_found.items():
                    f.write(f">{sp}\n{seq}\n")
            return True
    finally:
        if os.path.exists(tmp_q):
            try: os.remove(tmp_q)
            except: pass
    return False

def concatenate_genes():
    print("\nStarting concatenation into supermatrix...")
    species_data = {}
    fasta_files = [f for f in os.listdir(RESULTS_DIR) if f.endswith(".fasta")]
    
    if not fasta_files:
        print("No core genes found to concatenate.")
        return

    for f in fasta_files:
        for rec in SeqIO.parse(os.path.join(RESULTS_DIR, f), "fasta"):
            if rec.id not in species_data:
                species_data[rec.id] = ""
            species_data[rec.id] += str(rec.seq)
            
    with open(FINAL_MATRIX, "w") as out:
        for sp, seq in species_data.items():
            out.write(f">{sp}\n{seq}\n")
    print(f"DONE! Supermatrix saved to: {FINAL_MATRIX}")

if __name__ == "__main__":
    for d in [RESULTS_DIR, DB_DIR]:
        if not os.path.exists(d): 
            os.makedirs(d)

    if not os.path.exists(REF_FILE):
        print(f"ERROR: Cannot find {REF_FILE}")
        exit()

    ref_genes = list(SeqIO.parse(REF_FILE, "fasta"))
    genomes = [f for f in os.listdir(GENOME_DIR) if f.lower().endswith(".fa") and f != "reference.fa"]
    
    print(f"Found {len(genomes)} genomes. Building databases...")
    for g in genomes:
        db_path = os.path.join(DB_DIR, f"{os.path.splitext(g)[0]}_db")
        if not os.path.exists(db_path + ".nsq"):
            print(f" -> Building: {g}")
            make_cmd = [
                get_exe("makeblastdb"),
                "-in", os.path.join(GENOME_DIR, g),
                "-dbtype", "nucl",
                "-out", db_path
            ]
            subprocess.run(make_cmd, capture_output=True)

    print(f"Starting analysis on {len(ref_genes)} genes using {THREADS} threads...")
    mp.set_start_method('spawn', force=True)
    
    with mp.Pool(processes=THREADS) as pool:
        results = pool.map(partial(process_gene, all_genomes=genomes), ref_genes)
    
    total_found = sum(results)
    print(f"\nAnalysis Complete. Found {total_found} core genes.")
    
    if total_found > 0:
        concatenate_genes()