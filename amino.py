import os
import pandas as pd
from collections import Counter
from Bio import SeqIO

# --- CONFIGURATION ---
root_path = r'D:\New_project\Aminoacids'
output_file = "Providencia_AA_Final_Results.csv"

# These will be your Column Headers in Excel
amino_acids = ['Alanine_A', 'Cysteine_C', 'Aspartic_D', 'Glutamic_E', 'Phenylalanine_F', 
               'Glycine_G', 'Histidine_H', 'Isoleucine_I', 'Lysine_K', 'Leucine_L', 
               'Methionine_M', 'Asparagine_N', 'Proline_P', 'Glutamine_Q', 'Arginine_R', 
               'Serine_S', 'Threonine_T', 'Valine_V', 'Tryptophan_W', 'Tyrosine_Y']

# Mapping short codes to the full labels above
aa_map = {'A':'Alanine_A', 'C':'Cysteine_C', 'D':'Aspartic_D', 'E':'Glutamic_E', 
          'F':'Phenylalanine_F', 'G':'Glycine_G', 'H':'Histidine_H', 'I':'Isoleucine_I', 
          'K':'Lysine_K', 'L':'Leucine_L', 'M':'Methionine_M', 'N':'Asparagine_N', 
          'P':'Proline_P', 'Q':'Glutamine_Q', 'R':'Arginine_R', 'S':'Serine_S', 
          'T':'Threonine_T', 'V':'Valine_V', 'W':'Tryptophan_W', 'Y':'Tyrosine_Y'}

def calculate_frequencies():
    all_species_data = []
    
    # Get all .fa files
    files = [f for f in os.listdir(root_path) if f.lower().endswith(('.fa', '.fasta'))]
    
    print(f"Analyzing {len(files)} files...")

    for filename in files:
        file_path = os.path.join(root_path, filename)
        genome_counts = Counter()
        total_aa_count = 0
        species_name = os.path.splitext(filename)[0]
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                # Translate DNA to Protein (Table 11)
                prot_seq = record.seq.translate(table=11, to_stop=True)
                prot_str = str(prot_seq).upper()
                genome_counts.update(prot_str)
                total_aa_count += len(prot_str)
            
            if total_aa_count > 0:
                row = {"Species_Name": species_name}
                # Fill each specific Amino Acid column
                for code, label in aa_map.items():
                    count = genome_counts.get(code, 0)
                    row[label] = round((count / total_aa_count) * 100, 4)
                
                all_species_data.append(row)
                print(f"Processed: {species_name}")

        except Exception as e:
            print(f"Error in {filename}: {e}")

    if all_species_data:
        df = pd.DataFrame(all_species_data)
        # Force the column order: Species Name first, then A-Y
        column_order = ["Species_Name"] + amino_acids
        df = df[column_order]
        
        # Save with headers explicitly turned ON
        df.to_csv(output_file, index=False, header=True)
        print(f"\nSUCCESS: Open '{output_file}' in Excel to see the labels!")
    else:
        print("No files were processed. Check your file extensions!")

if __name__ == "__main__":
    calculate_frequencies()