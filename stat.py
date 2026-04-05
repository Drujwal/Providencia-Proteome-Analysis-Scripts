import pandas as pd
from scipy import stats
import os

# 1. Setup
# Ensure this path is correct for your local machine
os.chdir(r"D:\New_project\Test\Again")
df = pd.read_csv("Data.csv")

# Clean data: Remove outgroup and handle Phenotype labels
df = df[df['Species'].str.lower() != 'outgroup']
df['Phenotype'] = df['Phenotype'].astype(str).str.strip()

# Split groups (Resistant vs Sensitive)
res = df[df['Phenotype'] == 'R']
sen = df[df['Phenotype'] == 'S']

# 2. Properties to test (Added Aromatic and Aliphatic)
cols = [
    'Protein', 'MW', 'pI', 
    'Acidic', 'Basic', 
    'Hydrophobic', 'Hydrophilic', 
    'Aromatic', 'Aliphatic',  # Your new columns
    'Cheap', 'Expensive'
]

# 3. Perform Tests and Save Result
output_file = "stats_output.txt"
with open(output_file, "w") as f:
    f.write("STATISTICAL ANALYSIS: RESISTANT (R) vs SENSITIVE (S)\n")
    f.write(f"{'Property':<15} | {'Levene_P':<10} | {'MannWhit_P':<10} | {'Sig'}\n")
    f.write("-" * 65 + "\n")
    
    for c in cols:
        if c in df.columns:
            # Check if we have enough data in both groups
            if len(res[c]) > 0 and len(sen[c]) > 0:
                # Levene's (Equality of Variance - important for your 'Acidic' story)
                _, p_lev = stats.levene(res[c], sen[c])
                
                # Mann-Whitney U (Non-parametric test for medians/distribution)
                _, p_mw = stats.mannwhitneyu(res[c], sen[c], alternative='two-sided')
                
                # Significance stars based on Mann-Whitney P
                sig = "***" if p_mw < 0.001 else "**" if p_mw < 0.01 else "*" if p_mw < 0.05 else "ns"
                
                f.write(f"{c:<15} | {p_lev:<10.4f} | {p_mw:<10.4f} | {sig}\n")
                print(f"Calculated: {c}")
            else:
                print(f"Skipping {c}: Missing data in one phenotype group.")
        else:
            print(f"Skipping {c}: Column not found in CSV.")

print(f"\nSUCCESS: Results saved to {output_file} in your folder.")