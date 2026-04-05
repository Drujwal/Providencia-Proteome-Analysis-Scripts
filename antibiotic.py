import pandas as pd
from Bio import Phylo
import numpy as np
import os

# 1. Setup
os.chdir(r"D:\New_project\Test")

# 2. Load Tree and Data
tree = Phylo.read("Tree.nwk", "newick")
data = pd.read_csv("Data.csv")

# 3. Clean Data
data['Species'] = data['Species'].astype(str).str.strip()
status_dict = pd.Series(data.Status.values, index=data.Species).to_dict()

# 4. THE CASE-INSENSITIVE FIREWALL
tree_names = [leaf.name for leaf in tree.get_terminals() if leaf.name]

# This line checks for "outgroup" regardless of CAPITALIZATION
all_species = [name for name in tree_names if name in status_dict and name.lower() != "outgroup"]

print("--- VERIFICATION ---")
print(f"Total leaves in tree: {len(tree_names)}")
print(f"In-group species found: {len(all_species)}") 
# IF THIS SAYS 28, IT WORKED!
print("--------------------\n")

# 5. Analysis
resistant_species = [s for s in all_species if status_dict.get(s) == 1]
n_res = len(resistant_species)

def calculate_clumping_score(species_list):
    distances = []
    for i in range(len(species_list)):
        for j in range(i + 1, len(species_list)):
            distances.append(tree.distance(species_list[i], species_list[j]))
    return np.mean(distances) if distances else 0

real_score = calculate_clumping_score(resistant_species)
null_distributions = []

print(f"Running 10,000 permutations for {n_res} resistant strains...")
for _ in range(10000):
    random_sample = np.random.choice(all_species, n_res, replace=False)
    null_distributions.append(calculate_clumping_score(random_sample))

p_value = np.mean(np.array(null_distributions) <= real_score)

print(f"\n--- FINAL RESULTS ---")
print(f"Observed Mean Distance: {real_score:.6f}")
print(f"P-value: {p_value:.4f}")

if p_value < 0.05:
    print("RESULT: SIGNIFICANT CLUSTERING (Vertical Inheritance)")
else:
    print("RESULT: NOT SIGNIFICANT (Likely HGT/Convergent Evolution)")