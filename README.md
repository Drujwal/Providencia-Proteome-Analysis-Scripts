🧬 Providencia Proteome Analysis Scripts
📌 Overview
This repository contains scripts used for proteome-wide analysis of antibiotic resistance in Providencia.
The complete methodological details and results are described in the associated research article.

📂 Contents
scripts/
  core.py               # Core proteome extraction
  antibiotic.py         # Resistance mapping and permutation analysis
  amino.py              # Amino acid composition analysis
  proteingenome.py      # Protein count, pI, and molecular weight
  stat.py               # Statistical tests

reference_core_genes/
  reference.fa

⚙️ Usage
Scripts are designed to be run independently depending on the analysis step:
python core.py
python antibiotic.py
python amino.py
python proteingenome.py
python stat.py

📥 Data Availability
The full dataset used in this study is not included in this repository.
•	A reference core proteome file is provided for demonstration 
•	Additional data will be available upon reasonable request after publication 

📄 Associated Publication
Full methodological details, dataset description, and results are available in:
[Will be updated after Publication]

📬 Contact
For queries or data requests:
Dr. Ujwal Dahal
Email: dr.ujwaldahal@gmail.com

