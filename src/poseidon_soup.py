#!/usr/bin/env python3

import sys
from pathlib import Path
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
from Bio import SeqIO

# record paths
results = Path(sys.argv[1])
outdir = Path(sys.argv[2])
index_html = results / Path("html/full_aln/index.html")
aa_aln = results / Path(f"html/full_aln/data/{results.name}_aln.aa_ali.checked.fasta")
aa_repr = outdir / Path(f"{results.name}_aa_repr.fasta")
out_json = outdir / Path(f"{results.name}_poseidon.json")

# Parse PoSeiDon HTML file
handle = open(index_html)
soup = BeautifulSoup(handle, "html.parser")

# Parse amino acid alignment
faa = SeqIO.parse(aa_aln, "fasta")

# Record representative organism (Longest sequence)
representative_organism = ""
representative_seq = ""
for record in faa:
    if not representative_seq:
        representative_organism = str(record.id)
        representative_seq = str(record.seq)
    elif str(record.seq).count("-") < representative_seq.count("-"):
        representative_organism = str(record.id)
        representative_seq = str(record.seq)

# Write representative sequence to file
with open(aa_repr, "w") as file:
    file.write(f">{representative_organism}\n")
    file.write(representative_seq.replace("-", ""))

# Initialize final data frame
df = pd.DataFrame(columns=["Model", "Significance", "Significant Positions"])

# Iterate over all <tr> elements in HTML file
tr_list = soup.find_all("tr")
for tr in enumerate(tr_list, 1):
    str_tr = str(tr)   # convert to string

    # Check if it is a positive selection element
    if "| LRT" in str_tr:

        # Name of model
        model = str_tr.splitlines()[1].split("<summary>")[1].split(" | ")[0]

        # check significance of model
        if "is significant" in str_tr:
            significance = "significant"
        elif "is mixed!" in str_tr:
            significance = "mixed"
        else:
            significance = "not significant"

        # Find significant positively-selected positions
        raw_position = 0
        sig_positions = []
        for line in str_tr.splitlines():
            raw_position += 1
            if "span" in line:
                final_position = (raw_position - 3
                                  - representative_seq[:raw_position].count("-"))
                sig_positions.append(final_position)

        # Create dictionary to append as new row to data frame
        new_row = {"Model": model,
                   "Significance": significance,
                   "Significant Positions": sig_positions}
        df.loc[len(df)] = new_row

# Print output information
print("\n")
print(f"Description: {results.name}")
print("\n")
print(f"Representative Organism: {representative_organism.replace("E-", "Eurycea ")}")
print("\n")
print(f"Representative Amino Acid Sequence: {representative_seq.replace("-", "")}")
print("\n")
print("Poseidon Results:")
print(df)

# Sava data frame as JSON file
df.to_json(out_json, orient="records")
