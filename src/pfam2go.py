#!/usr/bin/env python3

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SearchIO

# record paths
hmmfile = Path(sys.argv[1])                                  # Path to HMMR output file
outdir = Path(sys.argv[2])                                   # Path to output directory
pfam2go = Path(sys.argv[3])                                  # Path to pfam2go file
outfile = outdir / Path(hmmfile.stem).with_suffix(".json")   # Path to output file

# Parse hmmscan output file
hmm_contents = SearchIO.parse(hmmfile, "hmmer3-text") 

# Initialize data frame
hmm_results = pd.DataFrame(columns=["Pfam Accession", "Hit ID",
                                    "Hit Description", "Ali from", "Ali to",
                                    "GO IDs", "GO Terms"])

# Iterate over each query result from hmmscan
for qresult in hmm_contents:

    # Iterate over each query hit
    for hit in qresult.hits:
        hit_id = hit.id                     # Domain ID
        hit_description = hit.description   # Domain Description

        # Iterate over each hit HSP (idrk what this means)
        for hsp in hit.hsps:
            ali_from = hsp.query_start + 1   # First residue of domain annotation
            ali_to = hsp.query_end           # Last residue of domain annotation

        # Parse pfam2go file
        go_ids = []     # Initialize list of GO IDs
        go_terms = []   # Initialize list of GO terms
        with open(pfam2go, "r") as f:
            for line in f.readlines():
                if (" " + hit_id + " ") in line:
                    pfam_accession = line.split(" ")[0].replace("Pfam:", "")   # Pfam accession number
                    go_ids.append(line.split(" > ")[1].split(" ; ")[1]         # Append GO ID to GO ID list
                                    .replace("\n", ""))
                    go_terms.append(line.split(" > ")[1].split(" ; ")[0]       # Append GO term to GO terms list
                                    .replace("GO:", ""))

        # Record data in dictionary to add to data frame
        new_row = {"Pfam Accession": pfam_accession,     # Pfam accession number
                   "Hit ID": hit_id,                     # Domain ID
                   "Hit Description": hit_description,   # Domain description
                   "Ali from": ali_from,                 # First residue of domain annotation
                   "Ali to": ali_to,                     # Last residue of domain annotation
                   "GO IDs": go_ids,                     # List of GO IDs
                   "GO Terms": go_terms}                 # List of GO terms

        # Append new row to data frame
        hmm_results.loc[len(hmm_results)] = new_row

# write data frame to JSON format 
hmm_results.to_json(outfile, orient="records")
