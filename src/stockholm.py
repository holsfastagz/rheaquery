#!/usr/bin/env python3

import sys
from pathlib import Path
import glob
from Bio import SeqIO

# record paths
aln = Path(sys.argv[1])                                 # Path to alignment
outdir = Path(sys.argv[2])                              # Path to output directory
outfile = outdir / Path(aln.stem).with_suffix(".sto")   # Output file path

# parse fasta alignment
fasta_records = SeqIO.parse(aln, "fasta")

# rewrite alignment as a stockholm file
SeqIO.write(fasta_records, outfile, "stockholm")
