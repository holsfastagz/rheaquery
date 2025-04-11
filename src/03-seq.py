#!/usr/bin/env python3

import sys
from pathlib import Path
import glob
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo

# record file paths
outdir = Path(sys.argv[1])                    # output directory
blast_out = Path(f"{sys.argv[1]}/blast")   # blast output directory
assemblies = Path(sys.argv[2])             # assemblies directory

# iterate over blast outfmt files
for outfmt_file in blast_out.glob("*.outfmt6"):

    # isolate species and gene names
    species = str(outfmt_file.name).split("_")[0]
    gene = str(outfmt_file.name).split("_")[1]

    print(f"Gathering sequences from {species}. . .")

    # read outfmt file and read to memory
    id_list = []
    with open(outfmt_file, "r") as f:
        outfmt_contents = f.readlines()

    # iterate over entry in outfmt and record id to list
    for entry in outfmt_contents:
        trinity_id = entry.split("\t")[1]
        id_list.append(trinity_id)

    # create output MSA file
    multi_fasta = outdir / "msa" / f"{species}_{gene}.multi.fasta"
    multi_fasta.parent.mkdir(parents=True, exist_ok=True)

    record_count = 0

    # iterate over each ID in outfmt file
    for id in id_list:
        # iterate over each sequence in assembly
        for record in SeqIO.parse(next(assemblies.glob(f"{species}*")), "fasta"):
            if id in record.description:
                record_count += 1
                
                # write header and sequence to output msa
                with multi_fasta.open("a") as f:
                    f.write(f">{record.description}\n")
                    f.write(f"{record.seq}\n")
    print(f"{record_count} sequences found (e-value <= 1e-10)")
