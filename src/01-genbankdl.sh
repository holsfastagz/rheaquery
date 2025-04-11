#!/bin/bash

accessions="$1"
outdir="$2/homo-seqs"

mkdir $outdir

for id in $(cat $accessions)
do
    [ -f $outdir/$id.fasta ] || \
    (
    echo "Downloading $id . . ."
    esearch -db protein -query $id \
    | efetch -format fasta > $outdir/$id.fasta
    )
done
