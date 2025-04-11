#!/bin/bash

assemblies=$1
outdir=$2/blast
homo_faas=$2/homo-seqs

mkdir $outdir

# Iterate over each gene
for homo_seq in $homo_faas/*
do
    # record the functional annotation from the fasta file of the
    # homologous protein
    annotation=$(head -1 $homo_seq | awk -F'[' '{print $1}' \
    | awk -F' ' '{for (i=2; i<=NF; i++) printf (i==2 ? "" : "_") $i; print ""}')


    [ -f $2/$annotation.fasta ] && rm $2/$annotation.fasta
    touch $2/$annotation.fasta

    for assembly in $assemblies/*.cds
    do
        # isolate species name
        species=$(printf "$assembly" | awk -F'/' '{print $(NF-0)}' \
        | awk -F'_' '{print $1}')

        # generate blastdb name
        blastdb_name="$assemblies/blastdbs/$species"

        # make blast database if it doesn't exist
        [ -f $assemblies/blastdbs/$species.nhr ] || \
        makeblastdb -in $assembly -dbtype nucl -out $blastdb_name

        # Generate blast output file name
        blast_out="$outdir/${species}_$annotation.outfmt6"

        # blast the homologous protein sequence against our
        # transcriptome if it has not already been done
        [ -f "$blast_out" ] || \
        (
        tblastn -query $homo_seq -db $blastdb_name -out $blast_out \
        -outfmt 6 -evalue 1e-3 &&
        echo -e "\nQuerying $annotation against $species assembly."
        )

        ## append species sequence to gene fasta file
        ##top_id=$(head -1 $blast_out | awk -F'\t' '{print $2}')
        top_id=$(awk '{sum[$2] += $4} END {for (s in sum) print s, sum[s]}' "$blast_out" | sort -k2,2nr | head -n 1 | awk '{print $1}')
        [ -z "$top_id" ] && echo "No sequence found." && continue
        out_seq=$(sed -n "/$top_id/,/^>/p" $assembly | sed "1d;\$d" | \
            tr -d '\n')
        echo "$out_seq"
        echo "Sequence length: ${#out_seq} bp" 
        printf ">$species\n" >> $2/$annotation.fasta
        printf "$out_seq\n" >> $2/$annotation.fasta
    done
done
