# rheaquery
A small pipeline for positive selection analysis in Biodiversity Discovery.

## Installation

Clone this repository onto your local machine or download and unzip the latest
release.

Run the installation script.

```sh
bash install.sh
```

Reinitialize your shell.

```sh
source ~/.bashrc
```

Rheaquery should now be installed. You can use `rheaquery -h` to make sure.

## Usage

Rheaquery's basic usage is:

```sh
rheaquery -m mode <options>
```

Rheaquery makes use of three software modules.

1. blast - identify homologous sequences in our assemblies based on a reference
from NCBI.

```sh
rheaquery -m blast -a accession.txt -f assemblies/ -o outdir/
```

2. pfam2go - annotate protein domains with GO terms using the pfam2go database,
with the results of hmmscan. 

```sh
rheaquery -m pfam2go -p gene_hmm.out -d pfam2go -o outdir/
```

3. soup - parse PoSeiDon results and record amino acid positions of positively-
selected codons in a representative sequence.

```sh
rheaquery -m soup -f gene_out/ -o outdir/
```
