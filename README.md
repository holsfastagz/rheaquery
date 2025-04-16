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

## Dependencies

- numpy
- pandas
- biopython
- beautifulsoup4

All dependencies can be downloaded using the following command:

```sh
pip install numpy pandas biopython beautifulsoup4
```

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

- `-a`: a text file containing the NCBI accession numbers of a homologous
protein.

- `-f`: a directory containing the transcriptome assemblies.

- `-o`: the output directory (default .).

2. pfam2go - annotate protein domains with GO terms using the pfam2go database,
with the results of hmmscan. 

```sh
rheaquery -m pfam2go -p gene_hmm.out -d pfam2go -o outdir/
```

- `-p`: the hmmscan output file.

- `-d`: the pfam2go database.

- `-o`: output directory (default .).

3. soup - parse PoSeiDon results and record amino acid positions of positively-
selected codons in a representative sequence.

```sh
rheaquery -m soup -f gene_out/ -o outdir/
```

- `-f`: PoSeiDon output directory for a gene, e.g. (results/pax6).

- `-o`: output directory (default .).
