# SubsetWT

This repository contains code to reproduce the experiments in the paper "Subset Wavelet Trees" for SEA 2023.

## Downloading the data

Downloading the data requires `curl`, and `fasterq-dump` from the [SRA toolkit](https://hpc.nih.gov/apps/sratoolkit.html).

### 3682 E. coli genomes

```
curl -O https://zenodo.org/record/6577997/files/coli3682_dataset.tar.gz
```

### 17,336,887 metagenomic reads
```
fasterq-dump ERR5035349
```

### 1,234,695 SARS-CoV-2 genomes. 

These genomes were all SARS-CoV-2 genomes available
at NCBI datasets at the time of writing the article. The collection has grown
significantly since that. The present collection genomes can be downloaded from https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/2697049/.

To obtain the same dataset as the one used in the paper, you need to extract the subset of genomes
that were found at the time of writing the article. The FASTA headers of the sequences
that were in the dataset are listed in `covid_dataset_fasta_headers.txt.gz` in this repository. 
You will need to write a script that extracts the seqences corresponding to those headers from the full database.

# Building

```
git submodule init
git submodule update
cd sdsl-lite
cd build
cmake .. -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc)
make
cd ../..
make
```

This creates an executable called `benchmark`.

# Running

The code takes in a plain-subsetwt SBWT file. There is one small example containing 3 E. coli genomes at `index.sswt`. To benchmark on it, run:

```
./benchmark index.sswt
```