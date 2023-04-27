# SubsetWT experiments

This repository contains instructions and code to reproduce the experiments in the paper "Subset Wavelet Trees" for SEA 2023.

# Downloading the data

Downloading the data requires `curl`, and `fasterq-dump` from the [SRA toolkit](https://hpc.nih.gov/apps/sratoolkit.html).

## 3682 E. coli genomes

```
curl -O https://zenodo.org/record/6577997/files/coli3682_dataset.tar.gz
```

## 17,336,887 metagenomic reads
```
fasterq-dump ERR5035349
```

# Building and running the experiments

First, pull the submodules with:

```
git submodule init
git submodule update
```

Then, go the SBWT submodule and build it using the instructions in the submodule. Then, compile the experiments with:

```
cd sdsl-lite
cd build
cmake .. -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc)
make
cd ../..
make microbenchmark kmer_search_benchmark
```

This creates an executables called `microbenchmark` and `kmer_search`. Both of these take as input a plain-matrix SBWT index. The SBWT index can be built by passing a list of filenames to the SBWT program, one filename per line. For example, if we have a directory `coli3682_dataset` containing a set of fasta-files, we can create the filename list by running `find coli3682_dataset/ -type f > list.txt`. The SBWT index can then be built with:

```
mkdir temp
./SBWT/build/bin/sbwt build -i list.txt -o index.sbwt -k 31 --add-reverse-complements --n-threads 4 --ram-gigas 8 --temp-dir temp
```

This will save the index to index.sbwt. This should not take more than a few hours for any of the datasets in the paper.

To run the microbenchmark of the paper on this index, run `./microbenchmark index.sbwt`. To run the k-mer search benchmark, run `./kmer_search index.sbwt queries.fna`, where queries.fna is the file containing the queries. In case of the metagenomic read set, we queried a file containing the first 25,000 reads in the dataset, and in case of the E. coli genomes, we queried the genome `coli3682_dataset/GCA_000005845.2_ASM584v2.fna`.



