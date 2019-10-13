# UMIC-seq
**UMI-linked consensus sequencing tool.**\
Scripts accompanying "Protein engineering with high-accuracy nanopore consensus sequences" (Zurek PJ, Knyphausen P, Neufeld K, Pushpanath A, Hollfelder F, *in preparation*) 

## Setup and installation
This python script is meant to be used as a stand-alone. Just download and copy it into your working directory or add it to your environment path.

### System requirements
You'll need Linux or MacOS for some of the dependencies (scikit-bio). It also runs perfectly well on the Windows Subsystem for Linux though, so that is an option for Windows users.\
An easy way to organise your python packages is with [conda](https://docs.conda.io/en/latest/miniconda.html).

### Python dependencies
Set up a new environment `conda create -n UMIC-seq python=3`, activate it `conda activate UMIC-seq` and add conda-forge to the channel list `conda config --add channels conda-forge`. Install the following packages with `conda install`*`packagename`*:\
*This will automatically install all further dependencies.*\
*Versions numbers are the ones the scripts were tested with, newer versions should work too.*
- python (version 3.7.4)
- biopython (version 1.74)
- scikit-bio (version 0.5.5)
- scikit-allel (version 1.2.1)

Alternatively, the conda environment specifically used in testing this script is provided here. Download UMIC-seq.yml and install with `conda env create -f UMIC-seq.yml`.


## Analysis workflow example

An example dataset in form of 100,000 basecalled reads in fastq format (example_randomreads.fastq.gz) is provided [externally](https://www.dropbox.com/s/d8tkadbvq95p06h/example_randomreads.fastq.gz?dl=1). Additionally, 5000 basecalled reads pre-enriched for 36 clusters are provided within this repository (example_clusterreads.fastq.gz).\
Also, provided are:
- barcodes.fasta: The three barcodes used to demultiplex the sample datasets.
- probe.fasta: The short sequence next to the UMI is provided for extraction on the sample datasets.


### Extraction of UMIs

The first step is to extract the UMI from the reads.
```
python UMIC-seq.py UMIextract --input example_reads.fastq --probe probe.fasta --umi_loc down --umi_len 65 --output ExtractedUMIs.fasta
```
Arguments:
- reads: Provide basecalled reads in fastq format.
- probe: A short sequence (eg 50 bp) adjacent to the UMI.
- umi_loc: Location of UMI in reference to the probe. Upstream (up) or downstream (down).
- umi_len: Length of the UMI to be extracted.
- output: Specify the name of the UMI output fasta file.

Optional:
- min_probe_score: Defaults to length of probe sequence. Minimal alignment score of probe for processing.


### Clustering approximation

Next, you might want to know what a suitable alignment score threshold for clustering the UMIs would be.
```
python UMIC-seq.py clustertest --input ExtractedUMIs.fasta --steps 20 70 10  --output UMIclustertest
```
Arguments:
- input: Fasta file of extracted UMIs.
- steps: Left border, right border and step width for sampled thresholds. Defaults to 20 70 10 (samples thresholds 20 30 40 .. 70).
- output: Prefix for output files.

Optional:
- samplesize: Defaults to 25. Number of clusters to be sampled for threshold approximation.
- threads: Number of threads to use for alignment processing. Defaults to CPU count.


### Full clustering

A full clustering of UMIs can be performed. Similar UMIs will be identified and their corresponding reads will be pooled in a fasta file.
```
python UMIC-seq.py clusterfull --input ExtractedUMIs.fasta --reads example_reads.fastq --aln_thresh 50 --size_thresh 50 --output UMIclusterfull
```
Arguments:
- input: Fasta file of extracted UMIs.
- reads: Fastq file of basecalled reads.
- aln_thresh: Alignment threshold for clustering. UMIs with alignment scores higher than aln_thresh will be clustered.
- size_thresh: Minimal size a cluster can have to be written to file.
- output: Folder name for output files.

Optional:
- threads: Number of threads to use for alignment processing. Defaults to CPU count.


## UMIC-seq_helper

A helper script is provided for some additional functionality.

### Demultiplexing

This is a very simple read demultiplexing script.
```
python UMIC-seq_helper.py demultiplex --barcodes barcodes.fasta --input filteredreads.fastq --threshs 22 29 --output demultiplexedreads
```
Arguments:
- barcodes: Fasta file with barcodes to be used for demultiplexing.
- input: Reads (fastq) to be demultiplexed.
- threshs: Two integer values, lower and upper bound, for alignment score. The script will also output some diagnostics to determine suitable thresholds. Alignment must have a greater score to the target barcode than the upper bound and a lower score to all other barcodes than lower bound.

Optional:
- threads: Number of threads to use for alignment processing. Defaults to CPU count.
- output: Name prefix for demultiplexed reads.


### Generate nanopolish shell script

To automate nanopolish on all generated clusters, the generateSH function can be used.
```
python UMIC-seq_helper.py generateSH --input clustefiles.txt --output run_polish.sh --arguments nanopolish_argumentlist.txt --keyword filename
```
Arguments:
- input: List of cluster file names. Can be generated by `ls *.fasta > clusterfiles.txt`.
- arguments: Text files with arguments to run. Example provided.
- keyword: Keyword in argumentlist to be replaced with entries from the cluster file names.

Optional:
- output: Name of output file.


### Filter mutations

Nanopolish outputs mutation calls in .vcf format. To filter mutations and generate the final sequences in fasta format, vcf2fasta can be used.
```
python UMIC-seq_helper.py vcf2fasta --input vcffiles_list.txt --min_suppfrac 0.6 --reference reference.fasta --output consensus_sequences.fasta 
```
Arguments:
- input: List of vcf file names. Can be generated by `ls -d "$PWD"*.vcf > vcffiles_list.txt`.
- min_suppfrac: Filtering mutations based on support fraction. Only calls with support fraction greater than min_suppfrac will be accepted.
- reference: Reference gene sequence.
- output: Name of output file.

