# UMIC-seq
**UMI-linked consensus sequencing tool.**\
Scripts accompanying "Protein engineering based on high fidelity nanopore sequencing using UMI-linked consensus reads" (Zurek PJ, Knyphausen P, Neufeld K, Pushpanath A, Hollfelder F, *in preparation*) 

## Setup and installation
This python script is meant to be used as a stand-alone. Just download and copy it into your working directory or add it to your environment path.

### System requirements
You'll need Linux or MacOS for some of the dependencies (scikit-bio). It also runs perfectly well on the Windows Subsystem for Linux though, so that is an option for Windows users.\
An easy way to organise your python packages is with [conda](https://docs.conda.io/en/latest/miniconda.html).

### Python dependencies
Set up a new environment `conda create -n UMIC-seq python=3` and install the following packages with `conda install`*`packagename`*:\
*Versions numbers are the ones the scripts were tested with, newer versions should work too.*
- python (version 3.6)
- biopython (version 1.73)
- scikit-bio (version 0.5.2)
- scikit-allel
- Matplotlib, numpy, scipy

Alternatively, the conda environment specifically used in testing this script is provided here. Download UMIC-seq.yml and install with `conda env create -f UMIC-seq.yml`.

Activate the new environment `conda activate UMIC-seq`.

## Analysis workflow example

An example dataset in form of 100,000 basecalled reads in fastq format (example_reads.fastq) are provided within this repository. Additionally a short sequence located next to the UMI is provided (probe.fasta).


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
- steps: Accepts left border, right border and step width for sampled thresholds. Defaults to 20 70 10 (samples thresholds 20 30 40 .. 70).
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


## iPython notebook

A notebook with the example in more detail will be made available.
