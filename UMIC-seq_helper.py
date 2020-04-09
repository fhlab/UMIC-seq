#Helper script for UMI-linked consensus sequencing
#Author: Paul Jannis Zurek, pjz26@cam.ac.uk
#21/08/2019 v 1.0
#07/02/2020 v 1.1
## Increased memory efficiency:
## - Using SeqIO.index
## - Using generators (and imap instead of map)
#09/04/2020 v1.1.1
## Improved writing to file speed for demultiplexer

import numpy as np
from Bio import SeqIO
from skbio.alignment import StripedSmithWaterman
import multiprocessing
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser(description="""Helper script for UMI-linked consensus sequencing.
                                 Author: Paul Zurek (pjz26@cam.ac.uk).
                                 Version 1.1.1""")
parser.add_argument('-T', '--threads', type=int, default=0, help='Number of threads to execute in parallel. Defaults to CPU count.')
parser.add_argument('-v', '--version', action='version', version='1.1.1')

subparsers = parser.add_subparsers(help='Select mode', dest='mode')

#Arguments for demultiplexing
demultiplex_parser = subparsers.add_parser('demultiplex', help='Demultiplexing of nanopore reads.')
demultiplex_parser.add_argument('-i', '--input', help='Provide basecalled reads in fastq format.', required=True)
demultiplex_parser.add_argument('-o', '--output', help='Prefix for demultiplexed files.', default='demultiplexed_', required=True)
demultiplex_parser.add_argument('--barcodes', help='Fasta file of barcodes.', required=True)
demultiplex_parser.add_argument('--threshs', help="""Lower and upper alignment score thresholds. 
Will default to the barcode length for lower bound and barcode length * 1.2 for upper bound, which seem reasonable starting points for demultiplexing of nanopore R9.4.1 reads.
Diagnostics to estimate better thresholds will be provided once the script is run""", type=int, nargs=2)

#Arguments for nanopolish script generation
genSH_parser = subparsers.add_parser('generateSH', help='Generate shell script for automated nanopolish execution.')
genSH_parser.add_argument('-i', '--input', help='Text file with list of cluster filenames.', required=True)
genSH_parser.add_argument('-o', '--output', help='Output file name.', default='run_nanopolish.sh', required=True)
genSH_parser.add_argument('--arguments', help='Text file with list of arguments to run.', required=True)
genSH_parser.add_argument('--keyword', help='Keyword to replace in the arguments list with the cluster file names from the input.', required=True)

#Arguments for mutation filter
vcf2fasta_parser = subparsers.add_parser('vcf2fasta', help='Filtering of mutations.')
vcf2fasta_parser.add_argument('-i', '--input', help='Text file with list of vcf files.', required=True)
vcf2fasta_parser.add_argument('-o', '--output', help='Name of the output file.', default='filtered_sequences.fasta', required=True)
vcf2fasta_parser.add_argument('--min_suppfrac', help='Support fraction threshold for filtering.', required=True, type=float, default=0.6)
vcf2fasta_parser.add_argument('--reference', help='Reference sequence fasta.', required=True)

#Parse arguments
args = parser.parse_args()
mode = args.mode
threads = args.threads
if threads == 0:
    threads = multiprocessing.cpu_count()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#DEMULTIPLEX

def align(rec_str):
    aln_scores = []
    for query in ssw_queries:
        aln = query(rec_str)
        aln_scores.append(aln.optimal_alignment_score)
    return aln_scores

def select_bc(scores):
    bc_id = np.argmax(scores)
    if bc_id >= n_barcodes:
        bc_id -= n_barcodes
    sorted_scores = sorted(scores, reverse=True)
    if sorted_scores[0] > highthresh and sorted_scores[1] < lowthresh:
        return bc_id
    else:
        return n_barcodes

if mode == 'demultiplex':
    input_file = args.input     #'./0_preproduction/UMIclustersample.fastq'   #Nanopore bc: 24 bp !
    output_name = args.output
    bc_file = args.barcodes   #'./1_demultiplex/barcodes.fasta'

    #Load barcodes
    print("Generating barcode queries...")
    barcodes = list(SeqIO.parse(bc_file, "fasta"))
    n_barcodes = len(barcodes)
    barcode_names = [bc.id for bc in barcodes] + ['none']
    n_assigns = n_barcodes + 1
    #Generate forward and reverse queries
    ssw_queries = [StripedSmithWaterman(str(bc.seq), score_only=True) for bc in barcodes] + \
                        [StripedSmithWaterman(str(bc.reverse_complement().seq), score_only=True) for bc in barcodes]
    
    #Set default thresholds
    if args.threshs is not None:
        lowthresh, highthresh = args.threshs    #22 (default to ~len(BC))  #29 (default to ~len(BC) * 1.2)    --> Histogram is provided to estimate useful thresholds
    else:
        len_barcode = len(barcodes[0])
        lowthresh = len_barcode
        highthresh = int(len_barcode * 1.2)
    
    #Number of sequences
    print("Loading records...")
    N_seq = len([1 for line in open(input_file)]) // 4
    print(f"{N_seq} entries")

    #Record sequence generator to save memory!
    records = (str(rec.seq) for rec in SeqIO.parse(input_file, 'fastq'))
    
    #Aligning
    print("Aligning...")
    pool = multiprocessing.Pool(processes=threads)
    score_lst = list(tqdm(pool.imap(align, records, chunksize=100), total=N_seq))
    #pool imap takes a generator and returns a generator (thus turning to list here)! saves memory and enables results to be used as they come (not needed here)!
    print("Assigning...")
    assigned_bcs = pool.map(select_bc, score_lst)
    pool.close()
    pool.join()

    #Sample some to be able to estimate good lower and upper thresholds:
    print("\nSampling some scores (to estimate suitable thresholds):")
    for i in range(250):
        print(f"Fwd scores: {score_lst[i][:n_barcodes]} Rev scores: {score_lst[i][n_barcodes:]} Assigned BC: {barcode_names[assigned_bcs[i]]}")

    #Print histogram of BC alignment scores to estimate good thresholds
    print("\nDistribution of all alignment scores:")
    hist, be = np.histogram([x for y in score_lst for x in y], bins=10)
    for i in range(len(hist)):
        print(f"Bin {round(be[i],1)}-{round(be[i+1],1)}:\t {hist[i]}")
    print("\nThresholds for assignment were set to:")
    print(f"Upper threshold: {highthresh} (highest score must be greater)")
    print(f"Lower threshold: {lowthresh} (second highest score must be lower)")

    #Calculate sequence count per barcode
    barcode_dist = []
    for i in range(n_assigns):
        count = assigned_bcs.count(i)
        barcode_dist.append(count)
    print('\nNumber of sequences assigned to barcodes:')
    print(f"Barcodes: {barcode_names}")
    print(f"Sequence count: {barcode_dist}")
    print(f"Relative (%): {[round(nbc/sum(barcode_dist)*100,2) for nbc in barcode_dist]}")

    
    #Generator for writing correct barcode files for memory efficiency
    print("\nWriting to file...")
    handles = []
    for bc_name in barcode_names:
        filename = f'{output_name}_{bc_name}.fastq'
        h = open(filename, "w")
        handles.append(h)

    i = 0
    for rec in SeqIO.parse(input_file, 'fastq'):
        txt = rec.format("fastq")
        bc_id = assigned_bcs[i]
        handles[bc_id].write(txt)
        i += 1

    for h in handles:
        h.close()
    print("Finished.")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CREATE SH

if mode == 'generateSH':
    input_file = args.input #'./4_clustering/BC01_CS-clusters/filenames.txt'
    output_file = args.output #'nanopolish.sh'  #optional
    keyword = args.keyword

    #Arguments
    arguments_file = args.arguments
    arguments = [line.rstrip('\n') for line in open(arguments_file)]

    #Get cluster filenames
    filenames = [line.rstrip('\n') for line in open(input_file)]
    print(f"Generating shell script for {len(filenames)} files, replacing '{keyword}'.")

    #Write arguments with cluster names replaces
    out_handle = open(output_file,'w')
    for fname in filenames:
        for arg in arguments:
            newline = arg.replace(keyword, fname)
            out_handle.write(newline)
            out_handle.write('\n')
            #print(newline)
    out_handle.close()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#VCF2FASTA

def mutate_reference(muts, refstr):
    ref = list(refstr)  #python strings are immutable. For convenience work with the reference as list and at the end convert back to string.
    #Checks
    #1: Positions within ref
    mut_positions = [x[1] for x in muts]
    if len(mut_positions) > 0:
        if max(mut_positions) > len(ref):
            raise Exception(f"Position not within reference: {mut_positions}")
    #2: WT position correct
    seq_wt = [x[0][0] for x in muts]
    ac_wt = [ref[x[1]-1] for x in muts]
    assert seq_wt == ac_wt, "Reference sequences is not the nanopolish reference"

    #Generate variant gene sequence
    variant = ref.copy()
    for mut in reversed(muts):
        pos = mut[1] - 1
        variant[pos] = mut[2]

    return "".join(variant)

  
if mode == 'vcf2fasta':
    from allel import read_vcf
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    #Arguments
    suppfrac = args.min_suppfrac #0.6
    input_file = args.input #'./5_nanopolished/vcffilenames.txt'
    reference_file = args.reference #'refQM.fasta'
    output_file = args.output #'consensus_sequences.fasta'

    #Load
    print("Loading vcf filenames and reference sequence")
    filenames = [line.rstrip('\n') for line in open(input_file)]  #Get this e.g. with ls -d "$PWD"*.vcf
    reference = SeqIO.read(reference_file, "fasta")
    refstr = str(reference.seq)

    records = []
    for fname in filenames:
        #Get mutations from vcf
        if "/" in fname:
            pos = fname.rfind("/")
            clusname = fname[pos+1:-4]
        else:
            clusname = fname[:-4]

        print(f"Loading {fname}       ", end="\r")
        vcf = read_vcf(fname, fields='*')
        muts = []
        if vcf is not None:
            for imut in range(len(vcf["variants/SupportFraction"])):
                ref = vcf["variants/REF"][imut]
                pos = vcf["variants/POS"][imut]
                alt = vcf["variants/ALT"][imut][0]
                if vcf["variants/SupportFraction"][imut] > suppfrac:
                    mut = [ref, pos, alt]
                    muts.append(mut)
        str_muts = ["".join(str(e) for e in v) for v in muts]
        print(f"Filtered mutations: {str_muts}")
        variant_sequence = mutate_reference(muts, refstr)
        #print(variant_sequence)
        record = SeqRecord(Seq(variant_sequence), id=f"consensus_{clusname}", description="")
        records.append(record)
    SeqIO.write(records, output_file, "fasta")
    print(f"Filtered consensus sequences written to {output_file}.")

