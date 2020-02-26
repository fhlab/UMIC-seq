#Main script for UMI-linked consensus sequencing
#Author: Paul Jannis Zurek, pjz26@cam.ac.uk
#21/08/2019 v1.0
#07/02/2020 v1.1
## Added quick stop to the clustering, to make it more time efficient. It'll stop clustering when essentially only outliers are left (low average clustersize).

from Bio import SeqIO
from Bio.Seq import Seq
from skbio.alignment import StripedSmithWaterman
import multiprocessing
import itertools
import os
import numpy as np
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="""Main script for UMI-linked consensus sequencing.
                                 Author: Paul Zurek (pjz26@cam.ac.uk).
                                 Version 1.1""")
#parser.add_argument('mode', help='Select mode', choices=('UMIextract', 'clustertest', 'clusterfull'))
parser.add_argument('-T', '--threads', type=int, default=0, help='Number of threads to execute in parallel. Defaults to CPU count.')
parser.add_argument('-v', '--version', action='version', version='1.0')

subparsers = parser.add_subparsers(help='Select mode', dest='mode')

#Arguments for UMIextract
extract_parser = subparsers.add_parser('UMIextract', help='Extraction of UMIs from reads.')
extract_parser.add_argument('-i', '--input', help='Provide basecalled reads in fastq format.', required=True)
extract_parser.add_argument('-o', '--output', help='Specify the name of the UMI output fasta file.', required=True)
extract_parser.add_argument('--probe', help='A short sequence (eg 50 bp) adjacent to the UMI in fasta format.', required=True)
extract_parser.add_argument('--umi_loc', help='Location of UMI in reference to the probe. Upstream (up) or downstream (down).', choices=('up', 'down'), required=True)
extract_parser.add_argument('--umi_len', help='Length of the UMI to be extracted.', type=int, required=True)
extract_parser.add_argument('--min_probe_score', help='Minimal alignment score of probe for processing. Defaults to length of probe sequence.', type=int, default=0)

#Arguments for clustertest
clustest_parser = subparsers.add_parser('clustertest', help='Test suitable alignment thresholds for clustering.')
clustest_parser.add_argument('-i', '--input', help='Fasta file of extracted UMIs.', required=True)
clustest_parser.add_argument('-o', '--output', help='Prefix for output files.', required=True)
clustest_parser.add_argument('--steps', help='Accepts left border, right border and step width for sampled thresholds. Defaults to 20 70 10 (samples thresholds 20 30 40 .. 70).', nargs=3, type=int, required=True)
clustest_parser.add_argument('--samplesize', help='Number of clusters to be sampled for threshold approximation. Defaults to 25.', type=int, default=25)

#Arguments for clusterfull
fullclus_parser = subparsers.add_parser('clusterfull', help='Full clustering of UMIs.')
fullclus_parser.add_argument('-i', '--input', help='Fasta file of extracted UMIs.', required=True)
fullclus_parser.add_argument('-o', '--output', help='Folder name for output files.', required=True)
fullclus_parser.add_argument('--reads', help='Fastq file of basecalled reads.', required=True)
fullclus_parser.add_argument('--aln_thresh', type=int, help='Alignment threshold for clustering. UMIs with alignment scores higher than aln_thresh will be clustered.', required=True)
fullclus_parser.add_argument('--size_thresh', type=int, help='Minimal size a cluster can have to be written to file.', required=True)
fullclus_parser.add_argument('--stop_thresh', type=int, default=5, required=False, help='Defaults to 5. Stops clustering if the average cluster size is smaller than this threshold. Essentially speeds up the clustering by dropping outliers. Set the threshold to 0 if you do not want the program to quit early!')
fullclus_parser.add_argument('--stop_window', type=int, default=20, required=False, help='Defaults to 20. Sets the number of clusters to be used to calculate average cluster size.')


#Parse arguments
args = parser.parse_args()
mode = args.mode
threads = args.threads
if threads == 0:
    threads = multiprocessing.cpu_count()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#UMI EXTRACTION
#Extracts UMI in correct orientation
#To optimize: Multiprocessing? Should be fast enough like this...

def prnt_aln(aln):
    print('Score: %d' % aln.optimal_alignment_score)
    print(aln.aligned_target_sequence)  #READ
    print(aln.aligned_query_sequence)   #Search
    s = aln.target_begin
    d = aln.target_end_optimal
    l = len(aln.target_sequence)
    print('Left Border: %d\t Right Border: %d \n' % (s, l-d))

def extract_left(aln):
    umi_end = aln.target_begin
    umi_begin = umi_end - umi_len
    return umi_begin, umi_end

def extract_right(aln):
    umi_begin = aln.target_end_optimal + 1
    umi_end = umi_begin + umi_len
    return umi_begin, umi_end


if mode == 'UMIextract':
    #Set up variables
    input_file = args.input
    probe_file = args.probe
    output_file = args.output
    umi_loc = args.umi_loc
    umi_len = args.umi_len

    proberec = SeqIO.read(probe_file, "fasta")
    probe_minaln_score = args.min_probe_score
    if probe_minaln_score == 0:
        probe_minaln_score = len(proberec.seq)

    probe_fwd = str(proberec.seq)
    probe_rev = str(proberec.reverse_complement().seq)
    query_fwd = StripedSmithWaterman(probe_fwd)
    query_rev = StripedSmithWaterman(probe_rev)

    umi_list = []
    count, bad_aln, short_aln, success = 0, 0, 0, 0
    for record in SeqIO.parse(input_file, "fastq"):
        rec_seq = str(record.seq)
        alnF = query_fwd(rec_seq)
        alnR = query_rev(rec_seq)
        scoreF = alnF.optimal_alignment_score
        scoreR = alnR.optimal_alignment_score
        #Check basic alignment score
        if scoreF > probe_minaln_score or scoreR > probe_minaln_score:  
            if scoreF > scoreR: #Target in fwd orientation
                #Get umi location
                if umi_loc == 'down':
                    umi_begin, umi_end = extract_right(alnF)   #FWD of downstream is right
                elif umi_loc == 'up':
                    umi_begin, umi_end = extract_left(alnF)   #FWD of upstream is left
                #append to UMI to record list
                if umi_end < len(alnF.target_sequence): #UMI could be out of bounds
                    umi = alnF.target_sequence[umi_begin:umi_end]
                    record.letter_annotations = {} #remove qality
                    record.seq = Seq(umi)
                    umi_list.append(record)
                    success += 1
                else:
                    short_aln += 1
            else: #Target in rev orientation
                #Get umi location
                if umi_loc == 'down':
                    umi_begin, umi_end = extract_left(alnR)   #REV of downstream is left
                elif umi_loc == 'up':
                    umi_begin, umi_end = extract_right(alnR)   #REV of upstream is right
                #append to UMI to record list
                if umi_begin > 0: #UMI could be out of bounds
                    umiR = alnR.target_sequence[umi_begin:umi_end]
                    record.letter_annotations = {} #remove qality
                    record.seq = Seq(umiR).reverse_complement()
                    umi_list.append(record)
                    success += 1
                else:
                    short_aln += 1
        else:
            bad_aln += 1
        count += 1

        if count % 1000 == 0:
            print("%d sequences analysed" % count, end='\r')

    print("%d sequences analysed" % count)
    print('\nUMIs extracted: %d' % success)
    print("Discarded: %.2f%%:" % ((bad_aln+short_aln)/count * 100))
    print("Bad alignment: %d" % bad_aln)
    print("Incomplete UMI: %d" % short_aln)

    SeqIO.write(umi_list, output_file, "fasta")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GENERAL CLLUSTERING FUNCTIONS

#Calculates SW alignment score between sequences i and j
def aln_score(query_nr, umi):
    aln = query_lst[query_nr](umi)
    score = aln.optimal_alignment_score
    return score

def simplesim_cluster(umis, thresh, max_clusters=0, save_scores=False, clussize_thresh=0, clussize_window=20):
    N_umis = len(umis)
    #Need to generate SSW query list, as I cannot pass those objects via pool.map
    global query_lst 
    query_lst = [StripedSmithWaterman(str(umi.seq), score_only=True) for umi in umis]

    clusnr = 0
    lab = [None for i in range(N_umis)]
    seq_index = list(range(N_umis))
    remaining_umis = [str(umi.seq) for umi in umis]
    clussize = []

    score_lst_lst = []
    pool = multiprocessing.Pool(processes=threads)
    #Goes through seq_index and only compares those indices and removes from there!
    while len(seq_index) > 0:
        #Calculate a list of similar umis to the first one left in seq_index 
        score_lst = pool.starmap(aln_score, zip(itertools.repeat(seq_index[0]), remaining_umis))
        if save_scores:
            score_lst_lst.append(score_lst)

        #Only propagate what is similar!
        #Go through the index list with similar sequences
        sim_lst = []    
        for i in reversed(range(len(seq_index))):
            if score_lst[i] > thresh:
                sim_lst.append(seq_index[i])
                lab[seq_index[i]] = clusnr
                del seq_index[i]
                del remaining_umis[i]
        
        #Threshold must be smaller than similarity to itself for the simplesim clusterer to work!! (could be a problem for non-relative metric, just the aln score)
        assert len(sim_lst) >= 1,  'Threshold too high'

        clussize.append(len(sim_lst))
        clusnr += 1
        
        if clussize_thresh > 0:
            #Check last few clusters to see if average is small
            if len(clussize) >= clussize_window:
                av_size = sum(clussize[-clussize_window:]) / clussize_window
                if av_size < clussize_thresh:
                    print("Ended early due to small average cluster size.")
                    break
                print("Cluster %d: %d entries.    \tSeq remaining: %d    \t Average size: %.2f" % (clusnr, len(sim_lst), len(seq_index), av_size), end='\r')
            else:
                print("Cluster %d: %d entries.    \tSeq remaining: %d    " % (clusnr, len(sim_lst), len(seq_index)), end='\r')
        else:
            print("Cluster %d: %d entries.    \tSeq remaining: %d    " % (clusnr, len(sim_lst), len(seq_index)), end='\r')
        
        if (clusnr >= max_clusters) and (max_clusters > 0):  #Sampling mode!
            break

    pool.close()
    if save_scores:
        return clusnr, clussize, lab, score_lst_lst
    return clusnr, clussize, lab



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CLUSTER THRESHOLD APPROXIMATION

#Analyse similarity within one cluster (used in averages_withincluster)
def within_cluster_analysis(clus_id, labels, umis, maxiter=200):
    #Get a list of the actual read ID for reads belonging to the given cluster
    read_ids = list(range(len(labels)))
    cluster_members = [i for i in read_ids if labels[i] == clus_id]
    if len(cluster_members) <= 1:  #if the threshold was so high that no clusters are formed I'll return a 0 as score for now.
        return 0.0, len(cluster_members)
    else:
        #Calculate alignment score for all combinations within this cluster
        #but not all all combinations, maximally for the first 100 sequences in the cluster, otherwise the calculations are too long
        if len(cluster_members) < maxiter:
            maxiter = len(cluster_members)
        count = 0
        total_score = 0
        for i in range(maxiter):
            #Align: Doesn't need to be multicore, as there are only few alignments to be done (nr = maxiter)
            #Main calculation burden is in the clustering, which already is multiprocessed
            query = StripedSmithWaterman(str(umis[cluster_members[i]].seq))
            for j in range(i+1, maxiter):
                aln = query(str(umis[cluster_members[j]].seq))
                score = aln.optimal_alignment_score
                total_score += score
                count += 1
        #print("Score: %.2f" % (total_score/count))
        return total_score/count, len(cluster_members)

#Calculates the within cluster similarity average for the first n clusters
def averages_withincluster(n_samples, labels, umis):
    similarity_lst = []
    length_lst = []
    for i in range(n_samples):
        sim, l = within_cluster_analysis(i, labels, umis)
        similarity_lst.append(sim)
        length_lst.append(l)
    return np.average(similarity_lst), np.average(length_lst), np.std(similarity_lst), np.std(length_lst)

#Plots similarity histogram
def similarity_histogram(score_lst_lst, outname):
    #Remove self-alignment
    for i in range(len(score_lst_lst)):
        maxind = score_lst_lst[i].index(max(score_lst_lst[i]))
        del score_lst_lst[i][maxind]
    #Plotting
    binwidth = 2
    f, axarr = plt.subplots(5, 5, sharex=True, sharey=True, figsize=(12,8))
    f.subplots_adjust(hspace = 0.2, wspace = 0.1)
    f.suptitle("Individual similarity histograms", fontsize=12)
    f.text(0.5, 0.04, 'Alignment score', ha='center')
    f.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
    i = 0
    for a in range(5):
        for b in range(5):
            axarr[a,b].hist(score_lst_lst[i], bins=np.arange(min(score_lst_lst[i]), max(score_lst_lst[i])+binwidth, binwidth))
            axarr[a,b].set_yscale('log')
            i += 1
    plt.savefig(output_name + "_similarityscores_hist.pdf", bbox_inches='tight')


#Function performing a threshold approximation for clustering
def threshold_approx(umis, ssize, left, right, step, outname):
    threshholds = list(range(left, right+step, step))
    similarities_lst, sim_std_lst = [], []
    sizes_lst, size_std_lst = [], []
    print(f"Starting threshold approximation ({threads} threads)")
    for thresh in threshholds:
        _, _, read_labels, scores_lst_lst = simplesim_cluster(umis, thresh, max_clusters=ssize, save_scores=True)
        average_sim, average_size, sim_std, size_std = averages_withincluster(ssize, read_labels, umis)
        similarities_lst.append(average_sim)
        sizes_lst.append(average_size)
        sim_std_lst.append(sim_std)
        size_std_lst.append(size_std)
        print(f"Threshold {thresh}: Average similarity {average_sim:.2f} and size {average_size:.0f}                ")

    if ssize >= 25:
        #Plot a similarity histogram
        similarity_histogram(scores_lst_lst, outname)

    #Plot threshold approximation
    fig, ax1 = plt.subplots(figsize=(6,6))
    ax2 = ax1.twinx()
    ax1.errorbar(threshholds, similarities_lst, yerr=sim_std_lst, color="k", marker="o", linestyle="-", linewidth=2, capsize=3, elinewidth=0.5)
    ax2.errorbar(threshholds, sizes_lst, yerr=size_std_lst, color="b", marker="o", linestyle="-", linewidth=2, capsize=3, elinewidth=0.5)
    ax1.set_xlabel("Alignment score threshold", fontsize=14)
    ax1.set_ylabel("Cluster similarity", color="k", fontsize=14)
    ax2.set_ylabel("Cluster size", color="b", fontsize=14)
    ax2.set_ylim(0, max(sizes_lst)+10)
    plt.savefig(outname + "_thresholdapproximation.pdf", bbox_inches='tight')


if mode == 'clustertest':
    #Setup
    input_file = args.input #'./2_UMIextract/UMIexCS_BC03.fasta'
    sample_size = args.samplesize  #25
    output_name = args.output  #'aproxCS_BC03'
    left, right, step = args.steps   #10, 60, 5

    #Load umis and start approximation
    umis = list(SeqIO.parse(input_file, "fasta"))
    threshold_approx(umis, sample_size, left, right, step, output_name)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FULL CLUSTERING

##Change max_clusters to eg 10 for testing, leave out/set to 0 for production
def cluster_sequences(umis, reads, aln_thresh, size_thresh, max_clusters=0, clussize_thresh=0, clussize_window=20):
    print("Beginning clustering...")
    clus_N, cluster_sizes, labels = simplesim_cluster(umis, aln_thresh, max_clusters=max_clusters, clussize_thresh=clussize_thresh, clussize_window=clussize_window) 

    #Printing some metrics
    print("\nClustering done!")
    print(f"Number of clusters: {clus_N}")
    #print(f"Number of clusters (TEST): {max(labels) + 1}")
    large_clus = [1 for size in cluster_sizes if size >= size_thresh]
    print(f"Clusters with >= {size_thresh} members: {sum(large_clus)}")
    seq_large_clus = sum([size for size in cluster_sizes if size >= size_thresh])
    rel_seq_large_clus = seq_large_clus / len(reads) * 100
    print(f"Total number of sequences in clusters with >= {size_thresh} members: {seq_large_clus} ({rel_seq_large_clus:.2f}%)")

    #Plot clustersizes histogram
    binwidth = 1
    """
    plt.figure()
    plt.hist(cluster_sizes, bins=np.arange(min(cluster_sizes), max(cluster_sizes)+binwidth, binwidth))
    plt.xlabel("Cluster size")
    plt.ylabel("Number of clusters")
    plt.savefig(output_folder + "clustersizes_clusters.pdf", bbox_inches='tight')
    """
    fig = plt.figure()
    hist_clussize = np.repeat(cluster_sizes, cluster_sizes)
    med = np.median(hist_clussize)
    plt.hist(hist_clussize, bins=np.arange(min(hist_clussize), max(hist_clussize)+binwidth, binwidth))
    textstr = f"Median: {med}"
    print(f"Median number of sequences per cluster: {med}")
    plt.text(0.85, 0.85, textstr, transform=fig.transFigure, ha='right')
    plt.xlabel("Cluster size")
    plt.ylabel("Number of sequences")
    plt.savefig(output_folder + "_clustersizes_sequences.pdf", bbox_inches='tight')
    print("Clustersize distributions plotted")

    count = 0
    for i in range(clus_N):
        if cluster_sizes[i] >= size_thresh:
            #Get reads according to UMI label
            count += 1
            clus_member_ids = [rec.id for rec, lab in zip(umis, labels) if lab == i]
            clus_members = [reads[seqid] for seqid in clus_member_ids]
            fname = output_folder + "/cluster_" + str(count) + ".fasta" #FASTA OR FASTQ?
            SeqIO.write(clus_members, fname, "fasta")
    print("Cluster files written")


if mode == 'clusterfull':
    #Setup
    input_UMIfile = args.input   #'./2_UMIextract/UMIexRS_BC03.fasta'
    input_READSfile = args.reads  #'./1_demultiplex/UMIrandomsample_BC03.fastq'
    aln_thresh = args.aln_thresh  #50
    size_thresh = args.size_thresh  #5  
    output_folder = args.output    #'./4_clustering/BC03_RS-clusters/'
    clussize_thresh = args.stop_thresh  #Defaults to 5. Set to 0 if you don't want it to stop clustering!
    clussize_window = args.stop_window  #defaults to 20. 
    
    #Make output directory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder) 

    #USE SeqIO.index FOR ALL! Should be decently memory efficient!
    reads = SeqIO.index(input_READSfile, "fastq")
    umis = list(SeqIO.parse(input_UMIfile, "fasta"))

    #NOW ENDS EARLY IF CLUSSIZE_WINDOW AND THRESH ARE SET!
    #Calculates average clusterisze over the last X clusters and stops clustering if this is below thresh!
    cluster_sequences(umis, reads, aln_thresh, size_thresh, clussize_thresh=clussize_thresh, clussize_window=clussize_window)


