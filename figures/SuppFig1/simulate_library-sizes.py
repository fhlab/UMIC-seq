# Author: Paul Zurek (pjz26@cam.ac.uk)
# Date: 20/07/2020
# # # # # # # # # # # # # # # #

from random import randint, random, shuffle
import numpy as np
from skbio.alignment import StripedSmithWaterman
from multiprocessing import Pool
from itertools import repeat, chain
from sklearn import metrics
import pandas as pd
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib import pyplot as plt
import seaborn as sns

sns.set()
sns.set_context("notebook", font_scale=1.2)
sns.set_style("ticks")#, {"xtick.direction": "in", "ytick.direction": "in", "xtick.top": "True"})


#UMIC-seq functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#QUICK AND GREEDY CD-HIT LIKE CLUSTERING: UMIC-seq

#Calculates SW alignment score between sequences i and j
def aln_score(query_nr, umi_nr): 
    aln = query_lst[query_nr](reads[umi_nr])
    score = aln.optimal_alignment_score
    return score

def simplesim_cluster(thresh, maxiter=100000, verbose=True):
    #Need to generate SSW query list, as I cannot pass those objects via pool.map
    global query_lst 
    query_lst = [StripedSmithWaterman(umi, score_only=True) for umi in reads]
    reads_N = len(reads)
    clusnr = 0
    lab = [-1 for i in range(reads_N)]
    seq_index = list(range(reads_N))
    clussize = []

    if verbose: print('Clustering with %d threads started' % threads)
    pool = Pool(processes=threads)
    #Goes through seq_index and only compares those indices and removes from there!
    while len(seq_index) > 0:
        #Calculate a list of similar umis to the first one left in seq_index      
        score_lst = pool.starmap(aln_score, zip(repeat(seq_index[0]), seq_index))
        
        #Only propagate what is similar!
        #Go through the index list with similar sequences
        sim_lst = []    
        for i in reversed(range(len(seq_index))):
            if score_lst[i] > thresh:
                sim_lst.append(seq_index[i])
                lab[seq_index[i]] = clusnr
                del seq_index[i]
        
        #Threshold must be smaller than similarity to itself for the simplesim clusterer to work!!
        assert len(sim_lst) >= 1,  'Threshold too high'

        clussize.append(len(sim_lst))
        clusnr += 1
        if verbose: print("Cluster %d: %d entries.    \tSeq remaining: %d    " % (clusnr, len(sim_lst), len(seq_index)), end='\r')
        if clusnr >= maxiter:
            if verbose: print("\nStopped clustering at %d" % clusnr)
            break

    pool.close()
    return clusnr, clussize, lab

#Save clusters as tsv file
def save_as_tsv(read_labels, readIDs, filename):
    clusnr = max(read_labels)  #number of clusters is the highest label number for UMIC-seq
    with open(filename, "w") as f :
        for clus in range(clusnr):  #for each cluster
            first = -1
            for i in range(len(read_labels)):  #for each read
                if read_labels[i] == clus:  #if read belongs to current cluster
                    if first == -1: #identify first occurence for tsv format
                        first = i
                    f.write(f"{readIDs[first]}\t{readIDs[i]}\n")   #write cluster in tsv: clusterparent-ID clustervariant-ID


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CALCULATE WITHIN CLUSTER SIMILARITIES OVER A RANGE OF THRESHOLDS FOR X FIRST CLUSTERS: THRESHOLD APPROXIMATION

def within_cluster_analysis(clus_id, labels, maxiter=100):
    #Get a list of the actual read ID for reads belonging to given cluster
    read_ids = list(range(len(labels)))
    cluster_members = [i for i in read_ids if labels[i] == clus_id]
    if len(cluster_members) <= 1:  #if the threshold was so high that no clusters are formed I guess I'll return 0 for similarity then?
        return 0.0, len(cluster_members)
    else:
        #Calculate alignment score for all combinations
        #but not all all combinations, maximally for the first 100 sequences in the cluster, otherwise the calculations are too long
        if len(cluster_members) < maxiter:
            maxiter = len(cluster_members)
        count = 0
        total_score = 0
        for i in range(maxiter):
            query = StripedSmithWaterman(reads[cluster_members[i]], score_only=True)   #SPEED UP HERE WITH MULTIPROCESSING!?
            for j in range(i+1, maxiter):
                aln = query(reads[cluster_members[j]])
                score = aln.optimal_alignment_score
                total_score += score
                count += 1
        #print("Score: %.2f" % (total_score/count))
        return total_score/count, len(cluster_members)

#Calculates the within cluster similarity average for the first n clusters
def averages_withincluster(n_samples, labels):
    similarity_lst = []
    length_lst = []
    for i in range(n_samples):
        sim, l = within_cluster_analysis(i, labels)
        similarity_lst.append(sim)
        length_lst.append(l)
    return np.average(similarity_lst), np.average(length_lst)

#Performs threshold approximation
def find_threshold(ssize = 50, left = 15, right = 40, step = 1):
    threshholds = list(range(left, right, step))
    similarities_lst = []
    sizes_lst = []
    for thresh in threshholds:
        _, _, read_labels = simplesim_cluster(thresh, ssize, False)   #ssize: how many clusters are sampled for approximation
        average_sim, average_size = averages_withincluster(ssize, read_labels)

        print("\nThreshold %d" % thresh)
        print("Similarity: %.2f" % average_sim)
        print("Size: %.2f" % average_size)

        similarities_lst.append(average_sim)
        sizes_lst.append(average_size)

    #Find good threshold
    cutoff = umi_length * 0.02
    good_thresh = threshholds[0]
    for i in range(1,len(similarities_lst)):
        if similarities_lst[i] - similarities_lst[i-1] < cutoff:
            if similarities_lst[i] > (umi_length * 0.5):  #base similarity required
                print("Suitable threshold saturation at thresh: %d (similarity of %.2f)" % (threshholds[i], similarities_lst[i]))
                good_thresh = threshholds[i-1]
                break

    #Plot threshold approximation
    fig, ax1 = plt.subplots(figsize=(6,6))
    ax2 = ax1.twinx()
    ax1.plot(threshholds, similarities_lst, color="k", marker="o", linestyle="-", linewidth=2)
    ax2.plot(threshholds, sizes_lst, color="b", marker="o", linestyle="-", linewidth=2)
    ax1.set_xlabel("Threshold", fontsize=14)
    ax1.set_ylabel("Cluster similarity", color="k", fontsize=14)
    ax2.set_ylabel("Cluster size", color="b", fontsize=14)
    plt.axvline(good_thresh, color="r", linewidth=1)
    plt.savefig("%sn%d_threshapprox.png" % (name, ssize), bbox_inches="tight")

    return good_thresh














#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GENERATION OF RANDOM UMIs AND READS THEREOF

#Generate random DNA sequence
def sample_umi(l):
    umi = ""
    bases = "GATC"
    for _ in range(l):
        umi += bases[randint(0,3)]
    return umi

#PARAMETERS
threads = 8
N_parents = [1000, 10000, 100000, 1000000]
clusterer = "MMSEQS"   #UMIC or MMSEQS
umi_length = 50  
off_accuracy = 0.92  #Average read accuracy (mean of normal distribution)
off_acc_scale = 0.01 #Standard deviation for normal distrubution for read accuracy
tripses = [1, 2, 3]   #triplicate numbers
off_N = 50      #Average number of offspring sequences
off_N_scale = 5   #Standard deviation for offspring number




#Generating random UMIs and reads
print("\n\n~~~~~~~~~~~~~~~~~~~~~~")
print("GENERATING READS")
for par_N in N_parents:
    print(f"Working with {par_N} parents!!")
    for trip in tripses:
        #Set parameters
        name = "%sx%s_ul%d_acc%.2f_%d_" % (par_N, off_N, umi_length, off_accuracy, trip) #For plotting
        print("\nWorking %s" % name)

        #N_par random UMI DNA sequences of length ul
        parents = []
        for i in range(par_N):
            parents.append(sample_umi(umi_length)) 

        #Generate offspring/"reads"
        reads = []    #List of sequenced UMIs/offspring
        read_parent_id = []  #List of sequence parents (true clusters)
        bases = "GATC"
        for i in range(par_N):
            rn_off = round(np.random.normal(loc=off_N, scale=off_N_scale))    #Sample number of reads of this UMI
            for j in range(rn_off):
                acc_off = (np.random.normal(loc=off_accuracy, scale=off_acc_scale))        #Sample error rate for this read
                err_randombase = (1 - acc_off) * 4 / 3   #Only 3/4 of random bases are mutations, so increase rate accordingly
                kid = ""
                for c in parents[i]:
                    if random() < err_randombase:
                        kid += bases[randint(0,3)]
                    else:
                        kid += c
                reads.append(kid)        #"Sequenced" UMIs: off_N times par_N read with error rate as defined
                read_parent_id.append(i) #Parent of each of the reads

        print("%d total sequence reads generated" % len(reads))

        #Shuffle these lists!
        reads_N = len(reads)
        shuffleind = list(range(reads_N))
        shuffle(shuffleind)
        zipper = sorted((zip(shuffleind, reads, read_parent_id)))
        z, reads, read_parent_id = (list(t) for t in zip(*zipper))

        #Write to file
        read_gen = (SeqRecord(Seq(read), id=f"{i}") for i, read in enumerate(reads))
        SeqIO.write(read_gen, f"{name}reads.fasta", "fasta")
        np.save(f"{name}parentIDs.npy", read_parent_id)




#Clustering
#UMIC-seq specific paramters
N_approxiter = 50   #50 clusters formed for threshold approximation
N_fulliter = 500    #500 clusters formed for calculating homogeneity
print("\n\n~~~~~~~~~~~~~~~~~~~~~~")
print("CLUSTERING")
for par_N in N_parents:
    for trip in tripses:
        name = "%sx%s_ul%d_acc%.2f_%d_" % (par_N, off_N, umi_length, off_accuracy, trip) #For plotting
        if clusterer == "MMSEQS":
            #Check mmseqs2 installation
            devnull = open(os.devnull, "w")   #write the which mmseqs command output to nothing
            mafft_path = subprocess.run(["which", "mmseqs"], stdout = devnull).returncode
            assert (mafft_path == 0), "No mmseqs installation found"
            devnull.close()
            #RUN MMSEQS2 WITH THE COMMAND
            ## mmseqs createdb ./bla.fasta blaDB
            ## mmseqs cluster blaDB blaclusDB tmpbla -c 0.8 --cov-mode 0 --min-seq-id 0.8 --seq-id-mode 2 --threads 8
            ## mmseqs createtsv blaDB blaDB blaclusDB clusters-bla.tsv
            subprocess.run(f"mmseqs createdb {name}reads.fasta {name}readDB".split())
            subprocess.run(f"mmseqs linclust {name}readDB {name}clusDB {name}tmp -c 0.8 --cov-mode 0 --min-seq-id 0.8 --seq-id-mode 2 --threads {threads}".split())
            subprocess.run(f"mmseqs createtsv {name}readDB {name}readDB {name}clusDB {name}MMSEQS-clusters.tsv".split())
        elif clusterer == "UMIC":
            #RUN UMIC-seq
            #Read in the UMI sequences
            records = SeqIO.parse(f"{name}reads.fasta", "fasta")
            reads = [str(rec.seq) for rec in records]
            #Find a good clustering threshold
            good_thresh = find_threshold(N_approxiter, 15, 40, 1)
            #Cluster
            clus_N, cluster_sizes, read_labels = simplesim_cluster(good_thresh, maxiter=N_fulliter)
            save_as_tsv(read_labels, [f"{i}" for i in range(len(read_labels))], f"{name}UMIC-clusters.tsv")   #readIDs are created just from the index (see above)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#HOMOGENEITY AND COMPLETENESS

homogeneity_lst_lst = []
completeness_lst_lst = []
for par_N in N_parents:
    homogeneity_lst = []
    completeness_lst = []
    for trip in tripses:
        name = "%sx%s_ul%d_acc%.2f_%d_" % (par_N, off_N, umi_length, off_accuracy, trip)
        print(f"Loading {name}")

        #Load in truth
        read_parent_id = np.load(f"{name}parentIDs.npy", allow_pickle=True)
        N_reads = len(read_parent_id)
        
        ### LOAD IN CLUSTER DATA
        clusters = np.loadtxt(f"{name}{clusterer}-clusters.tsv", delimiter="\t", dtype=str)
        print(f"{N_reads} truths, {len(clusters)} assignments")
        
        #Assign cluster label
        oldclus = ""
        clusterlst = [-1 for _ in range(N_reads)]
        counter = 0
        for i in range(len(clusters)):
            #Change the target sequence to the cluster number
            target = int(clusters[i][1])
            clusterlst[target] = counter
            #If new cluster: advance cluster number
            if clusters[i][0] != oldclus:  
                oldclus = clusters[i][0]
                counter += 1
        
        if clusterer == "UMIC":   #UMIC programmed here to stop early and just sample for accuracy!
            #Get only the real clusters that were assigned
            parent_classes = []
            for l in clusterlst:
                if l != -1:
                    parent_classes.append(read_parent_id[l])
            parent_classes = set(parent_classes)
            print(f"{counter} clusters generated from {len(parent_classes)} parental classes")

            #Calcualtes values if a label is assigned by the clusterer or when it should have been assigned! 
            #Now the completeness should be accurate!
            ad_label, ad_truth = [], []
            for i in range(len(clusterlst)):
                if (clusterlst[i] != -1) or (read_parent_id[i] in parent_classes):    
                    ad_label.append(clusterlst[i])
                    ad_truth.append(read_parent_id[i])
            comp = metrics.completeness_score(ad_truth, ad_label)

            #For homogeneity, the non-assigned "missed" classes should not be pooled. Thus only add classes that have been assigned. Check homogeneity for those!
            ad_label, ad_truth = [], []
            for i in range(len(clusterlst)):
                if clusterlst[i] != -1: 
                    ad_label.append(clusterlst[i])
                    ad_truth.append(read_parent_id[i])
            hom = metrics.homogeneity_score(ad_truth, ad_label)

        elif clusterer == "MMSEQS":
            comp = metrics.completeness_score(read_parent_id, clusterlst)
            hom = metrics.homogeneity_score(read_parent_id, clusterlst)

        print("Homogeneity: %0f" % hom)
        print("Completeness: %0f" % comp)
        homogeneity_lst.append(hom)
        completeness_lst.append(comp)

    print("\nDONE TRIPLICATES LOOP\n")
    homogeneity_lst_lst.append(homogeneity_lst)
    completeness_lst_lst.append(completeness_lst)


#Create pandas dataframe
n_np = len(N_parents)
n_trip = len(tripses)

d_trip = list(chain.from_iterable(repeat(tripses, n_np)))  #repeat this list with itertools.repeat (and flatten)
d_acc = list(repeat(off_accuracy, n_np * n_trip))
d_np = list(np.repeat(N_parents, n_trip))

d_hom = list(chain.from_iterable(homogeneity_lst_lst))
d_com = list(chain.from_iterable(completeness_lst_lst))



data = {"accuracy": d_acc,
"homogeneity": d_hom,
"completeness": d_com,
"triplicate": d_trip,
"n_parents": d_np
}

df = pd.DataFrame(data=data)
df.to_csv(f"librarysize_homogeneity_test-{n_np}-{n_trip}_{clusterer}.csv")




