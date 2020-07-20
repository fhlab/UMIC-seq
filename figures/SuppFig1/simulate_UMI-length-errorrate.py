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

threads = 8
umi_lengthses = [20, 30, 40, 50, 60, 70, 80, 90, 100]
off_accuracieses = [0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96]
error_rates = [1-e for e in off_accuracieses]
tripses = [1, 2, 3]


cluster_eff_lst_lst = []
homogeneity_lst_lst = []
completeness_lst_lst = []
for umi_length in umi_lengthses:
    print("\n\n~~~~~~~~~~~~~~~~~~~~~~")
    print("NOW WORKING ON UMI LENGTH %d" % umi_length)
    for trip in tripses:
        cluster_eff_lst = []
        homogeneity_lst = []
        completeness_lst = []
        for off_accuracy in off_accuracieses:
            #Set parameters
            par_N = 1000   #Number of parent sequences
            off_N = 50      #Average number of offspring sequences (mean of normal distribution)
            off_N_scale = 5   #Standard deviation for normal distrubution for offspring number
            name = "%sx%s_ul%d_acc%.2f_%d_" % (par_N, off_N, umi_length, off_accuracy, trip) #For plotting
            print("\nWorking %s" % name)

            #N_par random UMI DNA sequences of length ul
            parents = []
            for i in range(par_N):
                parents.append(sample_umi(umi_length)) 


            #Generate offspring/"reads"
            off_randombase = (1 - off_accuracy) * 4 / 3   #Only 3/4 of random bases are mutations, so adjust rate

            reads = []    #List of sequenced UMIs/offspring
            read_parent_id = []  #List of sequence parents (true clusters)
            bases = "GATC"
            for i in range(par_N):
                rn_off = round(np.random.normal(loc=off_N, scale=off_N_scale))
                for j in range(rn_off):
                    kid = ""
                    for c in parents[i]:
                        if random() < off_randombase:
                            kid += bases[randint(0,3)]
                        else:
                            kid += c
                    reads.append(kid)        #"Sequenced" UMIs: off_N times par_N read with error rate as defined
                    read_parent_id.append(i) #Parent of each of the reads

            print("%d total sequence reads generated" % len(reads))

            #Shuffle these lists! (probably not necessary, but why not)
            reads_N = len(reads)
            shuffleind = list(range(reads_N))
            shuffle(shuffleind)
            zipper = sorted((zip(shuffleind, reads, read_parent_id)))
            z, reads, read_parent_id = (list(t) for t in zip(*zipper))

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #HOMOGENEITY AND COMPLETENESS FOR FULL DATASET WITH PREVIOUSLY CALCULATED GOOD_THRESH!

            #Find suitable threshold
            print("\n\nThreshold approximation:")
            good_thresh = find_threshold(50, 15, 40, 1)
            print("\n\nFull clustering:")
            clus_N, cluster_sizes, read_labels = simplesim_cluster(good_thresh)


            print('\n\n~~~~~~~INFO~~~~~~~')
            hom = metrics.homogeneity_score(read_parent_id, read_labels)
            comp = metrics.completeness_score(read_parent_id, read_labels)
            print("Homogeneity: %0f" % hom)
            print("Completeness: %0f" % comp)
            homogeneity_lst.append(hom)
            completeness_lst.append(comp)

        print("\nDONE TRIPLICATES LOOP\n")

        cluster_eff_lst_lst.append(cluster_eff_lst)
        homogeneity_lst_lst.append(homogeneity_lst)
        completeness_lst_lst.append(completeness_lst)



#PANDAS DATAFRAME FROM THIS?!

n_ul = len(umi_lengthses)
n_trip = len(tripses)
n_acc = len(off_accuracieses)

sub_trip = np.repeat(tripses, n_acc)  #repeat each element of list with np.repeat
d_trip = list(chain.from_iterable(repeat(sub_trip, n_ul)))  #repeat this list with itertools.repeat (and flatten)

d_error = list(chain.from_iterable(repeat(error_rates, n_ul * n_trip)))
d_acc = list(chain.from_iterable(repeat(off_accuracieses, n_ul * n_trip)))
d_ul = list(np.repeat(umi_lengthses, n_acc * n_trip))

d_hom = list(chain.from_iterable(homogeneity_lst_lst))
d_com = list(chain.from_iterable(completeness_lst_lst))


data = {"error_rate": d_error,
"accuracy": d_acc,
"homogeneity": d_hom,
"completeness": d_com,
"triplicate": d_trip,
"umi_length": d_ul
}

df = pd.DataFrame(data=data)
df.to_csv(f"clusterefficiencies_{n_ul}length_{n_trip}trips_{n_acc}accuracies.csv")

