# Sequence Similarity Network as in Fig. 3
# Author: Paul Zurek (pjz26@cam.ac.uk)
# Date: 17/07/2020

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.manifold import TSNE


#Load in count data containing the mutations of each variant
countdata = pd.read_csv("Zurek_Extended_Data.csv")
#Generate list of mutations from the strings
mutations_lst = [muts.split(" ") for muts in countdata["Mutations"][1:]]
mutations_lst.insert(0,[]) #Handle WT

#### Generate tSNE
REBUILD = False   #Rebuild the tSNE from scratch or use the coordinates from Fig. 3?

if REBUILD:
    #Distance function for tSNE
    def list_distance(l1, l2):
        diffs = set(l1).symmetric_difference(set(l2))
        return len(diffs)

    #Calculate distance matrix
    print("Calculating distance matrix...")
    DM = [[0 for j in range(len(mutations_lst))] for i in range(len(mutations_lst))]
    for i in range(len(mutations_lst)):
        for j in range(len(mutations_lst)):
            DM[i][j] = list_distance(mutations_lst[i], mutations_lst[j])
        print("Row %d of %d" % (i+1, len(mutations_lst)), end="\r")
    DM = pd.DataFrame(DM)
    print("                         ")
    print("Distance matrix generated")

    #tSNE
    tsne = TSNE(n_components=2, metric="precomputed", verbose=1, perplexity=30, learning_rate=200)
    Xtsne = tsne.fit_transform(DM)
    np.save("tSNE_coordinates.npy", Xtsne)
    print("tSNE done and coordinates saved")
else:
    Xtsne = np.load("tSNE_coordinates.npy")



#Which round did the sequence first emerge
conditions = [(countdata["R1"] > 0), ((countdata["R2"] > 0) & (countdata["R1"] == 0)), ((countdata["R3"] > 0) & (countdata["R1"] == 0) & (countdata["R2"] == 0))]
choices = [1, 2, 3]
countdata["whichround"] = np.select(conditions, choices, default=0)
countdata["roundcolor"] = np.select(conditions, ["#2c7bb6","#ffffbf","#d7191c"], default="black")

#Calculate total count
countdata["totalcount"] = countdata["R1"] + countdata["R2"] + countdata["R3"]

#Plot tSNE by round
plt.figure(figsize=(7,6))
plt.scatter(Xtsne[:,0], Xtsne[:,1], s=6, c=countdata["roundcolor"], edgecolors="k", linewidths=0.3) 
plt.scatter(Xtsne[0,0], Xtsne[0,1], s=100, c="#2ca25f", edgecolors="k", linewidths=0.3)
plt.axis('off')
plt.savefig("tSNE_by-round.png", bbox_inches="tight")

### Finding interesting variants
#Find interesting variants by count
print("\n\nHigh count variants")
countdata["Mutations"][0] = "WT"  #Account for WT, otherwise will be dropped by .dropna
highvars = countdata.where(countdata["totalcount"] > 10).dropna()    
print(highvars)

#Find interesting variant by location
target = ((14,0), (26,10))   #(x1, y1) (x2, y2) of rectangle target box
print(f"\n\nLocation variants in box {target}")
for i in range(len(mutations_lst)):
    if target[0][0] < Xtsne[i,0] < target[1][0] and target[0][1] < Xtsne[i,1] < target[1][1]:
        if countdata["totalcount"][i] > 0:
            print("%d: %s" % (countdata["totalcount"][i], countdata["Mutations"][i]))


#### Coloring by founder mutation
founder_colors = ["white" for _ in range(len(mutations_lst))]
for i in range(len(mutations_lst)):         #Order here is important, stops with first condition that is true, so larger sets first
    mut = mutations_lst[i]
    if set(["A64E", "R102S", "D308V"]).issubset(mut):
        #print("%d: %s" % (countdata["totalcount"][i], " ".join(mut)))
        founder_colors[i] = "#ff7f00" #dark yellow
    elif set(["A64E", "R102S", "E323V"]).issubset(mut):
        #print("%d: %s" % (countdata["totalcount"][i], " ".join(mut)))
        founder_colors[i] = "#e31a1c" #dark red
    elif set(["A64E", "R102S"]).issubset(mut):
        #print("%d: %s" % (countdata["totalcount"][i], " ".join(mut)))
        founder_colors[i] = "#fb9a99" #light red
    elif set(["P119Q", "D308V"]).issubset(mut):
        #print("%d: %s" % (countdata["totalcount"][i], " ".join(mut)))
        founder_colors[i] = "#1f78b4" #dark blue
    elif set(["R102S", "E323V"]).issubset(mut):
        #print("%d: %s" % (countdata["totalcount"][i], " ".join(mut)))
        founder_colors[i] = "#fdbf6f" #light yelow
    elif set(["E323V", "P119Q"]).issubset(mut):
        #print("%d: %s" % (countdata["totalcount"][i], " ".join(mut)))
        founder_colors[i] = "#33a02c" #dark green
    elif set(["E323V"]).issubset(mut):
        #print("%d: %s" % (countdata["totalcount"][i], " ".join(mut)))
        founder_colors[i] = "#b2df8a" #light green
    elif set(["P119Q"]).issubset(mut):
        #print("%d: %s" % (countdata["totalcount"][i], " ".join(mut)))
        founder_colors[i] = "#a6cee3" #light blue
    elif set(["D308V"]).issubset(mut):
        #print("%d: %s" % (countdata["totalcount"][i], " ".join(mut)))
        founder_colors[i] = "#cab2d6" #purple


#More candidates: T123I, R102S, A64E

#Colors from https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=9
#a6cee3" #light blue
#1f78b4" #dark blue
#b2df8a" #light green
#33a02c" #dark green
#fb9a99" #light red
#e31a1c" #dark red
#fdbf6f" #light yelow
#ff7f00" #dark yellow
#cab2d6" #purple


sizes = [c+2 for c in countdata["totalcount"]]

plt.figure(figsize=(7,6))
plt.scatter(Xtsne[:,0], Xtsne[:,1], s=sizes, c=founder_colors, edgecolors="k", linewidths=0.3) 
plt.axis('off')
plt.savefig("tSNE_by-founder.png", bbox_inches="tight")
