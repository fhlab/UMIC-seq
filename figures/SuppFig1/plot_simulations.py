# Author: Paul Zurek (pjz26@cam.ac.uk)
# Date: 20/07/2020
# # # # # # # # # # # # # # # #

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

sns.set()
sns.set_context("notebook", font_scale=1.2)
sns.set_style("ticks")#, {"xtick.direction": "in", "ytick.direction": "in", "xtick.top": "True"})


### Supplementary Fig. 1 A
# Load data
df = pd.read_csv(f"clusterefficiency_UMI-length-errorrate.csv", index_col=0)
# Plot
fig, ax = plt.subplots(figsize=(6,6))
sns.lineplot(x="error_rate", y="homogeneity", hue="umi_length", data=df, palette=sns.color_palette("Greys", 9), legend=False, ci=False, linewidth=1)
sns.lineplot(x="error_rate", y="completeness", hue="umi_length", data=df, palette=sns.color_palette("Blues", 9), linewidth=3, ci="sd")
plt.xlabel("Error rate", fontsize=14)
plt.ylabel("Homogeneity / Completeness", fontsize=14)
handles, labels = ax.get_legend_handles_labels()
legend1 = plt.legend(handles=handles[1:], labels=labels[1:], title="UMI length", loc="lower left", fontsize=12, title_fontsize=12)   #Remove first item (misused seabron "title", and set real title to legend)

lines = plt.gca().get_lines()
include = [5, 14]   #grey: 5 or 6 and blue 14 or 15
#legend2 = plt.legend([lines[i] for i in include], [f"label test {str(i)}" for i in include], loc=4)
legend2 = plt.legend([lines[i] for i in include], ["Homogeneity", "Completeness"], loc="lower center", fontsize=12) # title_fontsize=12 title="Metric" 
plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)

plt.ylim([0.64,1.01])
plt.yticks([0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1])
plt.xlim([0.03,0.21])
plt.xticks([0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20])
plt.savefig("SuppFig1A.png", bbox_inches="tight")





### Supplementary Fig. 1 B
# Load data
MMseqs_df = pd.read_csv(f"clusterefficiency_library-sizes_MMseqs2.csv", index_col=0)
UMIC_df = pd.read_csv(f"clusterefficiency_library-sizes_UMIC-seq.csv", index_col=0)

#Group by library size to calculate means and standard deviation
MMseqs_means = MMseqs_df.groupby('n_parents').aggregate("mean")   
UMIC_means = UMIC_df.groupby('n_parents').aggregate("mean")   
MMseqs_std = MMseqs_df.groupby('n_parents').aggregate("std")   
UMIC_std = UMIC_df.groupby('n_parents').aggregate("std")  

n_par = [1000, 10000, 100000, 1000000]

hom_UMIC = UMIC_means["homogeneity"].tolist()
comp_UMIC = UMIC_means["completeness"].tolist()
hom_err_UMIC = UMIC_std["homogeneity"].tolist()
comp_err_UMIC = UMIC_std["completeness"].tolist()

hom_MMseqs = MMseqs_means["homogeneity"].tolist()
comp_MMseqs = MMseqs_means["completeness"].tolist()
hom_err_MMseqs = MMseqs_std["homogeneity"].tolist()
comp_err_MMseqs = MMseqs_std["completeness"].tolist()


plt.figure(figsize=(6,6))
plt.errorbar(n_par, hom_UMIC[:3] + hom_MMseqs[3:], yerr=hom_err_UMIC[:3] + hom_err_MMseqs[3:], c="#636363", linewidth=2, markersize=0, linestyle="dashed", capsize=4, capthick=2, zorder=1, label="Homogeneity")
plt.plot(n_par[:3], hom_UMIC[:3], c="#636363", linewidth=0, markersize=10, marker="o", markerfacecolor="#bdbdbd", markeredgecolor="#636363", markeredgewidth=2, zorder=2)
plt.plot(n_par[3:], hom_MMseqs[3:], c="#636363", linewidth=0, markersize=10, marker="s", markerfacecolor="#bdbdbd", markeredgecolor="#636363", markeredgewidth=2, zorder=3)

plt.errorbar(n_par, comp_UMIC[:3] + comp_MMseqs[3:], yerr=comp_err_UMIC[:3] + comp_err_MMseqs[3:], c="#3182bd", linewidth=2, markersize=0, linestyle="solid", capsize=4, capthick=2, zorder=4, label="Completeness")
plt.plot(n_par[:3], comp_UMIC[:3], c="#3182bd", linewidth=0, markersize=10, marker="o", markerfacecolor="#9ecae1", markeredgecolor="#3182bd", markeredgewidth=2, zorder=5)
plt.plot(n_par[3:], comp_MMseqs[3:], c="#3182bd", linewidth=0, markersize=10, marker="s", markerfacecolor="#9ecae1", markeredgecolor="#3182bd", markeredgewidth=2, zorder=6)

plt.legend(fontsize=12)
plt.xscale("log")
plt.ylim([0.79,1.01])
plt.yticks(np.arange(0.8, 1.01, 0.04))
plt.xlabel("Number of clusters (x 50 reads)", fontsize=14)
plt.ylabel("Homogeneity / Completeness", fontsize=14)
plt.savefig("SuppFig1B.png", bbox_inches="tight")

