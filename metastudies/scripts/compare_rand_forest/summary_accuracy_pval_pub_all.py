#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com

#This is a script for comparing different transformations
# --------------------------------------------------------------------------
print(f"""Running {__file__}.
This is a script for comparing random forest output with pvalues. 
Currently only works for the python output. It should summerize the data from
transformation_pvalue_plot.py. The dataframe that feeds into the boxplot
should be a dataframe with a column for each comparison/transformation dataset.
In the columns should be positive or negative pvalues. The pvalues will be 
positive if the difference between the column-tranformation are positive and
negative if the difference is negative.
""")
# --------------------------------------------------------------------------

#--------------------------------------------------------------------------
print("Loading external libraries.",flush = True)
# --------------------------------------------------------------------------
import os, sys
from scipy import stats
from statistics import mean
import math as math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import argparse
import random

# --------------------------------------------------------------------------
print("Reading commmandline input with optparse.", flush = True)
# --------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="This script runs a random forest test on various datasets.")
# parser.add_option("-f", "--file", dest="filename",
#                   help="write report to FILE", metavar="FILE")
parser.add_argument("-d", "--homedir",
                  default=os.path.expanduser(os.path.join("~", "git", "lognorm_vs_CODA")),
                  help="path to git balance treee exploration git repository", dest="homedir", metavar="homedir")
options, unknown = parser.parse_known_args()

# --------------------------------------------------------------------------
print("Establishing directory layout.", flush = True)
# --------------------------------------------------------------------------
home_dir = os.path.expanduser(options.homedir)
projects = ["Jones", "Vangay", "Zeller", "Noguera-Julian"]
output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(output_dir)
plot_pdf_fpath = os.path.join(output_dir, "all_summary_pval_pub.pdf")
plotdata_fpath = os.path.join(home_dir,"metastudies","output","all_summary_pvalue_plot_pub.csv")
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
ds_color = {
		"rarefied_DADA2" : "plum",
	'alr_DADA2' : 'white',
	"clr_DADA2" : 'white',
	"raw_DADA2" : 'white',
			'lognorm_DADA2' : "lime",
			# "lognorm_Silva_DADA2" : "lime",
			'Silva_DADA2' : 'white',
			'Silva_DADA2_blw.sqrt_enorm' : '#050598',
			'Shuffle1_PhILR_Silva_DADA2_blw.sqrt_enorm' : '#f7d8a0',
			'Shuffle2_PhILR_Silva_DADA2_blw.sqrt_enorm' : '#f7d8a0',
			'Shuffle3_PhILR_Silva_DADA2_blw.sqrt_enorm' : '#f7d8a0',
			'Filtered_Silva_DADA2' : 'white',
			'Filtered_Silva_DADA2_blw.sqrt_enorm' : '#050598',
			'Shuffle1_PhILR_Filtered_Silva_DADA2_blw.sqrt_enorm' : '#f7d8a0',
			'Shuffle2_PhILR_Filtered_Silva_DADA2_blw.sqrt_enorm' : '#f7d8a0',
			'Shuffle3_PhILR_Filtered_Silva_DADA2_blw.sqrt_enorm' : '#f7d8a0',
	    }

comp_ds = list(ds_color.keys())
print(comp_ds)
my_colors = list(ds_color.values())

train_percent = 0.75
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_pdf_fpath)
#set font sizes
plt.rc('font', size=15) 
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('axes', labelsize=20) 
plt.rc('axes', titlesize=50)
# --------------------------------------------------------------------------
print("Generating Data.", flush = True)
# --------------------------------------------------------------------------
plotdata = pd.DataFrame(columns=comp_ds, index=comp_ds)
for ds1 in comp_ds:
	for ds2 in comp_ds:
		if (ds1 != ds2):
			for project in projects:
				if (ds1 != ds2):
					# print(f"{ds1} {ds2}")
					ds1_means = []
					ds2_means = []
					for project in projects:
						# print(f"Adding project {project}")
						op_dir = os.path.join(home_dir, project, "output")
						result_fpath = os.path.join(op_dir, "tables", f"sklearn_random_forest_manual_{train_percent}train.csv")
						# print(result_fpath)
						my_table = pd.read_csv(result_fpath, sep=',', header=0)
						#table 1
						ds1_table = my_table.loc[my_table["dataset"] == ds1,]
						splits = ds1_table.columns[ds1_table.columns.str.startswith('split')].tolist()
						ds1_means.extend(ds1_table[splits].agg(mean, axis = 1).values)
						ds2_table = my_table.loc[my_table["dataset"] == ds2,]
						ds2_means.extend(ds2_table[splits].agg(mean, axis = 1).values)
					my_tpval = stats.wilcoxon(ds1_means, ds2_means,).pvalue
					ave_diff = sum(ds1_means) - sum(ds2_means) #ds1_ave - ds2_ave
					if ave_diff > 0:
						plotdata.loc[ds1,ds2] = math.log10(my_tpval)
					else:
						plotdata.loc[ds1,ds2] = -math.log10(my_tpval)
				else:
					plotdata.loc[ds1,ds2] = 0
plotdata.replace("propotion", "proportion")# because I mispelled it elsehere
#--------------------------------------------------------------------------
print("Generating graphic")
#--------------------------------------------------------------------------
fig = plt.figure(figsize=(16,12))
fig.suptitle(f"Metastudy {train_percent}training each dataset vs others by accuracy, Sklearn RF")
plt.subplots_adjust(bottom=0.8, left=0.8)
ax = fig.add_subplot(1,1,1)
bp = ax.boxplot(plotdata, patch_artist = True, labels=plotdata.columns,showfliers=False)
for patch, color in zip(bp['boxes'], my_colors):
	patch.set_facecolor(color)
ax.set_xticklabels(labels = plotdata.columns, rotation=90)
# plt.annotate(label, (x_lst[i], y_lst[i]))
plt.axhline(y = math.log10(0.05), color = 'r', label="Significantly worse accuracy below")
plt.axhline(y = 0, color = 'y', label="no difference", linestyle="--")
plt.axhline(y = -math.log10(0.05), color = 'g', label="Significantly better accuracy above")
# ax.fill_between(x=[0,], y1=math.log10(0.5), y2=5, color="green", alpha=0.1)
# ax.add_patch(Rectangle(xy=[0,math.log10(0.5)], width=len(plotdata.columns), height=5, facecolor = "green", fill=True))
# ax.add_patch(Rectangle(xy=[0,-math.log10(0.5)], width=len(plotdata.columns), height=-5, facecolor = "red", fill=True))
# plt.axvline(x=0, color='r', label="No difference", linestyle="--")
ax.set_xlabel(f"Wilcoxen pairwise pvalues")
ax.set_ylabel(f"log10 pvalue")
# ax.legend(loc="upper center", framealpha=0.1, prop={'size': 8})
for i in range(len(plotdata.columns)):
	y = plotdata.iloc[:,i]
	x = np.random.normal(1+i, 0.04, size=len(y))
	ax.plot(x, y, color="b",marker=".", linestyle = "None", alpha=0.5)
fig.tight_layout()
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )

# print("Making seperate legend.")
# fig.suptitle(f"Metastudy {train_percent}training {ds1} vs others by accuracy, Python only")
# ax = fig.add_subplot(1,1,1)
# for i in range(len(comp_ds)):
# 	ax.scatter(0, 0, s=70, label=comp_ds[i], color="black", marker=my_markers[i])
# for i in range(len(projects)):
# 	ax.scatter(0, 0, s=70, label=projects[i], color=my_colors, marker=my_markers[0])
# plt.axvline(x=0, color='r', label="No difference", linestyle="--")
# plt.axhline(y = math.log10(0.1), color = 'r', label="p=0.10")
# ax.legend(title="Legend", loc="center", mode="expand", framealpha=1)
# fig.tight_layout()
# print(f"Saving figure to pdf at {plot_pdf_fpath}", flush = True)
# pdf.savefig( fig )

print("Saving pdf", flush = True)
pdf.close()

plotdata.to_csv(plotdata_fpath)

print(f"{__file__} complete!")
