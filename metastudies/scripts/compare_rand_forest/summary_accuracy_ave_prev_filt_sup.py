#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com

#This is a script for comparing different transformations 
# --------------------------------------------------------------------------
print(f"""Running {__file__}.
This is a script for comparing random forest output with averages. 
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
plot_pdf_fpath = os.path.join(output_dir, "prev_filt_summary_ave_acc.pdf")
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
ds_color = {
	"rarefied_prev_filt_DADA2" : "plum",
	"alr_prev_filt_DADA2" : 'white',
	"clr_prev_filt_DADA2" : 'white',
	"raw_prev_filt_DADA2" : 'white',
			"propotions_prev_filt_DADA2" : "lime",
			"Heilinger_prev_filt_DADA2" : "lime",
			"lognorm_prev_filt_DADA2" : "lime",
			"Silva_prev_filt_DADA2" : 'white',
			"prev_filt90_Silva_DADA2_blw.sqrt_enorm" : '#050598',
			"Shuffle1_PhILR_prev_filt90_Silva_DADA2_blw.sqrt_enorm" : '#f7d8a0',
			"Shuffle2_PhILR_prev_filt90_Silva_DADA2_blw.sqrt_enorm" : '#f7d8a0',
			"Shuffle3_PhILR_prev_filt90_Silva_DADA2_blw.sqrt_enorm" : '#f7d8a0',
			"UPGMA_prev_filt_DADA2" : 'white',
			"prev_filt90_UPGMA_DADA2_blw.sqrt_enorm" : '#050598',
			"Shuffle1_PhILR_prev_filt90_UPGMA_DADA2_blw.sqrt_enorm" : '#f7d8a0',
			"Shuffle2_PhILR_prev_filt90_UPGMA_DADA2_blw.sqrt_enorm" : '#f7d8a0',
			"Shuffle3_PhILR_prev_filt90_UPGMA_DADA2_blw.sqrt_enorm" : '#f7d8a0',
			"IQtree_prev_filt_DADA2" : 'white',
			"prev_filt90_IQtree_blw.sqrt_enorm" : '#050598',
			"Shuffle1_PhILR_prev_filt90_IQtree_blw.sqrt_enorm" : '#f7d8a0',
			"Shuffle2_PhILR_prev_filt90_IQtree_blw.sqrt_enorm" : '#f7d8a0',
			"Shuffle3_PhILR_prev_filt90_IQtree_blw.sqrt_enorm" : '#f7d8a0',
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
median_props = {"color" : "red", "linewidth" : 3}
# --------------------------------------------------------------------------
print("Generating Data.", flush = True)
# --------------------------------------------------------------------------
all_means = {}
for ds1 in comp_ds:
	ds1_means = []
	print(f"Running {ds1}")
	for project in projects:
		# print(f"Adding project {project}")
		op_dir = os.path.join(home_dir, project, "output")
		result_fpath = os.path.join(op_dir, "tables", f"sklearn_random_forest_manual_{train_percent}train.csv")
		# print(result_fpath)
		my_table = pd.read_csv(result_fpath, sep=',', header=0)
		#table 1
		ds1_table = my_table.loc[my_table["dataset"] == ds1,]
		splits = ds1_table.columns[ds1_table.columns.str.startswith('split')].tolist()
		my_means = ds1_table[splits].agg(mean, axis = 1).values
		# my_means = [0 if x < 0 else x for x in my_means]
		ds1_means.extend(my_means)
	if len(ds1_means) < 1:
		print(f"There was a problem with {ds1}")
	all_means[ds1] = ds1_means
plotdata = pd.DataFrame(all_means)
plotdata.replace("propotion", "proportion")# because I mispelled it elsehere

print(f"My mean: {plotdata.mean()}")
#--------------------------------------------------------------------------
print("Generating graphic")
#--------------------------------------------------------------------------
fig = plt.figure(figsize=(24,17))
fig.suptitle(f"Metastudy {train_percent}training each dataset vs others by accuracy, Sklearn RF")
plt.subplots_adjust(bottom=0.8, left=0.8)
ax = fig.add_subplot(1,1,1)
# ax.boxplot(plotdata, labels=plotdata.columns, showfliers=False, medianprops=median_props)
ax.set_xticklabels(labels = plotdata.columns, rotation=90)
# plt.annotate(label, (x_lst[i], y_lst[i]))
ax.set_xlabel(f"Transformations")
ax.set_ylabel(f"Average accuracy for each metadata feature")
# ax.legend(loc="upper center", framealpha=0.1, prop={'size': 8})
bp = ax.boxplot(plotdata, patch_artist = True, labels=plotdata.columns,showfliers=False, medianprops=median_props)
for patch, color in zip(bp['boxes'], my_colors):
	patch.set_facecolor(color)
ax.axhline(y = plotdata.stack().median(), color = "brown", label="median")
fig.tight_layout()
print(f"Saving figure to pdf to {plot_pdf_fpath}", flush = True)
pdf.savefig( fig )

print(f"Closing pdf {plot_pdf_fpath}", flush = True)
pdf.close()

print(f"{__file__} complete!")
