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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import argparse
from scipy.stats import linregress

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
plot_pdf_fpath = os.path.join(output_dir, "r_sq_summary_pval_acc_vs_acc_python_by_transformation.pdf")
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
comp_ds = ['alr_DADA2', 'clr_DADA2', 'raw_DADA2', 'lognorm_DADA2', 'Silva_DADA2', \
	'Silva_DADA2_blw.sqrt_enorm', 'Shuffle1_PhILR_Silva_DADA2_blw.sqrt_enorm', \
	'Shuffle2_PhILR_Silva_DADA2_blw.sqrt_enorm', 'Shuffle3_PhILR_Silva_DADA2_blw.sqrt_enorm', \
	'Filtered_Silva_DADA2', 'Filtered_Silva_DADA2_blw.sqrt_enorm', \
	'Shuffle1_PhILR_Filtered_Silva_DADA2_blw.sqrt_enorm', \
	'Shuffle2_PhILR_Filtered_Silva_DADA2_blw.sqrt_enorm', \
	'Shuffle3_PhILR_Filtered_Silva_DADA2_blw.sqrt_enorm', \
	'Filtered_UPGMA_DADA2', 'Filtered_UPGMA_DADA2_blw.sqrt_enorm', \
	'Shuffle1_PhILR_Filtered_UPGMA_DADA2_blw.sqrt_enorm', \
	'Shuffle2_PhILR_Filtered_UPGMA_DADA2_blw.sqrt_enorm', \
	'Shuffle3_PhILR_Filtered_UPGMA_DADA2_blw.sqrt_enorm', \
	'Filtered_IQtree', 'Filtered_IQtree_blw.sqrt_enorm', \
	'Shuffle1_PhILR_Filtered_IQtree_blw.sqrt_enorm', \
	'Shuffle2_PhILR_Filtered_IQtree_blw.sqrt_enorm',\
	'Shuffle3_PhILR_Filtered_IQtree_blw.sqrt_enorm']

my_colors = ['white', 'white', 'white', "lime", 'white', '#050598', \
'#f7d8a0', '#f7d8a0', '#f7d8a0', 'white', '#050598', '#f7d8a0', \
'#f7d8a0', '#f7d8a0', 'white', '#050598', '#f7d8a0', '#f7d8a0', \
'#f7d8a0', 'white', '#050598', '#f7d8a0', '#f7d8a0', '#f7d8a0']

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
plotdata = pd.DataFrame(columns=comp_ds, index=comp_ds)
for ds1 in comp_ds:
	for ds2 in comp_ds:
		if (ds1 != ds2):
			for project in projects:
				if (ds1 != ds2):
					print(f"{ds1} {ds2}")
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
					slope, intercept, r_value, p_value, std_err = linregress(ds1_means, ds2_means)
					r_sq = r_value**2
					plotdata.loc[ds1,ds2] = r_sq
				else:
					plotdata.loc[ds1,ds2] = 0
#--------------------------------------------------------------------------
print("Generating graphic")
#--------------------------------------------------------------------------
fig = plt.figure(figsize=(11,11))
fig.suptitle(f"Metastudy {train_percent}training each dataset vs others by slope, Sklearn RF")
plt.subplots_adjust(bottom=0.8, left=0.8)
ax = fig.add_subplot(1,1,1)
ax.boxplot(plotdata, labels=plotdata.columns, showfliers=False, medianprops=median_props)
ax.set_xticklabels(labels = plotdata.columns, rotation=90)
# plt.annotate(label, (x_lst[i], y_lst[i]))
ax.set_xlabel(f"Tranformations")
ax.set_ylabel(f"Correlation coefficient")
ax.axhline(y = plotdata.stack().median(), color = "brown", label="median")
# ax.legend(loc="upper center", framealpha=0.1, prop={'size': 8})
for i in range(len(plotdata.columns)):
	y = plotdata.iloc[:,i]
	x = np.random.normal(1+i, 0.04, size=len(y))
	ax.plot(x, y, color=my_colors[i], marker=".", alpha=0.5)
fig.tight_layout()
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )

print("Saving pdf", flush = True)
pdf.close()

plotdata.to_csv(os.path.join(home_dir,"metastudies","output","summary_r_sq_plot.csv"))

print(f"{__file__} complete!")
