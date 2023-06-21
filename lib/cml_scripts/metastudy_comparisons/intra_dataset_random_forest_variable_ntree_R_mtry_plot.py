#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
# This is a script for comparing random forest output to pvalues
# --------------------------------------------------------------------------
print(f"Running {__file__}")
print("""This is a script for comparing random forest accuracy output to asv pvalues.""")
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
print("Loading external libraries.",flush = True)
# --------------------------------------------------------------------------
import os, sys
import time
from statistics import mean
from matplotlib import markers
import math as math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from sklearn.metrics import accuracy_score, roc_auc_score
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
plot_pdf_fpath = os.path.join(output_dir, "num_estimators_vs_accuracy_R_mtry.pdf")
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
comp_ds = ['alr_DADA2', 'clr_DADA2', 'Raw_DADA2', 'Filtered_IQtree', \
	'Filtered_IQtree_mean.descendants_enorm', 'Filtered_Silva_DADA2', \
	'Filtered_Silva_DADA2_mean.descendants_enorm', 'Filtered_UPGMA_DADA2', \
	'Filtered_UPGMA_DADA2_mean.descendants_enorm', 'lognorm_DADA2', 'Silva_DADA2', \
	'Silva_DADA2_mean.descendants_enorm']

my_markers = ["o", "s", "P", "v", "x"]

# --------------------------------------------------------------------------
print("Setting text sizes.", flush = True)
# --------------------------------------------------------------------------
plt.rc('font', size=15) 
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('axes', labelsize=20) 
plt.rc('axes', titlesize=50)
# --------------------------------------------------------------------------
print("Main loop.", flush = True)
# --------------------------------------------------------------------------
pdf = matplotlib.backends.backend_pdf.PdfPages(plot_pdf_fpath)
fig = plt.figure(figsize=(11,11))
fig.suptitle(f"Metastudy num estimators vs accuracys")
plt.subplots_adjust(bottom=0.8)
ax = fig.add_subplot(1,1,1)
for project in projects:
	my_results = {}#use structure {feature: [n_trees, mean_accuracy]}
	print(f"Adding project {project}")
	op_dir = os.path.join(home_dir, project, "output")
	result_fpath = os.path.join(op_dir, "tables", f"sklrn_randmforst_manual_0.75train_variable_ntree_R_mtry.csv")
	print(result_fpath)
	my_table = pd.read_csv(result_fpath, sep=',', header=0)
	my_table = my_table[my_table["dataset"]=="Raw_DADA2"]
	splits = my_table.columns[my_table.columns.str.startswith('split')].tolist()
	my_marker = my_markers[projects.index(project)]
	for meta_d in set(my_table["metadata"].values):
		new_table = my_table[my_table["metadata"]==meta_d]
		means = new_table[splits].agg(mean, axis = 1)
		assert not means.empty, f"is not in the table from {project}"
		my_labels = [f"{project}_{feat}" for feat in new_table["metadata"]]
		ax.scatter(new_table["n_trees"], means, s=70, label=meta_d, marker=my_marker)
		ax.plot(new_table["n_trees"], means)
ax.set_xlabel("Number of estimators")
ax.set_ylabel("Accuracy")
fig.tight_layout()
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )

# --------------------------------------------------------------------------
print("Making legend.", flush = True)
# --------------------------------------------------------------------------
# simply duplicating the above code to make legend
for project in projects:
	my_results = {}#use structure {feature: [n_trees, mean_accuracy]}
	print(f"Adding project {project}")
	op_dir = os.path.join(home_dir, project, "output")
	result_fpath = os.path.join(op_dir, "tables", f"sklrn_randmforst_manual_0.75train_variable_ntree_R_mtry.csv")
	print(result_fpath)
	my_table = pd.read_csv(result_fpath, sep=',', header=0)
	my_table = my_table[my_table["dataset"]=="Raw_DADA2"]
	splits = my_table.columns[my_table.columns.str.startswith('split')].tolist()
	my_marker = my_markers[projects.index(project)]
	for meta_d in set(my_table["metadata"].values):
		new_table = my_table[my_table["metadata"]==meta_d]
		means = new_table[splits].agg(mean, axis = 1)
		assert not means.empty, f"is not in the table from {project}"
		my_labels = [f"{project}_{feat}" for feat in new_table["metadata"]]
		ax.scatter(new_table["n_trees"].values, means, s=70)
		# ax.plot(new_table["n_trees"], means)
ax.set_xlabel("Number of estimators")
ax.set_ylabel("Accuracy")
ax.legend(title="Legend", loc="center", mode = "expand", framealpha=1)
fig.tight_layout()
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )

print("Saving pdf", flush = True)
pdf.close()
