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
parser.add_argument("-i", "--input_file_paradigm",
                  default="sklrn_randmforst_manual_0.75train_variable_ntree.csv",
                  help="""
									The file name that all of the input files share. They should all be in their
									perspective output/tables/ folders. 
									sklrn_randmforst_manual_0.75train_variable_ntree.csv is the default, the R
									version will be wide_random_forest_score_R_var_ntree20.csv.
									""", dest="input_file_paradigm", metavar="input_file_paradigm")
parser.add_argument("-l", "--label", default="sklearn_rf", help="Label for chart")
parser.add_argument("-c", "--col_name_line", default="dataset", help="Column name for lines")
parser.add_argument("-p", "--col_name_plot", default="metadata", help="Column name for each plot")

options, unknown = parser.parse_known_args()

# --------------------------------------------------------------------------
print("Establishing directory layout.", flush = True)
# --------------------------------------------------------------------------
home_dir = os.path.expanduser(options.homedir)
projects = ["Jones", "Vangay", "Zeller", "Noguera-Julian"]
output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(output_dir)
plot_pdf_fpath = os.path.join(output_dir, f"{options.label}_num_estimators_vs_accuracy_py.pdf")
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
fig.suptitle(f"Metastudy num estimators vs accuracy {options.label}")
plt.subplots_adjust(bottom=0.8)
ax = fig.add_subplot(1,1,1)
for project in projects:
	my_results = {}#use structure {feature: [n_trees, mean_accuracy]}
	print(f"Adding project {project}")
	op_dir = os.path.join(home_dir, project, "output")
	result_fpath = os.path.join(op_dir, "tables", options.input_file_paradigm)
	print(result_fpath)
	my_table = pd.read_csv(result_fpath, sep=',', header=0)
	my_table = my_table[my_table[options.col_name_line] =="Raw_DADA2"]
	splits = my_table.columns[my_table.columns.str.startswith('split')].tolist()
	my_marker = my_markers[projects.index(project)]
	for meta_d in set(my_table[options.col_name_plot].values):
		new_table = my_table[my_table[options.col_name_plot]==meta_d]
		means = new_table[splits].agg(mean, axis = 1)
		assert not means.empty, f"is not in the table from {project}"
		my_labels = [f"{project}_{feat}" for feat in new_table[options.col_name_plot]]
		ax.scatter(new_table["ntrees"], means, s=70, label=meta_d, marker=my_marker)
		ax.plot(new_table["ntrees"], means)
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
	result_fpath = os.path.join(op_dir, "tables", options.input_file_paradigm)
	print(result_fpath)
	my_table = pd.read_csv(result_fpath, sep=',', header=0)
	my_table = my_table[my_table[options.col_name_line] =="Raw_DADA2"]
	splits = my_table.columns[my_table.columns.str.startswith('split')].tolist()
	my_marker = my_markers[projects.index(project)]
	for meta_d in set(my_table[options.col_name_plot].values):
		new_table = my_table[my_table[options.col_name_plot]==meta_d]
		means = new_table[splits].agg(mean, axis = 1)
		assert not means.empty, f"is not in the table from {project}"
		my_labels = [f"{project}_{feat}" for feat in new_table[options.col_name_plot]]
		ax.scatter(new_table["ntrees"], means, s=70, label=meta_d, marker=my_marker)
		ax.plot(new_table["ntrees"], means)
ax.set_xlabel("Number of estimators")
ax.set_ylabel("Accuracy")
fig.tight_layout()
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )
ax.set_xlabel("Number of estimators")
ax.set_ylabel("Accuracy")
ax.legend(title="Legend", loc="center", mode = "expand", framealpha=1)
fig.tight_layout()
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )

print("Saving pdf", flush = True)
pdf.close()