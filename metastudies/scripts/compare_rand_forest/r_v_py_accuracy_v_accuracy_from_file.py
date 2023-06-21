#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
#This is a script for comparing R and Python random forest output
# --------------------------------------------------------------------------
print(f"Running {__file__}")
print("""
This is a script for comparing R and Python random forest output. 
Our hypthesis is that R does slightly better than Python.
""")
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
print("Loading external libraries.",flush = True)
# --------------------------------------------------------------------------
from operator import indexOf
import os, sys
import time
from statistics import mean
from turtle import color
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
# projects = ["Jones", "Vangay", "Zeller"]
output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(output_dir)
plot_pdf_fpath = os.path.join(output_dir, "inter_group_comp_R_v_Py_from_files.pdf")
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
comp_ds = ['raw_DADA2']

# comp_ds = ['alr_DADA2', 'clr_DADA2', 'Raw_DADA2', 'Filtered_IQtree', 'Filtered_Silva_DADA2',
#  'Filtered_UPGMA_DADA2', 'lognorm_DADA2', 'Silva_DADA2']

pdf = matplotlib.backends.backend_pdf.PdfPages(plot_pdf_fpath)
#set font sizes
plt.rc('font', size=15) 
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('axes', labelsize=20) 
plt.rc('axes', titlesize=50)
plt.rc('figure',titlesize=20)

for ds in comp_ds:
	print(ds)
	train_percent = 0.75
	pvals = {}
	r_score = {}
	py_score = {}
	features = []
	proj = []
	print(f"{ds}")
	for project in projects:
		print(f"Adding project {project}")
		op_dir = os.path.join(home_dir, project, "output")
		result_fpath = os.path.join(op_dir, "tables", f"sklearn_random_forest_from_file.csv")
		print(result_fpath)
		my_table = pd.read_csv(result_fpath, sep=',', header=0)
		#table 1
		my_table = my_table.loc[my_table["dataset"] == ds,]
		# print(ds1_table)
		splits = my_table.columns[my_table.columns.str.startswith('split')].tolist()
		# print(my_table[splits].to_string())
		means = my_table[splits].agg(mean, axis = 1)
		print(means)
		# if any(pd.isna(my_table[splits])):
		# 	print("Nan found")
		# 	print(my_table[splits])
		# 	print(project)
		# 	sys.exit()
		assert not means.empty, f"{ds} is not in the py table from {project}"
		for feat, ave in zip(list(my_table["metadata"].values) ,list(means)):
			print(f"Python: {ave} for {feat}")
			py_score[f"{project}_{feat}"] = [ave, project]
		#table 2
		op_dir = os.path.join(home_dir, project, "output")
		result_fpath = os.path.join(op_dir, "tables", f"wide_from_file_random_forest_score_R.csv")
		print("R tables")
		my_table = pd.read_csv(result_fpath, sep=',', header=0)
		# print(my_table.to_string())
		my_table = my_table.loc[my_table["trans_group"] == ds,]
		# print(my_table.to_string())
		splits = my_table.columns[my_table.columns.str.startswith('split')].tolist()
		means = my_table[splits].agg(mean, axis = 1)
		# assert not means.empty, f"{ds} is not in the table from {project}"
		# pval_fpath = os.path.join(op_dir, "tables", f"{project}_pValuesUnivariate_taxaVmetadata.csv")
		# my_table = pd.read_csv(pval_fpath, sep=',', header=0)
		print("R means")
		print(means)
		for feat, ave in zip(list(my_table["metadata_col"].values) ,list(means)):
			print(f"R: {ave} for {feat}")
			r_score[f"{project}_{feat}"] = [ave, project]
			# pvals[feat] = min(my_table.loc[(my_table["meta_name"]==feat) & (my_table["taxa_lev"] == "Genus"), "pval"].values) + 1e-100
	same_keys = set(r_score.keys()).intersection(set(py_score.keys()))
	r_score = {key:r_score[key] for key in same_keys}
	py_score = {key:py_score[key] for key in same_keys}
	ds1_lst = np.array(list(map(lambda x: x[0], list(py_score.values()))))
	ds2_lst = np.array(list(map(lambda x: x[0], list(r_score.values()))))
	if len(ds1_lst) > 0:
		print("build graphic")
		# a, b = np.polyfit(ds1_lst, ds2_lst, 1)
		fig = plt.figure(figsize=(11,11))
		fig.suptitle(f"Metastudy {train_percent}training {ds} Py vs R, accuracy vs accuracy")
		plt.subplots_adjust(bottom=0.8)
		ax = fig.add_subplot(1,1,1)
		# ax.scatter(ds1_lst, ds2_lst, label=list(r_score.keys()))
		ax.plot([0,1], [0,1], color = "r", label = "expected")
		# ax.plot(ds1_lst, a*ds1_lst+b, color = "green", label = "accuracy polyfit")
		ax.set_xlabel("Accuracy Python Rand Forest")
		ax.set_ylabel("Accuracy R Rand Forest")
		fig.tight_layout()
		my_projects = list(set(list(map(lambda x: x[1], list(r_score.values())))))
		my_markers = ["o", "s", "P", "v", "x"]
		for i, label in enumerate(list(r_score.keys())):
			my_proj = r_score[label][1]
			my_marker = my_markers[my_projects.index(my_proj)]
			ax.scatter(ds1_lst[i], ds2_lst[i], s=70, label=list(r_score.keys())[i], marker=my_marker)
			# plt.annotate(label, (ds1_lst[i], ds2_lst[i]))
		# ax.legend(title="Legend", loc="best", framealpha=0.1)
		print("Saving figure to pdf", flush = True)
		pdf.savefig( fig, bbox_inches='tight' )
		# sys.exit()
		# print(my_table.columns)
		# print(sorted(set(my_table["dataset"].values),key=lambda s: s.lower()))X
print("Making seperate legend.")
fig = plt.figure(figsize=(11,11))
fig.suptitle(f"Legend {train_percent}training {ds} Py vs R, accuracy vs accuracy")
plt.subplots_adjust(bottom=0.8)
ax = fig.add_subplot(1,1,1)
ax.plot([0,1], [0,1], color = "r", label = "expected")
# ax.plot(ds1_lst, a*ds1_lst+b, color = "green", label = "accuracy polyfit")
ax.set_xlabel("Accuracy Python Rand Forest")
ax.set_ylabel("Accuracy R Rand Forest")
fig.tight_layout()
my_projects = list(set(list(map(lambda x: x[1], list(r_score.values())))))
my_markers = ["o", "s", "P", "v", "x"]
for i, label in enumerate(list(r_score.keys())):
	my_proj = r_score[label][1]
	my_marker = my_markers[my_projects.index(my_proj)]
	ax.scatter(ds1_lst[i], ds2_lst[i], s=70, label=list(r_score.keys())[i], marker=my_marker)
	# plt.annotate(label, (ds1_lst[i], ds2_lst[i]))
ax.legend(title="Legend",  loc="center", framealpha=1, mode = "expand", markerscale=2)
print("Saving figure to pdf", flush = True)
pdf.savefig( fig, bbox_inches='tight' )

print("Saving pdf", flush = True)
pdf.close()
