#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
#This is a script for comparing R and Python random forest output
# --------------------------------------------------------------------------
print(f"""Running {__file__}.
This is a script for comparing R and Python random forest output. 
Our hypthesis is that R does slightly better than Python.""")
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
parser.add_argument("-x", "--x_axis_fname", default=None,
            			help="File name of the csv in the \{project\}/output/tables section in each project for the x axis")
parser.add_argument("-y", "--y_axis_fname", default=None,
            			help="File name of the csv in the \{project\}/output/tables section in each project for the y axis")
parser.add_argument("-a", "--x_axis_col_name", help="column name from the x axis csv", default="trans_group")
parser.add_argument("-b", "--y_axis_col_name", help="column name from the y axis csv", default="dataset")
parser.add_argument("-t", "--x_axis_metadata_col_name", help="metadata column name from the x axis csv", default="metadata_col")
parser.add_argument("-u", "--y_axis_metadata_col_name", help="metadata column name from the y axis csv", default="metadata")
parser.add_argument('-p', '--project_list', nargs='+', default=[],
										help="Project folders that have the x and y files.")
options, unknown = parser.parse_known_args()
parser.add_argument("-v", "--x_axis_label", help="label for the x axis", default=os.path.basename(options.x_axis_fname))
parser.add_argument("-w", "--y_axis_label", help="label for the x axis", default=os.path.basename(options.y_axis_fname))
options, unknown = parser.parse_known_args()
# --------------------------------------------------------------------------
print("Establishing directory layout.", flush = True)
# --------------------------------------------------------------------------
home_dir = os.path.expanduser(options.homedir)
projects = options.project_list
output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(output_dir)
plot_pdf_fpath = os.path.join(output_dir, f"inter_ds_comp_{options.x_axis_label}_vs_{options.y_axis_label}.pdf")
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
comp_ds = ['alr_DADA2', 'clr_DADA2', 'Raw_DADA2', 'Filtered_IQtree', \
	'Filtered_IQtree_blw.sqrt_enorm', 'Filtered_Silva_DADA2', \
	'Filtered_Silva_DADA2_blw.sqrt_enorm', 'Filtered_UPGMA_DADA2', \
	'Filtered_UPGMA_DADA2_blw.sqrt_enorm', 'lognorm_DADA2', 'Silva_DADA2', \
	'Silva_DADA2_blw.sqrt_enorm']

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
	x_score = {}
	y_score = {}
	features = []
	proj = []
	print(f"{ds}")
	for project in projects:
		print(f"Adding project {project}")
		op_dir = os.path.join(home_dir, project, "output")
		result_fpath = os.path.join(op_dir, "tables", options.x_axis_fname)
		print(result_fpath)
		my_table = pd.read_csv(result_fpath, sep=',', header=0)
		#table 1
		# print(f"{result_fpath}")
		my_table = my_table.loc[my_table[options.x_axis_col_name] == ds,]
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
		assert not means.empty, f"{ds} is not in the {options.x_label} table from {project}"
		for feat, ave in zip(list(my_table[options.x_axis_metadata_col_name].values) ,list(means)):
			print(f"{options.x_axis_label}: {ave} for {feat}")
			x_score[f"{project}_{feat}"] = [ave, project]
		#table 2
		op_dir = os.path.join(home_dir, project, "output")
		result_fpath = os.path.join(op_dir, "tables", options.y_axis_fname)
		print(f"{options.y_axis_label} tables")
		my_table = pd.read_csv(result_fpath, sep=',', header=0)
		# print(my_table.to_string())
		my_table = my_table.loc[my_table[options.y_axis_col_name] == ds,]
		# print(my_table.to_string())
		splits = my_table.columns[my_table.columns.str.startswith('split')].tolist()
		means = my_table[splits].agg(mean, axis = 1)
		# assert not means.empty, f"{ds} is not in the table from {project}"
		for feat, ave in zip(list(my_table[options.y_axis_metadata_col_name].values) ,list(means)):
			print(f"{options.y_axis_label}: {ave} for {feat}")
			y_score[f"{project}_{feat}"] = [ave, project]
			# pvals[feat] = min(my_table.loc[(my_table["meta_name"]==feat) & (my_table["taxa_lev"] == "Genus"), "pval"].values) + 1e-100
	same_keys = set(y_score.keys()).intersection(set(x_score.keys()))
	x_score = {key:x_score[key] for key in same_keys}
	y_score = {key:y_score[key] for key in same_keys}
	x_lst = np.array(list(map(lambda x: x[0], list(x_score.values()))))
	y_lst = np.array(list(map(lambda x: x[0], list(y_score.values()))))
	if len(x_lst) > 0:
		print("build graphic")
		fig = plt.figure(figsize=(11,11))
		fig.suptitle(f"Metastudy {train_percent}training {ds} {options.x_axis_label} vs {options.y_axis_label}")
		plt.subplots_adjust(bottom=0.8)
		ax = fig.add_subplot(1,1,1)
		# ax.scatter(x_lst, y_lst, label=list(y_score.keys()))
		ax.plot([0,1], [0,1], color = "r", label = "expected")
		# ax.plot(x_lst, a*x_lst+b, color = "green", label = "accuracy polyfit")
		ax.set_xlabel(f"{options.x_axis_label} accuracy")
		ax.set_ylabel(f"{options.y_axis_label} accuracy")
		fig.tight_layout()
		my_projects = list(set(list(map(lambda x: x[1], list(y_score.values())))))#pulling second element from each dict.value
		my_markers = ["o", "s", "P", "v", "x"]
		for i, label in enumerate(list(y_score.keys())):
			my_proj = y_score[label][1]
			my_marker = my_markers[my_projects.index(my_proj)]
			ax.scatter(x_lst[i], y_lst[i], s=70, label=list(y_score.keys())[i], marker=my_marker)
			# plt.annotate(label, (x_lst[i], y_lst[i]))
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
fig.tight_layout()
my_projects = list(set(list(map(lambda x: x[1], list(y_score.values())))))
for i, label in enumerate(list(y_score.keys())):
	my_proj = y_score[label][1]
	my_marker = my_markers[my_projects.index(my_proj)]
	ax.scatter(0, 0, s=70, label=list(y_score.keys())[i], marker=my_marker)
	# plt.annotate(label, (x_lst[i], y_lst[i]))
ax.legend(title="Legend",  loc="center", framealpha=1, mode = "expand", markerscale=2)
print("Saving figure to pdf", flush = True)
pdf.savefig( fig, bbox_inches='tight' )

print("Saving pdf", flush = True)
pdf.close()
