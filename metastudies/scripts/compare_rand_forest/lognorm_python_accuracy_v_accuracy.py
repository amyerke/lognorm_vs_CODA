#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
# This is a script for comparing random forest output to pvalues
# The solution for non-repeating colors: https://stackoverflow.com/questions/53199728/how-can-i-stop-matplotlib-from-repeating-colors
# ds1_score data architecture: ds1_score[f"{project}_{feat}"] = [ave, project]

print(f"""Running {__file__}.
This is a script for comparing random forest output with pvalues.
""")
# --------------------------------------------------------------------------
print("Loading external libraries.",flush = True)
# --------------------------------------------------------------------------
import os
from statistics import mean
from random import sample, seed
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from scipy.stats import linregress
import argparse

# --------------------------------------------------------------------------
print("Reading commmandline input with optparse.", flush = True)
# --------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="This script runs a random forest test on various datasets.")
# parser.add_option("-f", "--file", dest="filename",
#                   help="write report to FILE", metavar="FILE")
parser.add_argument("-d", "--homedir",
                  default=os.path.expanduser(os.path.join("~", "git", "balance_tree_exploration")),
                  help="path to git balance treee exploration git repository", dest="homedir", metavar="homedir")
options, unknown = parser.parse_known_args()

# --------------------------------------------------------------------------
print("Establishing directory layout.", flush = True)
# --------------------------------------------------------------------------
home_dir = os.path.expanduser(options.homedir)
projects = ["Jones", "Vangay", "Zeller", "Noguera-Julian"]
output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(output_dir)
plot_pdf_fpath = os.path.join(output_dir, "lognorm_accuracy_vs_accuracy_python_by_transformation.pdf")
r_sq_boxplot_pdf_fpath = os.path.join(output_dir, "lognorm_rsq_accuracy_v_accuracy_py_by_transf.pdf")
train_percent = 0.75
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
font1 = {'family':'serif','color':'blue','size':20}
font2 = {'family':'serif','color':'darkred','size':15}
# comp_ds = ['alr_DADA2', 'clr_DADA2', 'raw_DADA2', 'Filtered_IQtree_blw.sqrt_enorm', \
# 	"lognorm_Silva_DADA2", 'Filtered_Silva_DADA2', 'Filtered_Silva_DADA2_blw.sqrt_enorm', \
# 	'Filtered_UPGMA_DADA2_blw.sqrt_enorm', 'Silva_DADA2', 'Silva_DADA2_blw.sqrt_enorm']
comp_ds = ['alr_DADA2', 'clr_DADA2', 'raw_DADA2', 'Filtered_IQtree_blw.sqrt_enorm', \
	'Filtered_Silva_DADA2', 'Filtered_Silva_DADA2_blw.sqrt_enorm', \
	'Filtered_UPGMA_DADA2_blw.sqrt_enorm', 'Silva_DADA2', 'Silva_DADA2_blw.sqrt_enorm']


pdf = matplotlib.backends.backend_pdf.PdfPages(plot_pdf_fpath)
#set font sizes
plt.rc('font', size=15) 
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('axes', labelsize=17) 
plt.rc('axes', titlesize=25)
# plt.rc('image', cmap='tab20c')
# plt.set_cmap("tab20")
r_sq = []

fig = plt.figure(figsize=(15,15))
num_rows = 3
num_cols = 3
ax_count = 1
my_markers = ["$A$","$B$","$C$","$D$","$E$","$F$","$G$","$H$","$I$","$J$","$K$","$L$",\
	"$M$","$N$","$O$","$P$","$Q$","$R$","$S$","$T$","$U$","$V$","$W$","$X$","$Y$","$Z$"]
my_colors = ["green", "blue", "orange", "red"]
ln_folder = os.path.join(output_dir, "lognormVs")#place for tables of ds2s vs lognorm
ds1 = "lognorm_DADA2"
for d2 in range(len(comp_ds)):
	ds2 = comp_ds[d2]
	train_percent = 0.75
	ds1_score = []
	ds1_feature = []
	ds1_project = []
	ds1_proj_feat = []
	ds2_score = []
	ds2_feature = []
	ds2_project = []
	ds2_proj_feat = []
	num_colors = list()
	print(f"Running datasets: {ds1} and {ds2}")
	for project in projects:
		print(f"Adding project {project}")
		main_output_label = f"sklearn_random_forest_manual_{train_percent}train"
		op_dir = os.path.join(home_dir, project, "output")
		result_fpath = os.path.join(op_dir, "tables", f"{main_output_label}.csv")
		print(result_fpath)
		my_table = pd.read_csv(result_fpath, sep=',', header=0)
		#table 1
		ds1_table = my_table.loc[my_table["dataset"] == ds1,]
		splits = ds1_table.columns[ds1_table.columns.str.startswith('split')].tolist()
		means = ds1_table[splits].agg(mean, axis = 1)
		num_colors.append(len(means))#stores the num features per proj -used to prevent repeat colors in plot
		assert not means.empty, f"{ds1} is not in the table from {project}"
		for feat, ave in zip(list(ds1_table["metadata"].values) ,list(means)):
			ds1_feature.append(feat)
			ds1_project.append(project)
			ds1_proj_feat.append(f"{project}_{feat}")
			if ave >= 0:
				ds1_score.append(ave)
			else:
				ds1_score.append(0)
		#table 2
		ds2_table = my_table.loc[my_table["dataset"] == ds2,]
		splits = ds2_table.columns[ds2_table.columns.str.startswith('split')].tolist()
		means = ds2_table[splits].agg(mean, axis = 1)
		assert not means.empty, f"{ds2} is not in the table from {project}"
		# print(my_table)
		for feat, ave in zip(list(ds2_table["metadata"].values) ,list(means)):
			ds2_feature.append(feat)
			ds2_project.append(project)
			ds2_proj_feat.append(f"{project}_{feat}")
			if ave >= 0:
				ds2_score.append(ave)
			else:
				ds2_score.append(0)
	assert ds1_proj_feat == ds2_proj_feat, f"{ds1} and {ds2} features are not the same"
	assert ds1_feature == ds2_feature, f"{ds1} and {ds2} features are not the same"
	print(f"Saving table to 'output_dir/lognormVs'.")
	pre_df_zip = list(zip(ds1_proj_feat, ds1_feature,ds1_score, ds1_project, ds2_feature, ds2_score, ds2_project))
	my_df = pd.DataFrame(pre_df_zip, columns=["proj_feat", "ln_feature","ln_score","ln_project",f"{ds2}_feature", f"{ds2}_score", f"{ds2}_project"])
	if not os.path.exists(ln_folder):
		os.makedirs(ln_folder)
	my_df.to_csv(os.path.join(ln_folder, f"lnVs_{ds2}.csv"), index=False)
	print("Building graphic.")
	slope, intercept, r_value, p_value, std_err = linregress(ds1_score, ds2_score)
	r_sq.append(r_value**2)
	a, b = np.polyfit(ds1_score, ds2_score, 1)
	# fig.suptitle(f"Metastudy {train_percent}training {ds1} vs {ds2} by accuracy, Python only")
	plt.subplots_adjust(bottom=0.8)
	ax = fig.add_subplot(num_rows,num_cols, ax_count)
	my_projects = list(set(ds1_project))
	flag_old_label = ds1_project[0]
	feature_counter = 0
	for i in range(len(ds1_score)):
		my_proj = ds1_project[i]
		if flag_old_label != my_proj:
			flag_old_label = my_proj
			feature_counter = 0
		# print(f"{my_proj} {ds1_score[i]} {ds2_score[i]}, {ds1_proj_feat[i]}, {ds2_proj_feat[i]}, featcount: {feature_counter}")
		my_color = my_colors[my_projects.index(my_proj)]
		ax.scatter(ds1_score[i], ds2_score[i], s=100, color=my_color, label=ds1_feature[i], marker=my_markers[feature_counter])
		# ax.text(ds1_score[i], ds2_score[i], ds1_proj_feat[i], color=my_color, fontsize=5)
		feature_counter += 1
	ax.locator_params(axis='both', nbins=10)
	ax.plot([0,1], [0,1], color = "r", label = "expected")
	ax.plot(np.array(ds1_score), a*np.array(ds1_score)+b, color = "green", label = "accuracy polyfit")
	ax.text(0.1, 0.9, f"r squared: {round(r_value**2, 4)}", fontsize=20)
	ax.set_xlabel(f"{ds1}")
	ax.set_ylabel(f"{ds2}")
	ax.set_title("Scores")
	ax_count += 1
	# plt.set_cmap("tab20")
	# if ax_count > 9:
		# # fig.tight_layout(pad = 1, h_pad=1, w_pad=2)
		# # plt.subplots_adjust(left=0.1)
		# # print("Saving figure to pdf", flush = True)
		# # pdf.savefig( fig )
		# print("Making seperate legend.")
		# fig = plt.figure(figsize=(11,16))
		# ax = fig.add_subplot(1,1,1)
		# ax.legend(title="",  loc="center", framealpha=1, mode = "expand", markerscale=2)
		# print("Saving figure to pdf", flush = True)
		# pdf.savefig( fig )
	# ax.legend(title="Legend", loc="upper left", framealpha=1)
fig.tight_layout(pad = 1, h_pad=1, w_pad=2)
plt.subplots_adjust(left=0.1)
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )

print("Making seperate legend.")
fig = plt.figure(figsize=(11,25))
ax = fig.add_subplot(1,1,1)
flag_old_label = ds1_project[0]
feature_counter = 0
print("building graphic")
for i in range(len(ds1_score)):
	my_proj = ds1_project[i]
	if flag_old_label != my_proj:
		flag_old_label = my_proj
		feature_counter = 0
	# print(f"{my_proj} {ds1_score[i]} {ds2_score[i]}, {ds1_proj_feat[i]}, {ds2_proj_feat[i]}, featcount: {feature_counter}")
	my_color = my_colors[my_projects.index(my_proj)]
	ax.scatter(ds1_score[i], ds2_score[i], s=100, color=my_color, label=ds1_feature[i], marker=my_markers[feature_counter])
	ax.text(ds1_score[i], ds2_score[i], ds1_feature[i], color=my_color, fontsize=5)
	ax.legend(title="",  loc="center", framealpha=1, mode = "expand", markerscale=2)
	feature_counter += 1
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )
print(f"Saving pdf to {plot_pdf_fpath}", flush = True)
pdf.close()

print(f"Accuracy vs accuracy R squared plots")
pdf = matplotlib.backends.backend_pdf.PdfPages(r_sq_boxplot_pdf_fpath)
fig = plt.figure(figsize=(11,16))
fig.suptitle(f"Accuracy vs accuracy R squared")
plt.subplots_adjust(bottom=0.8)
ax = fig.add_subplot(1,1,1)
bp = ax.boxplot(r_sq)
fig.tight_layout()
pdf.savefig( fig, bbox_inches='tight')

print(f"Saving pdf to {r_sq_boxplot_pdf_fpath}", flush = True)
pdf.close()

first_df = "placeholder"
print(f"Joining tables from {ln_folder}.")
for d2 in range(len(comp_ds)):
	print(f"counter = {d2}")
	ds2 = comp_ds[d2]
	print(f"DS2: {ds2}")
	file_pth = os.path.join(ln_folder, f"lnVs_{ds2}.csv")
	if d2 == 0:
		first_df = pd.read_csv(file_pth)
		# print(first_df)
	else:
		my_df = pd.read_csv(file_pth)
		my_cols = [s for s in my_df.columns if "ln_" in s]
		my_df = my_df.drop(my_cols, axis=1)
		print(my_df)
		first_df = pd.merge(first_df, my_df, on="proj_feat")

first_df.to_csv(os.path.join(ln_folder, f"final_ln.csv"), index=False)
# print(first_df)
my_cols = [s for s in first_df.columns if "_score" in s]
my_cols.insert(0, "proj_feat")
print(my_cols)
first_df = first_df.filter(my_cols)
print(first_df)
first_df.to_csv(os.path.join(ln_folder, f"final_ln_scoreOnly.csv"), index=False)

print(f"{__file__} complete!")