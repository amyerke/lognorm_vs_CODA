#!/usr/bin/env python3
print(f"""
Starting: {__file__}
Author: Aaron Yerke (aaronyerke@gmail.com)
For determining if PhILR weighting schemes improve any ML techniques.
This was useful: https://machinelearningmastery.com/compare-machine-learning-algorithms-python-scikit-learn/
Compare Algorithms
This script creates several plots and tables.
The first boxplot shows each algorithm by PhILR weighting scheme.
Then, this information is condensed into a table.
Another table is created to summerize the data in another way.
Finally, line plots of the means of the first plot are created to further condense the data.
""")

# --------------------------------------------------------------------------
print("Loading external libraries.")
# --------------------------------------------------------------------------

from cProfile import label
import math
import os, sys
from matplotlib import markers
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import matplotlib.colors as mcolors
from pandas.api.types import is_string_dtype
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
import argparse

# --------------------------------------------------------------------------
print("Running optparse.")
# --------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Process some integers.')
# parser.add_option("-f", "--file", dest="filename",
#                   help="write report to FILE", metavar="FILE")
parser.add_argument("-m", "--metadata_cols",
                  action="store_false", dest="meta_col",
                  help="Metadata columns to analyse")
parser.add_argument("-d", "--homedir",
                  default=os.path.expanduser(os.path.join("~", "git", "lognorm_vs_CODA")),
                  help="path to git balance treee exploration git repository", dest="homedir", metavar="homedir")
parser.add_argument("-p", "--project", default="string",
                  help="project folder", metavar="project")
parser.add_argument("-a", "--use_all_meta", default=False,
                  help="use all metadata", metavar="use_all_meta")
parser.add_argument("-f", "--metada_fn", default=False, dest="meta_fn",
                  help="Name of file at the top of the project folder to use as metadata.", 
									metavar="meta_fn")
parser.add_argument("-l", "--delimiter", default="\t",
                  help="File delimiting symbol for metadata. Default is tab.",
									metavar="delim", dest="delim")
parser.add_argument("-i", "--meta_index_col", default=0,
                  help="Name of column to use as row name for metadata",
                  metavar="meta_index_col", dest="meta_index_col")

options, unknown = parser.parse_known_args()

# --------------------------------------------------------------------------
print("Establishing directory layout.")
# --------------------------------------------------------------------------

home_dir = options.homedir
project = options.project
output_dir = os.path.join(home_dir, project, "output")
assert os.path.exists(output_dir)

# --------------------------------------------------------------------------
print("Establishing other constants.")
# --------------------------------------------------------------------------

output_label = "sklearn_ml_acc"
philr_groups = ["Silva_DADA2", "Filtered_Silva_DADA2", "Filtered_UPGMA_DADA2", "Filtered_IQtree"]
philr_part_weights = ["uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts"]
philr_ilr_weights = ["uniform","blw","blw.sqrt","mean.descendants"]
scoring = "accuracy"
col_names = ["metadata", "ilr_weight", "part_weight", "model", "split1", "split2", "split3", "split4", "split5", "split6", "split7", "split8", "split9", "split10"]
model_symbols = ["|", "<", "*", "x", "o", "^", ","]

for philr_group in philr_groups:
	main_output_label = f"{output_label}_{philr_group}"
	philr_dir = os.path.join(output_dir, "tables", philr_group)
	if not os.path.exists(philr_dir):
		print(f"""{philr_dir} does not exist.
			Use lib/cml_scripts/mlm_metrics/PhILR_random_trees_and_counts_tables.R 
			to create it.""")
		sys.exit()

	# --------------------------------------------------------------------------
	print("Importing metadata to working env.")
	# --------------------------------------------------------------------------

	meta_df = pd.read_csv(os.path.join(home_dir, project, options.meta_fn), \
		sep=options.delim, header=0, index_col=options.meta_index_col)

	if options.use_all_meta.lower() == "true":
		metad_cols = range(len(meta_df.columns))
	else:
		metad_cols = options.meta_col

	# --------------------------------------------------------------------------
	print(f"Establishing other constants {philr_group}.")
	# --------------------------------------------------------------------------

	result_fpath = os.path.join(output_dir, "tables", f"{main_output_label}_{project}_data.csv")
	pdf_fpath = os.path.join(output_dir, "graphics", f"bp_{main_output_label}_{project}.pdf")
	sum_pdf_fpath = os.path.join(output_dir, "graphics", f"sum_{main_output_label}_{project}.pdf")
	algo_table_fpath = os.path.join(output_dir, "tables", f"algo_{main_output_label}_{project}.csv")
	weight_fpath = os.path.join(output_dir, "tables", f"weight_{main_output_label}_{project}.csv")
	partw_fpath = os.path.join(output_dir, "graphics", f"pw_{main_output_label}_{project}.png")
	ilrw_fpath = os.path.join(output_dir, "graphics", f"iw_{main_output_label}_{project}.png")
	comb_weights_fpath = os.path.join(output_dir, "graphics", f"comb_weights_{main_output_label}_{project}.png")

	# prepare configuration for cross validation test harness
	seed = 7
	# --------------------------------------------------------------------------
	print(f"Preparing models for {philr_group}.")
	# --------------------------------------------------------------------------
	models = []
	models.append(('LR', LogisticRegression(max_iter=1000)))
	models.append(('LDA', LinearDiscriminantAnalysis()))
	models.append(('KNN', KNeighborsClassifier()))
	models.append(('DTREE', DecisionTreeClassifier()))
	models.append(('RF', RandomForestClassifier()))
	models.append(('GausNB', GaussianNB()))
	models.append(('SVM', SVC()))

	# --------------------------------------------------------------------------
	print(f"Building {scoring} scores for each weighing scheme.\
		Results found at {result_fpath} for {philr_group}.")
	# --------------------------------------------------------------------------
	with open(result_fpath, "w+") as fl:
		fl.write(",".join(col_names))
		fl.write("\n")
		for pw in philr_part_weights:
			for iw in philr_ilr_weights:
				table_fn = f"{philr_group}_{iw}_{pw}.csv"
				my_df = pd.read_csv(os.path.join(philr_dir, table_fn), sep=',', header=0, index_col=0)
				meta_df = meta_df.loc[list(my_df.index.values)]#drops rows from metadata that aren't in my_df
				assert list(my_df.index.values) == list(meta_df.index.values) #making sure I'm indexing correctly
				for meta_c in metad_cols:
					m_c = list(meta_df.columns)[meta_c]
					# print(m_c)
					spetz_var = meta_df[m_c]#metadata var to test
					# print(spetz_var.dtype)
					# assert is_string_dtype(spetz_var)
					# if spetz_var.dtype.name == "object":
					if is_string_dtype(spetz_var) == True and spetz_var.isnull().sum() < 5:
						print("evaluate each model in turn.")
						for name, model in models:
							kfold = model_selection.KFold(n_splits=10, random_state=seed, shuffle=True)
							cv_results = model_selection.cross_val_score(model, my_df, spetz_var, cv=kfold, scoring=scoring)
							# result_str = np.array2string(cv_results, separator=",",suffix="/n")
							result_str = ",".join(map(str, cv_results.tolist()))
							msg = f"{m_c},{iw},{pw},{name},{result_str}\n"
							fl.write(msg)
	# --------------------------------------------------------------------------
	print(f"Finished recording accuracy. Heading towards boxplot creation for {philr_group}.")
	# --------------------------------------------------------------------------
	#Setup for building boxplots
	result_df = pd.read_csv(result_fpath, sep=',', header=0)
	print(result_df.head())
	algos = list(set(result_df.loc[:,"model"]))
	algos.sort()
	print(algos)
	metadata_cats = list(set(result_df["metadata"]))
	print("metada_cats")
	print(metadata_cats)
	num_cols = 2
	num_rows = abs(-len(algos)//num_cols)
	print(result_df.head())

	f_header = ["feature", "f_mean","f_sd", "top_algo"]
	algo_table  = pd.DataFrame(columns = f_header)

	pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_fpath)
	# sub_plot_counter = 0
	for meta_c in metadata_cats:
		meta_result_df = pd.DataFrame(result_df[result_df["metadata"] == meta_c])
		flat_num_only = pd.DataFrame(meta_result_df.iloc[:,5:]).to_numpy().flatten()
		f_mean = np.nanmean(flat_num_only)
		f_sd = np.std(flat_num_only)
		print("f_mean", str(f_mean), type(f_mean))
		f_sd = meta_result_df.iloc[:,5:].values.std()
		f_max = max(meta_result_df)
		f_min = min(meta_result_df)
		f_q1 = np.quantile(flat_num_only, .25)
		f_q2 = np.quantile(flat_num_only, .50)
		f_q3 = np.quantile(flat_num_only, .75)
		#for boxplot
		fig = plt.figure(figsize=(11,11))
		fig.suptitle(f"{project} PhILR weighting ML accuracy {meta_c} for {philr_group}.", fontsize=20)
		plt.subplots_adjust(bottom=0.8)

		# for algo table
		algo_table_line = [str(meta_c), str(f_mean), str(f_sd)]
		algo_means = dict()
		summary_table = dict()

		if not pd.isna(f_mean):#don't proced if the data for the metadata is wonky
			fig_means = dict() #to hold data for summary boxplot figure
			for al in range(len(algos)):
				algo = algos[al]
				# print(algo)
				fig_df = pd.DataFrame(meta_result_df[meta_result_df["model"] == algo])
				# print(fig_df.head())
				plot_data = fig_df.iloc[:,5:].transpose()
				# print(f"df shape: {fig_df.shape[0]} {fig_df.shape[1]}")
				#side loop to get means of each
				fig_means[algo] = list(plot_data.mean(axis=0))
				#for boxplot
				ax = fig.add_subplot(num_rows,num_cols, al +1)
				ax.set_ylabel(scoring)
				ax.boxplot(plot_data)
				ax.title.set_text(f"{algo} by weighting scheme")
				ax.axhline(np.nanmean(plot_data), c="r", linestyle="dashed", label="Metadata mean")
				ax.axhline(f_mean, c="g", linestyle = ("-."), label = "Algorithm mean", marker = "_")
				ax.locator_params(axis='y', tight=True, nbins=4)
				new_labs = [f"{x}\n{y}" for x,y in zip(fig_df.loc[:,"ilr_weight"].values, fig_df.loc[:,"part_weight"].values)]
				# ax.set_xticklabels(fig_df.loc[:,"ilr_weight"].tolist(), rotation=90)
				ax.set_xticklabels(new_labs, rotation=90)
				ax.tick_params(axis='x', which='major', labelsize=6)
				ax.set_ylim([0.5, 1])#ax.set_ylim([f_mean - (1.5 * f_sd), 1])
				#for algo table
				algo_means[algo] = np.nanmean(plot_data) - f_mean

			ax = fig.add_subplot(num_rows,num_cols, al +2)
			ax.get_xaxis().set_visible(False)
			ax.get_yaxis().set_visible(False)
			ax.axhline(0, c="r", linestyle="dashed", label="Algorithm mean")
			ax.axhline(0, c="g", linestyle = ("-."), label = "Metadata mean", marker = "_")
			ax.legend(title="Legend", loc="center", framealpha=1, mode = "expand", markerscale=2)

			#for boxplot
			fig.tight_layout()
			pdf.savefig( fig )
			print(pd.DataFrame.from_dict(fig_means))
			print(list(fig_means.keys()))

			#for algo table
			algo_means["top_algo"] = max(algo_means, key=algo_means.get)
			algo_means["feature"] = meta_c
			algo_means["f_mean"] = f_mean
			algo_means["f_sd"] = f_sd
			algo_table = algo_table.append(algo_means, ignore_index=True)

	print(f"Saving pdf for {philr_group}")
	pdf.close()

	algo_table = algo_table.reindex(columns=f_header)
	algo_table = algo_table.round(decimals = 3)
	print(f"Saving algo table for or {philr_group}")
	algo_table.to_csv(algo_table_fpath, index = False)

	# --------------------------------------------------------------------------
	print(f"Building another summary table to show best weights for each feature and algo or {philr_group}.")
	# --------------------------------------------------------------------------
	w_header = ["feature", "algo", "algo_mean", "algo_sd", "top_pw", "top_ilr"]
	# For each feature, for each algo, top_pw, top_ilr
	weight_dict = dict(
					feature = [],
					algo_name = [],
					algo_comb_sd = [],
					algo_comb_mean = [],
					algo_top_mean =[],
					top_pw = [],
					top_ilr = [])

	summary_fig_data = { "feature" : [], "algo" : [], "weigh_scheme": [], "mean":[]}
	for meta_c in metadata_cats:
		f_mean = np.nanmean(result_df.iloc[:,4:])
		meta_result_df = pd.DataFrame(result_df[result_df["metadata"] == meta_c])
		f_mean = np.nanmean(meta_result_df.iloc[:,5:])
		f_sd = np.std(meta_result_df.iloc[:,5:], ddof=1)
		algo_means = dict()
		if not pd.isna(f_mean) and f_mean != 0:#don't proced if the data for the metadata is wonky
			for al in range(len(algos)):
				algo = algos[al]
				# print(algo)
				fig_df = pd.DataFrame(meta_result_df[meta_result_df["model"] == algo])
				plot_data = fig_df.iloc[:,5:]
				algo_mean = np.nanmean(plot_data)
				all_means = fig_df.mean(axis=1)# find all means of "plot data"
				print("all means")
				print(len(all_means))
				top_mean = max(all_means)#pick pw and il from highest means
				al_sd = plot_data.values.std()
				meta_result_df.iloc[:,5:].values.std()
				for ind, my_mean in enumerate(all_means):
					if my_mean == top_mean:
						top_index = ind
						# print(ind, my_mean)
				my_pw = list(fig_df["part_weight"])[top_index]
				my_iw = list(fig_df["ilr_weight"])[top_index]
				weight_dict["feature"].append(meta_c)
				weight_dict["algo_name"].append(algo)
				weight_dict["algo_comb_mean"].append(algo_mean)
				weight_dict["algo_comb_sd"].append(al_sd)
				weight_dict["algo_top_mean"].append(top_mean)
				weight_dict["top_pw"].append(my_pw)
				weight_dict["top_ilr"].append(my_iw)

	weight_table = pd.DataFrame(weight_dict)
	weight_table.to_csv(weight_fpath, index = False)
	# --------------------------------------------------------------------------
	print(f"Summary table to be found at: {weight_fpath} for {philr_group}.")
	# --------------------------------------------------------------------------

	# --------------------------------------------------------------------------
	print(f"Building summary scatterplot for {philr_group}.")
	# --------------------------------------------------------------------------
	num_rows = 5
	max_plots_per_page = 8
	pdf = matplotlib.backends.backend_pdf.PdfPages(sum_pdf_fpath)
	fig = plt.figure(figsize=(11,12))
	fig.suptitle(f"{project} Algorithm by PhILR weighting {scoring} mean or {philr_group}", fontsize = 20)
	# plt.subplots_adjust(bottom=0.8)
	sub_plot_counter = 0
	page_counter = 1
	for met in range(0,len(metadata_cats)):
		print(f"met val: {met}, numrows {num_rows}")
		meta_c = metadata_cats[met]
		meta_result_df = pd.DataFrame(result_df[result_df["metadata"] == meta_c])
		fig_means = dict() #to hold data for summary boxplot figure
		flat_num_only = pd.DataFrame(meta_result_df.iloc[:,5:]).to_numpy().flatten()
		f_mean = float(np.nanmean(flat_num_only))
		f_sd = float(np.std(flat_num_only))
		y_tick_interval = ((1.5 * f_sd) * 2) / 4
		if not pd.isna(f_mean) and f_mean > 0.5:#don't proced if the data for the metadata is wonky
			my_range = np.arange(round(f_mean - (1.5 * f_sd),1), round(f_mean + (1.5 * f_sd),1), round(y_tick_interval, 2))
			sub_plot_counter += 1
			if sub_plot_counter % max_plots_per_page == 0:
				pdf.savefig( fig )
				fig = plt.figure(figsize=(11,12))
				fig.suptitle(f"{project}, Algorithm by PhILR weighting {scoring} ")
			for al in range(len(algos)):
				algo = algos[al]
				# print(algo)
				fig_df = pd.DataFrame(meta_result_df[meta_result_df["model"] == algo])
				# print(fig_df.head())
				plot_data = fig_df.iloc[:,5:].transpose()
				print(f"df shape: {fig_df.shape[0]} {fig_df.shape[1]}")
				#side loop to get means of each
				fig_means[algo] = list(plot_data.mean(axis=0))
			my_means = pd.DataFrame.from_dict(fig_means)
			new_labs = [f"{x}\n{y}" for x,y in zip(fig_df.loc[:,"ilr_weight"].values, fig_df.loc[:,"part_weight"].values)]
			ax = fig.add_subplot(num_rows, num_cols, sub_plot_counter)
			ax.set_ylabel(scoring)
			# ax.scatter(range(0,my_means.shape[0]), my_means.iloc[:,0], label=my_means.columns[0])
			for co in range(0,my_means.shape[1]):
				# ax.scatter(range(0,my_means.shape[0]), my_means.iloc[:,co], label=my_means.columns[co], marker=model_symbols[co])
				ax.plot(range(0,my_means.shape[0]), 
								my_means.iloc[:,co], 
								marker=model_symbols[co], 
								label=my_means.columns[co], 
								linewidth=1, 
								markersize=2)
			# ax.title.set_text(f"Algorithm by weight, {scoring} mean, {meta_c}")
			ax.set_title(f"Algorithm by weight, {scoring} mean, {meta_c}", fontsize=10)
			ax.set_xticks(ticks=range(0,my_means.shape[0]), labels=fig_df.loc[:,"ilr_weight"].tolist(), rotation=90)
			ax.set_xticklabels(new_labs, rotation=90)
			ax.tick_params(axis='x', which='major', labelsize=6)
			# ax.set_yticks(ticks=my_range)
			ax.set_yticks(ticks=[0.55,0.65,0.75,0.85,0.95])

	# --------------------------------------------------------------------------
	print(f"Creating legend for {sum_pdf_fpath}.")
	# --------------------------------------------------------------------------
	ax = fig.add_subplot(num_rows, num_cols, sub_plot_counter + 1)
	for co in range(0,my_means.shape[1]):
		# ax.scatter(0,0, label=my_means.columns[co], marker=model_symbols[co],)
		ax.plot(0,0, label=my_means.columns[co], marker=model_symbols[co])
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)
	ax.legend(title="Algorithm", loc="center", framealpha=1, mode = "expand", markerscale=2)
	fig.tight_layout()
	pdf.savefig( fig )

	# --------------------------------------------------------------------------
	print(f"Saved summary scatterplot to {sum_pdf_fpath}.")
	# --------------------------------------------------------------------------
	pdf.close()

# --------------------------------------------------------------------------
print(f"{__file__} complete!")
# --------------------------------------------------------------------------
