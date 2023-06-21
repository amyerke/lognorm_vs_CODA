#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
# This is a script for assessing the various compositional data 
# transformations against the random forest
# --------------------------------------------------------------------------
print("Loading external libraries.",flush = True)
# --------------------------------------------------------------------------
import os, sys
import time
from statistics import mean
from matplotlib import markers
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import matplotlib.colors as mcolors
from pandas.api.types import is_string_dtype
from requests import head
from sklearn import model_selection
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score, roc_auc_score
import argparse
import random

# --------------------------------------------------------------------------
print("Reading commmandline input with optparse.", flush = True)
# --------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="This script runs a random forest test on various datasets.")
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
parser.add_argument("-f", "--metada_fn", default="string", dest="meta_fn",
                  help="Name of file at the top of the project folder to use as metadata.", 
									metavar="meta_fn")
parser.add_argument("-l", "--delimiter", default="\t",
                  help="File delimiting symbol for metadata. Default is tab.",
									metavar="delim", dest="delim")
parser.add_argument("-i", "--meta_index_col", default=0,
                  help="Name of column to use as row name for metadata",
                  metavar="meta_index_col", dest="meta_index_col")
parser.add_argument("-t", "--training", default=0.9,
                  help="Percentating of table to use for training. The rest will be used for testing.",
                  metavar="training", dest="training")
parser.add_argument("-n", "--n_tree_step", default=100, type=int,
                  help="Step of increase in estimators.")
parser.add_argument("-s", "--n_tree_min", default=100, type=int,
                  help="How many n_tree to start.")
parser.add_argument("-r", "--n_tree_max", default=100, type=int,
                  help="How many n_tree to in last step.")
# parser.add_argument("-y", "--m_try", default=1/3, type=float,
#                   help="Max features for RFC, should be fraction.")

options, unknown = parser.parse_known_args()

# --------------------------------------------------------------------------
print("Defining functions", flush = True)
# --------------------------------------------------------------------------
def add_PhILR_dfs_to_table(lst, \
	root_folder, \
	base_fn, \
	# philr_part_weights = ["anorm","enorm"], \
	# philr_ilr_weights = ["blw.sqrt","mean.descendants"], \
	philr_part_weights = ["enorm"], \
	philr_ilr_weights = ["blw.sqrt"], \
	color = "w"):
	if not os.path.exists(root_folder):
		print(f"{root_folder} does not exist. Use PhILR_random_trees_and_counts_tables.R to create it.", flush = True)
		sys.exit()
	for pw in philr_part_weights:
		for iw in philr_ilr_weights:
			my_label = f"{base_fn}_{iw}_{pw}"
			table_fn = f"{my_label}.csv"
			# my_df = pd.read_csv(os.path.join(root_folder, table_fn), sep=',', header=0, index_col=0)
			lst.append((my_label, (os.path.join(root_folder, table_fn), ','), color))
	return lst

def add_random_tree_PhILRs_to_table(lst, \
	root_folder, \
	base_fn, \
	# philr_part_weights = ["anorm","enorm"], \
	# philr_ilr_weights = ["blw.sqrt","mean.descendants"], \
	philr_part_weights = ["enorm"], \
	philr_ilr_weights = ["blw.sqrt"], \
	color = "w", \
	num_rand_trees = 10):
	print(f"Adding random trees from {base_fn}.")
	if not os.path.exists(root_folder):
		print(f"{root_folder} does not exist. Use PhILR_random_trees_and_counts_tables.R to create it.", flush = True)
		sys.exit()
	for rand in range(1,num_rand_trees+1):
		for pw in philr_part_weights:
			for iw in philr_ilr_weights:
				my_label = f"Shuffle{rand}_PhILR_{base_fn}_{iw}_{pw}"
				table_fn = f"{my_label}.csv"
				my_df = pd.read_csv(os.path.join(root_folder, table_fn), sep=',', header=0, index_col=0)
				lst.append((my_label, (os.path.join(root_folder, table_fn), ','), color))
	return lst

def df_factory(my_path, my_sep):
	# Use for building df from "tables" list in random forest loop
	try:
		df = pd.read_csv(my_path, sep=my_sep, header=0, index_col=0)
		return df
	except Exception as e:
		print(f"An exception occurred during creation of dataframe from {my_path}", flush = True)
		print(e, flush=True)
		sys.exit(f"There was a problem loading {my_path}")

# --------------------------------------------------------------------------
print("Establishing directory layout.", flush = True)
# --------------------------------------------------------------------------
home_dir = os.path.expanduser(options.homedir)
project = options.project
output_dir = os.path.join(home_dir, project, "output")
assert os.path.exists(output_dir)

# --------------------------------------------------------------------------
print("Establishing other constants", flush = True)
# --------------------------------------------------------------------------
seed = 7
scoring = "Accuracy"
train_percent = options.training
main_output_label = f"sklrn_randmforst_manual_{train_percent}train_variable_ntree"
#info for random forest classification
result_fpath = os.path.join(output_dir, "tables", f"{main_output_label}.csv")
col_names = ["dataset", "metadata", "color", "n_trees"]
num_iterations = 20
col_names = col_names + [f"split{x}" for x in range(num_iterations)]
print(col_names)
#info for random forest feature importance
feature_pdf_fpath = os.path.join(output_dir, "graphics", f"feature_imp_{main_output_label}.pdf")
#info for boxplot
boxplot_pdf_fpath = os.path.join(output_dir, "graphics", f"bp_{main_output_label}.pdf")
time_dic = {}

# --------------------------------------------------------------------------
print("Importing data to working env.", flush = True)
# --------------------------------------------------------------------------
meta_df = pd.read_csv(os.path.expanduser(os.path.join(home_dir, project, str(options.meta_fn))), \
	sep=options.delim, header=0, index_col=options.meta_index_col)
print(meta_df.dtypes, flush=True)
meta_df = meta_df.select_dtypes(["object", "category", "string"])
metad_cols = range(len(meta_df.columns))
#Commented code for scrambling the rownames
# # meta_df.head
# # meta_df = meta_df.sample(frac=1)
# # meta_df.head

# ----------------------------------------------------------------------------
print("Setting up tables to feed the random forest model.", flush = True)
# --------------------------------------------------------------------------
tables = []
tables.append(("Raw_DADA2",(os.path.join(output_dir, "tables", "ForwardReads_DADA2.txt"),"\t"),"r"))
# tables.append(("HashSeq", (os.path.join(output_dir,  "hashseq", "hashseq.csv"),","), "r"))
tables.append(("lognorm_DADA2", (os.path.join(output_dir, "tables", "lognorm_dada2.csv"), ","), "y"))
# tables.append(("lognorm_HashSeq", (os.path.join(output_dir,"tables", "lognorm_hashseq.csv"), ","), "y"))
tables.append(("alr_DADA2", (os.path.join(output_dir, "tables", "alr_asv.csv"), ","), "g"))
# tables.append(("alr_HashSeq", (os.path.join(output_dir,"tables", "alr_hashseq.csv"), ","), "g"))
tables.append(("clr_DADA2", (os.path.join(output_dir, "tables", "clr_asv.csv"), ","), "g"))
# tables.append(("clr_HashSeq", (os.path.join(output_dir,"tables", "clr_hashseq.csv"), ","), "m"))
tables.append(("Silva_DADA2", (os.path.join(output_dir,"tables", "Silva_DADA2", "Silva_DADA2.csv"), ","), "#64baeb"))

# --------------------------------------------------------------------------
print(f"Running random forest model to find {scoring}.", flush = True)
# --------------------------------------------------------------------------
pdf = matplotlib.backends.backend_pdf.PdfPages(feature_pdf_fpath)
with open(result_fpath, "w+") as fl:
	fl.write(",".join(col_names))
	fl.write("\n")
	for nt in range(options.n_tree_min, options.n_tree_max, options.n_tree_step):
		for meta_c in metad_cols:
			m_c = list(meta_df.columns)[meta_c]
			# meta_df = meta_df.loc[list(my_table.index.values)
			for name, table_info, color in tables:
				my_table = df_factory(table_info[0], table_info[1])
				my_accuracy = [0] * num_iterations
				random.seed(10)
				for i in range(num_iterations):
					rand_int = random.randint(0, 1000)
					respns_var = meta_df.loc[list(my_table.index.values),m_c]#metadata var to test
					non_nan_vals = np.invert(pd.isna(respns_var.values))#removing nan values
					my_table = my_table[non_nan_vals]#removing nan values
					respns_var = respns_var[non_nan_vals]#removing nan values
					pred_train, pred_test, resp_train, resp_test = model_selection.train_test_split(my_table, respns_var, train_size=float(train_percent), random_state=rand_int, shuffle=True) 
					if is_string_dtype(respns_var) == True:
						clf = RandomForestClassifier(max_depth=2, random_state=0, n_estimators=nt)
						clf.fit(pred_train, resp_train)
						resp_pred = clf.predict(pred_test)
						# if len(set(respns_var)) <= 2:
						# 	resp_pred = list(x[0:len(resp_pred[1])-1] for x in resp_pred)
						# my_accuracy[i] = roc_auc_score(y_true = resp_test, y_score=resp_pred, multi_class="ovr", average="weighted")
						my_accuracy[i] = accuracy_score(resp_test, resp_pred)

				final_acc = ",".join(map(str, my_accuracy))
				print(final_acc)
				msg = f"{name},{m_c},{color},{nt}, {final_acc}\n"
				print(msg)
				fl.write(msg)
			print(f"{name}, {m_c} mean is: {mean(my_accuracy)}", flush = True)
			#Section for making figure regarding feature importance
			# feature_names = my_table.columns
			# start_time = time.time()
			# importances = clf.feature_importances_
			# std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)
			# elapsed_time = time.time() - start_time
			# print(f"Time to compute the importances: {elapsed_time:.3f} seconds for {name}{i}")
			# forest_importances = pd.Series(importances, index=feature_names)
			# # forest_importances = forest_importances.sort_values()¶[0:10]
			# fig, ax = plt.subplots()
			# forest_importances.plot.bar(yerr=std, ax=ax)
			# ax.set_title("Feature importances using MDI")
			# ax.set_ylabel("Mean decrease in impurity")
			# fig.tight_layout()
			# pdf.savefig( fig )
			#End Section for making figure regarding feature importance
print("Finished recording accuracy.", flush = True)
# print("Saving pdf", flush = True)
# pdf.close()

# --------------------------------------------------------------------------
print(f"Building scatterplot PDF.", flush = True)
# --------------------------------------------------------------------------
result_df = pd.read_csv(result_fpath, sep=',', header=0, index_col=0)
metadata_cats = list(set(result_df["metadata"]))
my_datasets = list(set(result_df["datasets"]))
splits = result_df.columns[result_df.columns.str.startswith('split')].tolist()
num_cols = 2
num_rows = abs(-len(tables)//num_cols)

pdf = matplotlib.backends.backend_pdf.PdfPages(boxplot_pdf_fpath)
for meta_c in metadata_cats:
	fig = plt.figure(figsize=(11,11))
	fig.suptitle(f"{project} sklearn randforest manual {train_percent}training {scoring} {meta_c} estimators vs accuracy")
	plt.subplots_adjust(bottom=0.8)
	ax = fig.add_subplot(1,1,1)
	meta_result_df = pd.DataFrame(result_df[result_df["metadata"] == meta_c])
	for ds in my_datasets:
		ds_result = meta_result_df[meta_result_df["dataset"] == ds]
		for rw in range(len(ds_result.index)):
			ax.scatter(ds_result.loc[rw,"ntrees"], mean(ds_result.loc[rw,splits]))
			# plt.annotate(f"{ds1}_{label}", (pval_log[i], ds1_lst[i]))

	#for boxplot
	fig.tight_layout()
	pdf.savefig( fig )

print("Saving pdf", flush = True)
pdf.close()
