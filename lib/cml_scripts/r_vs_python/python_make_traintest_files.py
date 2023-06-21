#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
# This is a script for assessing the various compositional data 
# transformations against the random forest
# --------------------------------------------------------------------------
print("Loading external libraries.", flush = True)
# --------------------------------------------------------------------------
import os, sys
from statistics import mean
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from requests import head
from sklearn import model_selection
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
parser.add_argument("-p", "--project", default="Jones",
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
parser.add_argument("-t", "--training", default=0.75,
                  help="Percentating of table to use for training. The rest will be used for testing.",
                  metavar="training", dest="training")

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
result_table_dir = os.path.join(output_dir, "tables", "r_v_py_train_test_tables")
if not os.path.exists(result_table_dir):
	os.mkdir(result_table_dir)

# --------------------------------------------------------------------------
print("Establishing other constants", flush = True)
# --------------------------------------------------------------------------
scoring = "Accuracy"
train_percent = float(options.training)
main_output_label = f"sklearn_random_forest_manual_{train_percent}train"
#info for random forest classification
result_fpath = os.path.join(output_dir, "tables", f"{main_output_label}.csv")
col_names = ["dataset", "metadata", "color"]
num_iterations = 20
col_names = col_names + [f"split{x}" for x in range(num_iterations)]
print(col_names)

# --------------------------------------------------------------------------
print("Importing data to working env.", flush = True)
# --------------------------------------------------------------------------
meta_df = pd.read_csv(os.path.expanduser(os.path.join(home_dir, project, str(options.meta_fn))), \
	sep=options.delim, header=0, index_col=options.meta_index_col)
print(meta_df.dtypes, flush=True)
metad_cols = range(len(meta_df.columns))
#Commented code for scrambling the rownames
# # meta_df.head
# # meta_df = meta_df.sample(frac=1)
# # meta_df.head

# ----------------------------------------------------------------------------
print("Setting up tables to feed the random forest model.", flush = True)
# --------------------------------------------------------------------------
tables = []
tables.append(("raw_DADA2",(os.path.join(output_dir, "tables", "ForwardReads_DADA2.txt"),"\t"),"white"))

colors = list([sublist[-1] for sublist in tables])
print(colors)
labels = list([sublist[0] for sublist in tables])
print(labels)
# --------------------------------------------------------------------------
print(f"Building train/testing files.", flush = True)
# --------------------------------------------------------------------------

file_names = []
for meta_c in metad_cols:
	m_c = list(meta_df.columns)[meta_c]
	# meta_df = meta_df.loc[list(my_table.index.values)
	for name, table_info, color in tables:
		my_table = df_factory(table_info[0], table_info[1])
		my_accuracy = [0] * num_iterations
		random.seed(98)
		for i in range(num_iterations):
			rand_int = random.randint(0, 1000)
			respns_var = meta_df.loc[list(my_table.index.values),m_c]#metadata var to test
			non_nan_vals = np.invert(pd.isna(respns_var.values))#removing nan values
			my_table = my_table[non_nan_vals]#removing nan values
			respns_var = respns_var[non_nan_vals]#removing nan values
			pred_train, pred_test, resp_train, resp_test = model_selection.train_test_split(my_table, respns_var, train_size=float(train_percent), random_state=rand_int, shuffle=True)
			#file name: metadatafeature_iteration_num_response/predictor_train/test
			my_filename = f"{m_c}(_){i}(_)pred(_)train.csv".replace("/","∕")
			file_names.append(my_filename)
			pred_train.to_csv(os.path.join( result_table_dir, my_filename))
			my_filename = f"{m_c}(_){i}(_)pred(_)test.csv".replace("/","∕")
			file_names.append(my_filename)
			pred_test.to_csv( os.path.join(result_table_dir, my_filename))
			my_filename = f"{m_c}(_){i}(_)resp(_)train.csv".replace("/","∕")
			file_names.append(my_filename)
			resp_train.to_csv( os.path.join(result_table_dir, my_filename))
			my_filename = f"{m_c}(_){i}(_)resp(_)test.csv".replace("/","∕")
			file_names.append(my_filename)
			resp_test.to_csv( os.path.join(result_table_dir, my_filename))
print(f"""
File names:
""")
print(file_names)
# --------------------------------------------------------------------------
print(f"{__file__} complete!")
# --------------------------------------------------------------------------s
