#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
# Script for making figure showing selected transformations from unfiltered table
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
import argparse

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
	num_rand_trees = 3):
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
train_percent = float(options.training)
main_output_label = f"sklearn_random_forest_manual_{train_percent}train"
#info for random forest classification
result_fpath = os.path.join(output_dir, "tables", f"{main_output_label}.csv")
col_names = ["dataset", "metadata", "color"]
num_iterations = 20
col_names = col_names + [f"split{x}" for x in range(num_iterations)]
print(col_names)
#info for random forest feature importance
feature_pdf_fpath = os.path.join(output_dir, "graphics", f"seminar_feature_imp_{main_output_label}.pdf")
#info for boxplot
boxplot_pdf_fpath = os.path.join(output_dir, "graphics", f"seminar_bp_no_filt{main_output_label}.pdf")

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

# tables.append(("rarefied_DADA2", (os.path.join(output_dir, "tables", "1000_rrarefy_asv.csv"), ","), "purple"))
# tables.append(("lognorm_HashSeq", (os.path.join(output_dir,"tables", "lognorm_hashseq.csv"), ","), "lime"))
tables.append(("alr_DADA2", (os.path.join(output_dir, "tables", "alr_asv.csv"), ","), "white"))
# tables.append(("alr_HashSeq", (os.path.join(output_dir,"tables", "alr_hashseq.csv"), ","), "g"))
tables.append(("clr_DADA2", (os.path.join(output_dir, "tables", "clr_asv.csv"), ","), "white"))
# tables.append(("clr_HashSeq", (os.path.join(output_dir,"tables", "clr_hashseq.csv"), ","), "m"))
tables.append(("raw_DADA2",(os.path.join(output_dir, "tables", "ForwardReads_DADA2.txt"),"\t"),"white"))
# tables.append(("HashSeq", (os.path.join(output_dir,  "hashseq", "hashseq.csv"),","), "r"))
tables.append(("proportions_DADA2", (os.path.join(output_dir, "tables", "propotions_asv.csv"), ","), "lime"))
# tables.append(("propotions_DADA2", (os.path.join(output_dir, "tables", "propotions_asv.csv"), ","), "lime"))
tables.append(("Heilinger_DADA2", (os.path.join(output_dir, "tables", "heilinger_asv.csv"), ","), "lime"))
tables.append(("lognorm_DADA2", (os.path.join(output_dir, "tables", "lognorm_dada2.csv"), ","), "lime"))
# tables.append(("lognorm_Silva_DADA2", (os.path.join(output_dir, "tables", "lognorm_Silva.csv"), ","), "lime"))
tables.append(("Silva_DADA2", (os.path.join(output_dir,"tables", "Silva_DADA2", "Silva_DADA2.csv"), ","), "white"))
tables = add_PhILR_dfs_to_table(tables, os.path.join(output_dir, "tables", "Silva_DADA2"), "Silva_DADA2", color = "#050598")
tables = add_random_tree_PhILRs_to_table(tables, os.path.join(output_dir, "tables", "Silva_DADA2"), "Silva_DADA2", color = "#f7d8a0", num_rand_trees=3)
tables.append(("Filtered_Silva_DADA2", (os.path.join(output_dir,"tables", "Filtered_Silva_DADA2", "Filtered_Silva_DADA2.csv"), ","), "white"))
tables = add_PhILR_dfs_to_table(tables, os.path.join(output_dir, "tables", "Filtered_Silva_DADA2"), "Filtered_Silva_DADA2", color = "#050598")
tables = add_random_tree_PhILRs_to_table(tables, os.path.join(output_dir, "tables", "Filtered_Silva_DADA2"), "Filtered_Silva_DADA2", color = "#f7d8a0", num_rand_trees=3)
tables.append(("Filtered_UPGMA_DADA2", (os.path.join(output_dir,"tables", "Filtered_UPGMA_DADA2", "Filtered_UPGMA_DADA2.csv"), ","), "white"))
tables = add_PhILR_dfs_to_table(tables, os.path.join(output_dir, "tables", "Filtered_UPGMA_DADA2"), "Filtered_UPGMA_DADA2", color = "#050598")
tables = add_random_tree_PhILRs_to_table(tables, os.path.join(output_dir, "tables", "Filtered_UPGMA_DADA2"), "Filtered_UPGMA_DADA2", color = "#f7d8a0", num_rand_trees=3)
tables.append(("Filtered_IQtree", (os.path.join(output_dir,"tables", "Filtered_IQtree", "Filtered_IQtree.csv"), ","), "white"))
tables = add_PhILR_dfs_to_table(tables, os.path.join(output_dir, "tables", "Filtered_IQtree"), "Filtered_IQtree", color = "#050598")
tables = add_random_tree_PhILRs_to_table(tables, os.path.join(output_dir, "tables", "Filtered_IQtree"), "Filtered_IQtree", color = "#f7d8a0", num_rand_trees=3)

colors = list([sublist[-1] for sublist in tables])
print(colors)
labels = list([sublist[0] for sublist in tables])
print(labels)

# --------------------------------------------------------------------------
print(f"Building boxplot PDF.", flush = True)
# --------------------------------------------------------------------------
#Setup for building boxplots
result_df = pd.read_csv(result_fpath, sep=',', header=0, index_col=0)
metadata_cats = list(set(result_df["metadata"]))
num_cols = 2
num_rows = abs(-len(tables)//num_cols)
result_df = result_df.loc[labels]#remove rows that don't come from our custom table

colors = list([sublist[-1] for sublist in tables])
print(colors)
pdf = matplotlib.backends.backend_pdf.PdfPages(boxplot_pdf_fpath)
for meta_c in metadata_cats:
	fig = plt.figure(figsize=(11,11))
	fig.suptitle(f"{project} sklearn random forest manual {train_percent}training {scoring} {meta_c}")
	plt.subplots_adjust(bottom=0.8)
	meta_result_df = pd.DataFrame(result_df[result_df["metadata"] == meta_c])
	# flat_num_only = pd.DataFrame(meta_result_df.iloc[:,5:]).to_numpy().flatten()
	plot_data = meta_result_df.iloc[:,2:].transpose()
	ax = fig.add_subplot(1,1,1)
	bp = ax.boxplot(plot_data, patch_artist = True, labels=plot_data.columns)
	for patch, color in zip(bp['boxes'], colors):
		patch.set_facecolor(color)
	ax.axhline(np.nanmean(plot_data), c="brown", linestyle="dashed")
	ax.set_xticklabels(labels = plot_data.columns, rotation=90)
	ax.tick_params(axis='x', which='major', labelsize=10)
	#for boxplot
	fig.tight_layout()
	pdf.savefig( fig )

print(f"Saving pdf to f{boxplot_pdf_fpath}", flush = True)
pdf.close()

# --------------------------------------------------------------------------
print(f"{__file__} complete!")
# --------------------------------------------------------------------------s
