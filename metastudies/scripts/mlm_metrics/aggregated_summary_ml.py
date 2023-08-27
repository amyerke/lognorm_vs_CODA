#!/usr/bin/env python3
print(f"""
Starting: {__file__}
Author: Aaron Yerke (aaronyerke@gmail.com)
For determining if PhILR weighting schemes improve any ML techniques.
This was useful: https://machinelearningmastery.com/compare-machine-learning-algorithms-python-scikit-learn/
Compare Algorithms
This script creates aggregated summary mlm metrics pdf that combines all projects.
""")

# --------------------------------------------------------------------------
print("Loading external libraries.")
# --------------------------------------------------------------------------

from cProfile import label
import os, sys
from statistics import mean
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from pandas.api.types import is_string_dtype
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
                  help="path to git repository", dest="homedir", metavar="homedir")
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
projects = ["Jones", "Zeller","Vangay", "Noguera-Julian"]
proj_metaname = ["patient_metadata.tsv", "patient_metadata.csv", "patient_metadata.tsv", "patient_metadata.tsv"]
proj_meta_delim = ["\t",",","\t","\t"]
proj_meta_colname = ["Run", "Run", "run_accession", "Run"]
final_output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(final_output_dir)

# --------------------------------------------------------------------------
print("Establishing other constants.")
# --------------------------------------------------------------------------

output_label = "sklearn_ml_acc"
philr_groups = ["Silva_DADA2", "Filtered_Silva_DADA2", "Filtered_UPGMA_DADA2", "Filtered_IQtree",
		"prev_filt90_Silva_DADA2", "prev_filt90_UPGMA_DADA2", "prev_filt90_IQtree"]
philr_part_weights = ["uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts"]
philr_ilr_weights = ["uniform","blw","blw.sqrt","mean.descendants"]
scoring = "accuracy"
model_symbols = ["|", "<", "*", "x", "o", "^", ","]

"""
To make plot of aggregated summary plot
Combine result_df results for all philr combinations, ie: uniform vs uniform for each seperate MLA
The results look something like this:

metadata ilr_weight part_weight  model    split1    split2    split3    split4    split5    split6    split7    split8    split9   split10
ETHNICITY    uniform     uniform     LR  0.690909  0.781818  0.654545  0.818182  0.685185  0.666667  0.703704  0.722222  0.796296  0.759259
ETHNICITY    uniform     uniform    LDA  0.763636  0.745455  0.636364  0.818182  0.740741  0.685185  0.777778  0.759259  0.814815  0.814815
"""
# --------------------------------------------------------------------------
print(f"Aggregating data for aggregated summary scatterplot.")
# --------------------------------------------------------------------------
# Nested dict for storing aggragated variables for final plot
# keys will be {model}: {ilr_weight}_{part weight} : []
# and added in summary aggregate plot section
weights_models = {}

for proj in range(len(projects)):
	project = projects[proj]
	output_dir = os.path.join(home_dir, project, "output")
	for philr_group in philr_groups:
		main_output_label = f"{output_label}_{philr_group}"

		philr_dir = os.path.join(output_dir, "tables", philr_group)
		if not os.path.exists(philr_dir):
			print(f"""{philr_dir} does not exist.
				Use lib/cml_scripts/mlm_metrics/PhILR_random_trees_and_counts_tables.R 
				to create it.""")
			sys.exit()
		result_fpath = os.path.join(output_dir, "tables", f"{main_output_label}_{project}_{philr_group}_data.csv")
		result_df = pd.read_csv(result_fpath, sep=',', header=0)
		for line in range(result_df.shape[0]):
			#Dict for storing aggragated variables for final plot
			# Nested dict for storing aggragated variables for final plot
			# keys will be {model}: [{ilr_weight}_{part weight}:[],]
			my_model = result_df.loc[line, "model"]
			my_part = result_df.loc[line, "part_weight"]
			my_ilr = result_df.loc[line, "ilr_weight"]
			my_data = list(result_df.iloc[line, 4:])
			key_label = f"{my_part}_{my_ilr}"
			if not my_model in weights_models.keys():
				weights_models[f"{my_model}"] = {}
			if key_label in weights_models[f"{my_model}"]:
				weights_models[f"{my_model}"][f"{key_label}"].extend(my_data)
			else:
				weights_models[my_model][f"{key_label}"] = my_data

agg_sum_pdf_fpath = os.path.join(home_dir, "metastudies", "output", f"agg_summary_mla_philr.pdf")
pdf = matplotlib.backends.backend_pdf.PdfPages(agg_sum_pdf_fpath)
fig = plt.figure(figsize=(7,6))
fig.suptitle(f"Aggregated MLA by PhILR weighting scores")
ax = fig.add_subplot()
plt.subplots_adjust(bottom=0.3)
ax.set_ylabel("Score")
for num, model_key in enumerate(weights_models.keys()):
	print(model_key)
	weight_keys = sorted(weights_models[model_key].keys())
	means = []
	stds = []
	labels = []
	for weight_key in weight_keys:
		my_data = weights_models[model_key][weight_key]
		means.append(np.nanmean(weights_models[model_key][weight_key]))
		stds.append(np.nanstd(weights_models[model_key][weight_key]))
		labels.append(weight_key.replace("_", "\n"))
	ax.plot(labels,
					means,
					label=model_key,
					marker=model_symbols[num], 
					linewidth=1, 
					markersize=2)
	# ax.set_title(f"", fontsize=10)
	ax.set_xticks(ticks=range(0,len(means)), labels=labels, rotation=90)
	ax.set_yticks(ticks=[0.55,0.65,0.75,0.85,0.95])
	ax.tick_params(axis='x', which='major', labelsize=8)
	ax.legend()

pdf.savefig( fig )
pdf.close()

# --------------------------------------------------------------------------
print(f"{__file__} complete!")
# --------------------------------------------------------------------------