
# --------------------------------------------------------------------------
print("Loading external libraries.")
# --------------------------------------------------------------------------
import os, sys
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
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
import argparse

# --------------------------------------------------------------------------
print("Establishing directory layout.")
# --------------------------------------------------------------------------
home_dir = os.path.expanduser(os.path.join("~", "git", "lognorm_vs_CODA"))
project = "Jones"
output_dir = os.path.join(home_dir, project, "output")
assert os.path.exists(output_dir)

# --------------------------------------------------------------------------
print("Establishing other constants")
# --------------------------------------------------------------------------
seed = 7
scoring = "accuracy"
main_output_label = "sklearn_random_forest_cross_val_10fold"
result_fpath = os.path.join(output_dir, "tables", f"{main_output_label}_{project}_data.csv")
num_iterations = 10
col_names = ["dataset", "metadata"]
col_names = col_names + [f"split{x}" for x in range(num_iterations)]
print(col_names)
pdf_fpath = os.path.join(output_dir, "graphics", f"bp_{main_output_label}_{project}.pdf")

# --------------------------------------------------------------------------
print("Importing data to working env.")
# --------------------------------------------------------------------------
meta_df = pd.read_csv(os.path.join(home_dir, project, "patient_metadata.tsv"), \
	sep="\t", header=0, index_col="Run")
metad_cols = range(len(meta_df.columns))
hashseq_df = pd.read_csv(os.path.join(output_dir, "hashseq", "hashseq.csv"), sep=",", header=0, index_col=0)
asv_table = pd.read_csv(os.path.join(output_dir, "tables", "ForwardReads_DADA2.txt"), sep="\t", header=0, index_col=0)
clr_table = pd.read_csv(os.path.join(output_dir, "tables", "clr_asv.csv"), sep=",", header=0, index_col=0)
alr_table = pd.read_csv(os.path.join(output_dir, "tables", "alr_asv.csv"), sep=",", header=0, index_col=0)
ln_table = pd.read_csv(os.path.join(output_dir, "tables", "lognorm_asv.csv"), sep=",", header=0, index_col=0)
ln_hs_tab = pd.read_csv(os.path.join(output_dir,"tables", "lognorm_hashseq.csv"), sep=",", header=0, index_col=0)
HashSeq_clr = pd.read_csv(os.path.join(output_dir,"tables", "clr_hashseq.csv"), sep=",", header=0, index_col=0)
HashSeq_alr = pd.read_csv(os.path.join(output_dir,"tables", "alr_hashseq.csv"), sep=",", header=0, index_col=0)

# meta_df = meta_df.loc[list(asv_table.index.values)]#drops rows from metadata that aren't in asv_table
# if all(meta_df.index == hashseq_df.index):
# 	print("dataframes are the same.")

# --------------------------------------------------------------------------
print("Setting up tables to feed the random forest model.")
# --------------------------------------------------------------------------
tables = []
tables.append(("HashSeq", hashseq_df))
tables.append(("Raw_DADA2", asv_table))
tables.append(("lognorm_DADA2", ln_table))
tables.append(("lognorm_HashSeq", ln_hs_tab))
tables.append(("alr_DADA2", alr_table))
tables.append(("alr_HashSeq", HashSeq_alr))
tables.append(("clr_DADA2", clr_table))
tables.append(("clr_HashSeq", HashSeq_clr))

philr_part_weights = ["anorm","enorm"]
philr_ilr_weights = ["blw.sqrt","mean.descendants"]
silva_philr_dir = os.path.join(output_dir, "tables", "silva_philr_weights")
if not os.path.exists(silva_philr_dir):
  print(f"{silva_philr_dir} does not exist. Use silva_philr_weights.R to create it.")
  sys.exit()
for pw in philr_part_weights:
	for iw in philr_ilr_weights:
		table_fn = f"ref_tree_cln_{iw}_{pw}.csv"
		my_df = pd.read_csv(os.path.join(silva_philr_dir, table_fn), sep=',', header=0, index_col=0)
		my_label = f"ref_tree_philr_{iw}_{pw}"
		tables.append((my_label, my_df))

# meta_df.head
# meta_df = meta_df.sample(frac=1)
# meta_df.head

# --------------------------------------------------------------------------
print(f"Running random forest model to find {scoring}.")
# --------------------------------------------------------------------------
with open(result_fpath, "w+") as fl:
	fl.write(",".join(col_names))
	fl.write("\n")
	for name, my_table in tables:
		for meta_c in metad_cols:
			m_c = list(meta_df.columns)[meta_c]
			respns_var = meta_df.loc[list(my_table.index.values),m_c]#metadata var to test
			if is_string_dtype(respns_var) == True and respns_var.isnull().sum() < 5:
				kfold = model_selection.KFold(n_splits=10, random_state=seed, shuffle=True)
				cv_results = model_selection.cross_val_score(RandomForestClassifier(), my_table, respns_var, cv=kfold, scoring=scoring)
				# result_str = np.array2string(cv_results, separator=",",suffix="/n")
				result_str = ",".join(map(str, cv_results.tolist()))
				msg = f"{name},{m_c},{result_str}\n"
				fl.write(msg)
				print(f"{name}, {m_c} mean is: {mean(cv_results)}")
print("Finished recording accuracy.")

# --------------------------------------------------------------------------
print(f"Building boxplot PDF.")
# --------------------------------------------------------------------------
#Setup for building boxplots
result_df = pd.read_csv(result_fpath, sep=',', header=0)
metadata_cats = list(set(result_df["metadata"]))
num_cols = 2
num_rows = abs(-len(result_df.dataset.unique())//num_cols)

pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_fpath)
for meta_c in metadata_cats:
	fig = plt.figure(figsize=(11,11))
	fig.suptitle(f"{project} random forest {scoring} {meta_c}")
	plt.subplots_adjust(bottom=0.8)
	meta_result_df = pd.DataFrame(result_df[result_df["metadata"] == meta_c])
	flat_num_only = pd.DataFrame(meta_result_df.iloc[:,5:]).to_numpy().flatten()
	f_mean = np.nanmean(flat_num_only)
	plot_data = meta_result_df.iloc[:,2:].transpose()
	ax = fig.add_subplot(1,1, 1)
	# print(head(plot_data))
	ax.boxplot(plot_data)
	# ax.axhline(np.nanmean(plot_data), c="r", linestyle="dashed")
	ax.axhline(f_mean, c="g", linestyle = ("-."))
	ax.set_xticklabels(meta_result_df["dataset"].tolist(), rotation=90)
	ax.tick_params(axis='x', which='major', labelsize=10)
	#for boxplot
	fig.tight_layout()
	pdf.savefig( fig )

print("Saving pdf")
pdf.close()
