#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
# This is a script for assessing the various compositional data 
# transformations against the random forest
# --------------------------------------------------------------------------
print("Loading external libraries.", flush = True)
# --------------------------------------------------------------------------
import imp
import os, sys
from statistics import mean
import math
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import accuracy_score, roc_auc_score, r2_score
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
parser.add_argument("-p", "--project", default="Jones",
                  help="project folder", metavar="project")

options, unknown = parser.parse_known_args()

# --------------------------------------------------------------------------
print("Defining functions", flush = True)
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
print("Establishing directory layout.", flush = True)
# --------------------------------------------------------------------------
home_dir = os.path.expanduser(options.homedir)
project = options.project
output_dir = os.path.join(home_dir, project, "output")
assert os.path.exists(output_dir)
input_dir = os.path.join(output_dir, "tables", "r_v_py_train_test_tables")
if not os.path.exists(input_dir):
	os.mkdir(input_dir)
os.chdir(input_dir)

# --------------------------------------------------------------------------
print("Establishing other constants", flush = True)
# --------------------------------------------------------------------------
scoring = "Accuracy"
main_output_label = f"sklearn_random_forest_from_file"
#info for random forest classification
result_fpath = os.path.join(output_dir, "tables", f"{main_output_label}.csv")
col_names = ["dataset", "metadata", "color"]
num_iterations = 20
col_names = col_names + [f"split{x}" for x in range(num_iterations)]
print(col_names)

# --------------------------------------------------------------------------
print("Importing data to working env.", flush = True)
# --------------------------------------------------------------------------


tables = []
tables.append(("raw_DADA2",(os.path.join(output_dir, "tables", "ForwardReads_DADA2.txt"),"\t"),"white"))

colors = list([sublist[-1] for sublist in tables])
print(colors)
labels = list([sublist[0] for sublist in tables])
print(labels)
# --------------------------------------------------------------------------
print(f"Reading train/testing files.", flush = True)
# --------------------------------------------------------------------------
my_files = os.listdir(input_dir)
print(f"number of files: {len(my_files)}")

metadata = []
iteration_min = float("nan")
iteration_max = float("nan")
for filen in my_files:
	if filen.endswith(".csv"):
		my_splits = filen.split("(_)")
		print(my_splits)
		meta = my_splits[0]
		iteration = int(my_splits[1])
		print( iteration)
		if not meta in metadata:
			metadata.append(meta)
		# if iteration_min is NaN, then make iteration, otherwise leave it iteration_min
		iteration_min = iteration if math.isnan(iteration_min) else iteration_min
		iteration_max = iteration if math.isnan(iteration_max) else iteration_max
		iteration_min = iteration if iteration < iteration_min else iteration_min
		iteration_max = iteration if iteration > iteration_max else iteration_max
print(metadata)
print(f"min{iteration_min} max{iteration_max}")

with open(result_fpath, "w+") as fl:
	fl.write(",".join(col_names))
	fl.write("\n")
	for m_c in metadata:
		my_accuracy = [0] * len(range(iteration_min,iteration_max + 1))
		for i in range(iteration_min,iteration_max + 1):
			print(f"{m_c}, {i}")
			pred_train = pd.read_csv(f"{m_c}(_){i}(_)pred(_)train.csv", header=0, index_col=0)
			pred_test = pd.read_csv(f"{m_c}(_){i}(_)pred(_)test.csv", header=0, index_col=0)
			resp_train = pd.read_csv(f"{m_c}(_){i}(_)resp(_)train.csv", header=0, index_col=0)
			resp_test = pd.read_csv(f"{m_c}(_){i}(_)resp(_)test.csv", header=0, index_col=0)
			new_m_c = m_c.replace("∕","/")
			if is_numeric_dtype(resp_test[new_m_c]) == True:
				print("going to RandomForestRegressor()")
				clf = RandomForestRegressor()
				clf.fit(pred_train, resp_train[new_m_c])
				resp_pred = clf.predict(pred_test)
				my_score = r2_score(resp_test[new_m_c], resp_pred, sample_weight=None)
			else:
				clf = RandomForestClassifier()
				clf.fit(pred_train, resp_train[new_m_c])
				resp_pred = clf.predict(pred_test)
				my_score = clf.score(pred_test, resp_test[new_m_c], sample_weight=None)
			my_accuracy[i] = my_score
			print(my_accuracy)
		final_acc = ",".join(map(str, my_accuracy))
		print(final_acc)
		msg = f"raw_DADA2,{m_c},black,{final_acc}\n"
		print(msg)
		fl.write(msg)
	print(f"{m_c} mean is: {mean(my_accuracy)}", flush = True)

# ----------------------------------------------------------------------------
print("Feeding tables to the random forest model.", flush = True)
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
print(f"{__file__} complete!")
# --------------------------------------------------------------------------s
