#!/usr/bin/env python3

"""
Author: Aaron Yerke (aaronyerke@gmail.com)
For determining if PhILR weighting schemes improve any ML techniques. 
"""

# --------------------------------------------------------------------------
print("Loading external libraries.")
# --------------------------------------------------------------------------
from cgi import print_directory
import os, sys
from posixpath import sep
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pycaret as pC
from pycaret.classification import *

print("Establishing directory layout.")
home_dir = os.path.join("~", "git", "lognorm_vs_CODA")
project = "Jones"
output_dir = os.path.join(home_dir, project, 'output')

print("Establishing other constants.")
main_output_label = "testing_pycaret"

print("Importing data to working env.")
my_df = pd.read_csv(os.path.join(output_dir, "tables", "philr_denovo_tree_UPGMA_1000_sd_filtered.csv"),\
   sep=',', header=0, index_col=0)
meta_df = pd.read_csv(os.path.join(home_dir, project, "patient_metadata.tsv"),\
   sep='\t', header=0, index_col=0)

print("Dropping any extra values from metadata.")
meta_df = meta_df.loc[list(my_df.index.values)]
list(my_df.index.values) == list(meta_df.index.values)

spetz_var = meta_df["type"]

new_df = my_df.join(spetz_var)

s = pC.classification.setup(new_df,"type")

print("made models")

best = pC.compare_models(include = ["lr", "svm", "rt", "knn"])
print("compared models")
print(best)
best_df = pd.DataFrame(best)
best_df.dataFrame.to_csv(os.path.join(output_dir, "tables", main_output_label + ".csv"))

print("Exiting python script")
