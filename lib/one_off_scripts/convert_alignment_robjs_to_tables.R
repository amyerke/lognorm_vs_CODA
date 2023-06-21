#!/usr/bin/env Rscript
# Author: Aaron Yerke

print("Load libraries.")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DECIPHER", quietly = TRUE)) BiocManager::install("DECIPHER")
library("DECIPHER")

print("Establish directory layout.")
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Jones"
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

print("Creating table version of alignment.")
alignment <- readRDS(file.path(output_dir, "r_objects","ForwardReads_DADA2_alignment.rds"))
writeXStringSet(alignment, 
            file = file.path(output_dir, "trees", "ForwardReads_DADA2_taxonomy.aln"))

print("Reached end of script.")


