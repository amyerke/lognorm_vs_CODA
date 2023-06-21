#!/usr/bin/env Rscript
# Author: Aaron Yerke
# Script for turning the alignment.rds from my dada2 script into a .aln file
# This is sort of an adhoc fix so that I don't have to re-run the dada2 
# script again

if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git',"lognorm_vs_CODA"), 
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default=NULL, 
                        help="project folder", metavar="character")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)


print("Load libraries.")
library("DECIPHER")

print("Establish directory layout.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

tree_dir <- file.path(output_dir, "trees")

if (! file.exists(tree_dir)){
    dir.create( tree_dir )
}

print("Creating table version of alignment.")
alignment <- readRDS(file.path(output_dir, "r_objects","ForwardReads_DADA2_alignment.rds"))
writeXStringSet(alignment, 
            file = file.path(tree_dir, "ForwardReads_DADA2_taxonomy.aln"))


