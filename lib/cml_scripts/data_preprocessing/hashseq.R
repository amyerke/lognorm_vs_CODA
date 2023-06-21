# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for running hashseq https://github.com/FarnazFouladi/HashSeq

rm(list = ls()) #clear workspace

options(java.parameters="-Xmx500g")

print("Installing dependencies.")
# install.packages("devtools")
if (!requireNamespace("HashSeq", quietly = TRUE)){
  devtools::install_github("FarnazFouladi/HashSeq")
}
library("HashSeq")
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

print("Loading cml options")
option_list <- list(
  optparse::make_option(c("-i", "--input_dir"), type="character", 
              default=file.path('~','git',"lognorm_vs_CODA"), 
              help="dataset dir path", metavar="home dir"),
  optparse::make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="output_dir", metavar="project"),
  optparse::make_option(c("-t", "--threshold"), type="numeric", default=0,
              help="hashseq threshold")); 

opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

hash_in <- path.expand(opt$input_dir)
print(paste(hash_in, "\nexists:", file.exists(hash_in)))
hash_out <- path.expand(opt$output_dir)
print(paste(hash_out, "\nexists:", file.exists(hash_in)))
file.exists(hash_out)

print("Running hashseq")

HashSeq::inferTrueSequences(inputDir = hash_in, 
                            outputDir = hash_out,
                            abundanceThreshold = as.numeric(opt$threshold))

print("R Script complete.")

