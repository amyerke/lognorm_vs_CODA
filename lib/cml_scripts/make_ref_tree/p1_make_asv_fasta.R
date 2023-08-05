#!/usr/bin/env Rscript
#script for creating fasta for blasting from DADA2 results

####-read in external libraries---------------------------------------####
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

####-cml argument processing------------------------------------------####
option_list <- list(
  make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git',"lognorm_vs_CODA"), 
              help="dataset dir path", metavar="character"),
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder name in homedir", metavar="character"),
  make_option(c("-i", "--input_rds"), type="character", default="ForwardReads_DADA2.rds", 
              help="input rds in r_objects folder", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="dada2seqs.fasta", 
              help="name of output file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

####-Read in constants------------------------------------------------####
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
o_file_name <- opt$output
seqtab <- readRDS(file.path( output_dir, "r_objects", opt$input_rds))
print("Varibles are imported")

####-Creating fasta---------------------------------------------------####
setwd(file.path( output_dir, "tree_process_blast"))
file.create(o_file_name)

for (i in 1:ncol(seqtab)){
  seq <- colnames(seqtab)[i]
  write(paste0(">", seq, "\n", seq, "\n"),
    file=o_file_name,
    append=TRUE)
}
print("Fasta file created")