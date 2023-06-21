#!/usr/bin/env Rscript
#script for creating fasta for blasting from DADA2 results
##-cml argument processing------------------------------------------##
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

option_list <- list(
  make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git',"lognorm_vs_CODA"), 
              help="dataset dir path", metavar="character"),
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder name in homedir", metavar="character")
  # make_option(c("-o", "--output"), type="character", default=NULL, 
  #             help="project folder name in homedir", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

#creating fasta for blasting from DADA2 results
home_dir <- opt$homedir
project <- opt$project
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, project, 'output')
o_file_name <- "dada2seqs.fasta"

con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2.rds"))
seqtab <- readRDS(con)

setwd(file.path( output_dir, "tree_process_blast"))
file.create(o_file_name)

for (i in 1:ncol(seqtab)){
  seq <- colnames(seqtab)[i]
  write(paste0(">", seq, "\n", seq, "\n"), 
    file=o_file_name, 
    append=TRUE)
}
getwd()
print("Fasta file created")