#!/usr/bin/env Rscript
# This script generates images that can be used to finetune the length of 
# It expects the dir layout to be as though you ran "create new subproject_dirs.R"
# the truncLen parameter of the filterAndTrim DADA2 function.

# rm(list = ls()) #clear workspace

##-cml argument processing------------------------------------------##
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

option_list <- list(
  make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git',"lognorm_vs_CODA"), 
              help="dataset dir path", metavar="home dir"),
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder", metavar="project"),
  make_option(c("-f", "--r1_pattern"), type="character", default="_R1.fastq.gz", 
              help="pattern for forward filenames for paired-end reads", 
              metavar="forward read suffix"),
  make_option(c("-r", "--r2_pattern"), type="character", 
              default="_R2.fastq.gz", 
              help="pattern for reverse filenames for paired-end reads", 
              metavar="reverse read suffix")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

##-Load Depencencies------------------------------------------------##
library(dada2); packageVersion("dada2")

##-Establish directory layout---------------------------------------##
home_dir <- opt$homedir
project <- opt$project
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
f_path <- file.path(home_dir, project, "downloaded_seqs") # CHANGE ME to the directory containing the fastq files after unzipping.list.files(f_path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(f_path, pattern=opt$r1_pattern, full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualF <- dada2::plotQualityProfile(fnFs[1:2])
message("created forward plotQualityProfile")

dir.create(file.path(output_dir, "graphics"), showWarnings = FALSE)
setwd(file.path(output_dir, "graphics"))

png(filename="plotQualF.png")
plot(plotQualF)
dev.off()
message("saved forward plotQualityProfile")

tryCatch({
  fnRs <- sort(list.files(f_path, pattern=opt$r2_pattern, full.names = TRUE))
  plotQualR <- dada2::plotQualityProfile(fnRs[1:2])
  message("Created reverse plotQualityProfile")
  png(filename="plotQualR.png")
  plot(plotQualR)
  dev.off()
  message("Saved reverse plotQualityProfile")
}, warning = function(w) {
  print('Opps, a warning is thrown')
  message(w)
}, error = function(e) {
  print('Opps, an error is thrown')
  message(cond)
})

message("Reached end of script")
