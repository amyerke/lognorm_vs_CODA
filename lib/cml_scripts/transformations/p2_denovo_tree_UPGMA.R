#!/usr/bin/env Rscript
# Author: Aaron Yerke (aaronyerke@gmail.com)
# This is a pipeline that was created to build denovo UPGMA trees.
# Many resources were used to create this pipeline:
#   Processing the sequences through DADA2, making trees with phangorn, and putting them into phyloseq objects:
#     This is by far the biggest source:
#     https://github.com/spholmes/F1000_workflow/blob/master/MicrobiomeWorkflow/MicrobiomeWorkflowII.Rmd
#     or (in article format)
#     https://f1000research.com/articles/5-1492/v2
rm(list = ls()) #clear workspace
##-cml argument processing------------------------------------------##
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git',"lognorm_vs_CODA"), 
              help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default="Jones", 
              help="project folder", metavar="character"),
  optparse::make_option(c("-m", "--metadata"), type="character", default="~/git/lognorm_vs_CODA/Jones/patient_metadata.tsv",
              help="metadata file path with filename", metavar="character"),
  optparse::make_option(c("-l", "--metadata_delim"), type="character", default="\t",
              help="metadata file deliminator", metavar="character"),
  optparse::make_option(c("-r", "--metadata_rowname"), type="character", default="Run",
              help="metadata file row to use for row names", metavar="character"),
  optparse::make_option(c("-s", "--outputfilesuffix"), type="character", default="",
                        help="output_file_suffix", metavar="character"),
  optparse::make_option(c("-f", "--filter_level"), type="numeric", default="0",
                        help="taxonimic level for making otu table 1-6", metavar="numeric"),
  optparse::make_option(c("-i", "--input_table"), type="character", default="ForwardReads_DADA2.rds",
                        help="input file name, must be in project/r_objects/"),
  optparse::make_option(c("-x", "--output_file_prefix"), type="character", default="",
                        help="output_file_prefix", metavar = "outputPrefix"),
  optparse::make_option(c("-a", "--alignment_file"), type="character", default="ForwardReads_DADA2_alignment.rds",
                        help="alignment file, should be an R object in output/r_objects/",
                        metavar = "algnmnt")
  ); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

####-load other dependencies------------------------------------------####
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager", type = "source", 
                   repos = "http://archive.linux.duke.edu/cran/")}
if (!requireNamespace("phangorn", quietly = TRUE)){
  install.packages("phangorn",type = "source", 
                   repos = "http://archive.linux.duke.edu/cran/")}
library("phangorn")
if (!requireNamespace("phyloseq", quietly = TRUE)){BiocManager::install("phyloseq")}
library("phyloseq")
if (!requireNamespace("DECIPHER", quietly = TRUE)){BiocManager::install("DECIPHER")}
library("DECIPHER")

print("external libraries loaded")

####-Establish directory layout---------------------------------------####
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')

print("Established directory layout")

####-Import R objects and data preprocessing--------------------------####
seqtab <- readRDS(file.path( output_dir, "r_objects", opt$input_table))
print(paste("Loaded seqtab."))
taxTab <- readRDS(file.path( output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))
alignment <- readRDS(file.path(output_dir, "r_objects",opt$alignment_file))
print("Imported R objects")

####-import tables----------------------------------------------------####
myMeta <- read.table(opt$metadata,
                    sep=opt$metadata_delim,
                    header=TRUE,
                    row.names = opt$metadata_rowname,
                    check.names = FALSE,
                    stringsAsFactors=FALSE)

print("Imported tables")

####-Build tree------------------------------------------------------####
phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
dm <- phangorn::dist.ml(phangAlign)#create distance matrix
treeNJ <- phangorn::upgma(dm) #make tree
fit <- phangorn::pml(treeNJ, data=phangAlign)#fit model
# fitGTR <- update(fit, k=4, inv=0.2)#fit model with updated parameters
# fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "stochastic", control = pml.control(trace = 0))
print("phangorn completed")

ps <- phyloseq::phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(myMeta),
               tax_table(taxTab),
               phy_tree(treeNJ))
# ps
print("Created ps")

#examine tree
library("ape")

pdf(file.path(output_dir, "graphics", paste0("upgma_denovo","_2",opt$outputfilesuffix, ".pdf")))

plot_tree(ps, "treeonly", nodeplotblank, ladderize="left")

dev.off()

saveRDS(ps, file.path(output_dir, "r_objects", paste0(opt$output_file_prefix, "denovo_tree_UPGMA_phyloseq_obj.rds")))

print("script complete")

