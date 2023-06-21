#!/usr/bin/env Rscript
# Author: Aaron Yerke (aaronyerke@gmail.com)
# This is a pipeline that was created to run the philr transform.
# Many resources were used to create this pipeline:
#   Processing the sequences through DADA2, making trees with phangorn, and putting them into phyloseq objects:
#     This is by far the biggest source:
#     https://github.com/spholmes/F1000_workflow/blob/master/MicrobiomeWorkflow/MicrobiomeWorkflowII.Rmd
#     or (in article format)
#     https://f1000research.com/articles/5-1492/v2
#     
#   DADA2:
#     https://benjjneb.github.io/dada2/tutorial.html
#   
#   Using philr:
#     https://bioconductor.org/packages/release/bioc/html/philr.html
#     Code example:
#       https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R
# rm(list = ls()) #clear workspace

if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git',"lognorm_vs_CODA"), 
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default=NULL, 
                        help="project folder", metavar="character"),
  optparse::make_option(c("-o", "--output_prefix"), type="character", default=NULL,
                        help="the file name output prefix to differentiate it from other otu tables", 
                        metavar="character"),
  optparse::make_option(c("-s", "--outputfilesuffix"), type="character", default="",
                        help="output_file_suffix", metavar="character"),
  optparse::make_option(c("-l", "--filter_level"), type="numeric", default=NULL,
                        help="min seq depth", metavar="numeric"),
  optparse::make_option(c("-f", "--phyloseq_obj"), type="character", 
                        default="denovo_tree_UPGMA_phyloseq_obj.rds",
                        help="phyloseq_obj file name without dir found in the r_objec", 
                        metavar="character")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

##-load other dependencies------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
  BiocManager::install("philr")
  BiocManager::install("ape")
}

library("phyloseq")
library("ape")
library(philr); packageVersion("philr")

##-Establish directory layout---------------------------------------##
home_dir <- opt$homedir
project <- opt$project
#home_dir = file.path('cloud','project')
output_dir = file.path(home_dir, project, 'output')
setwd(file.path(home_dir,project))

##-Import R objects and data----------------------------------------##
con <- gzfile(file.path( output_dir, "r_objects", opt$phyloseq_obj))
ps <- readRDS(con)
close(con)

seq_depth_df <- read.table(file.path(output_dir, "tables", "total_seq_depth.csv"),
                        sep = ",",
                        header = T)
print(length(row.names(seq_depth_df)))
print(dim(seq_depth_df))
seq_depth_df <- seq_depth_df[seq_depth_df[,1] > as.numeric(opt$filter_level),]
print(dim(seq_depth_df))
print(seq_depth_df[,1] > as.numeric(opt$filter_level))

# con <- gzfile(file.path( output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
# ps <- readRDS(con)
# close(con)

##-philr munging----------------------------------------------------##
ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
ps <- filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)

if (opt$filter_level > 0){
  ps <- phyloseq::prune_samples(row.names(seq_depth_df), ps)
}

ps <- transform_sample_counts(ps, function(x) x+1)

##-Test for reqs----------------------------------------------------##
print("rooted tree?")
is.rooted(phy_tree(ps)) # Is the tree Rooted?
print('All multichotomies resolved?')
is.binary.tree(phy_tree(ps)) # All multichotomies resolved?

## ---- message=FALSE, warning=FALSE-----------------------------------------
phy_tree(ps) <- ape::makeNodeLabel(phy_tree(ps), method="number", prefix='n')

name.balance(phy_tree(ps), tax_table(ps), 'n1')

##-philr transform--------------------------------------------------##
ps.philr <- philr(ps@otu_table, ps@phy_tree,
                  part.weights='enorm.x.gm.counts',
                  ilr.weights='blw.sqrt')

print("philR transform complete, next step is saving output")

##-Saving output----------------------------------------------------##
# commented out to activate ref_tree_philr
saveRDS(ps,
        file.path(output_dir, "r_objects", paste0("phyloseq_obj", opt$outputfilesuffix, ".rds")))
saveRDS(ps.philr,
        file.path(output_dir, "r_objects", paste0("philr_phyloseq_obj", opt$outputfilesuffix, ".rds")))
write.table(ps.philr,
            file.path(output_dir, "tables", paste0("philr", opt$outputfilesuffix, ".csv")),
            sep = ",")

# commented out to activate denovo tree
# saveRDS(ps,
#         file.path(output_dir, "r_objects", "ref_tree_ps_philr_transform.rds"))
# 
# write.table(ps.philr,
#             file.path(output_dir, "tables", "ref_tree_ps_philr_transform.csv"),
#             sep = ",")

# write.table(metadata,
#             file.path(output_dir, "tables", paste0("ps_sample_data", opt$outputfilesuffix, ".csv")),
#             sep = ",")
# 
# write.table(ps@otu_table,
#             file.path(output_dir, "tables", paste0(opt$output_prefix, "asv_table", opt$outputfilesuffix,".csv", )),
#             sep = ",")