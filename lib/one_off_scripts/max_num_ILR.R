# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script estimating max possible ILRs in dataset using the formula
# (2J−2)!/(2^(J−1) * (J−1)!) as found in "The isometric logratio 
# transformation in compositional data analysis: a practical evaluation
rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Jones"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Import tables and data preprocessing-----------------------------##
asv_table <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds"))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
asv_tax <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)


all_possble_ilrs <- function(num_balances) {
  # the formula is: (2J−2)!/(2^(J−1) * (J−1)!)
  return( (factorial(2*num_balances-2)/(2^(num_balances - 1) * factorial(num_balances*1) )))
}

num_taxa = c(5, 10, 20, 50, 53, 100)

for (ta in num_taxa){
  print(paste("number of ILRs for", ta, ":", all_possble_ilrs(ta)))
}

# [1] "number of ILRs for 5 : 21"
# [1] "number of ILRs for 10 : 3445942.5"
# [1] "number of ILRs for 20 : 4.10039726631895e+20"
# [1] "number of ILRs for 50 : 5.50584270656709e+74"
# [1] "number of ILRs for 53 : 5.34948196789134e+80"
# [1] "number of ILRs for 100 : Inf"
