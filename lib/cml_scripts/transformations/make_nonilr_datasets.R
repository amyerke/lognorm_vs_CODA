# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making noncoda and coda tables to run through my random forest
# Requires raw data2 output. It also makes a 90% filtered asv table and makes
# the same transformations using this table.
# I.e. #creating an asv dataset where asvs that are in less than 10% of the
# samples are removed

rm(list = ls()) #clear workspace

####-Functions--------------------------------------------------------####
print("Loading custom functions.")
raw_ps_to_clean_ps <- function(ps) {
  #requires ape, phyloseq and philr_tutorial_normalization 
  #performs philr tutorial normalization on phyloseq obj
  clean_otu = data.frame(ps@otu_table@.Data)
  clean_otu = philr_tutorial_normalization(clean_otu)
  ps_clean = phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                 phy_tree(ps@phy_tree),
                                 tax_table(ps@tax_table), 
                                 sample_data(ps@sam_data))
  return(ps_clean)
}

####-Load Depencencies------------------------------------------------####
options(repos = list(CRAN="http://cran.rstudio.com/"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("compositions", quietly = TRUE)) BiocManager::install("compositions")
library("compositions")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
if (!requireNamespace("vegan", quietly = TRUE)) BiocManager::install("vegan")
library("vegan")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
library("phyloseq")
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")}
library("optparse")
print("finished loading libraries")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git',"lognorm_vs_CODA"), 
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default="Vangay", 
                        help="project folder", metavar="character"),
  optparse::make_option(c("-i", "--initial_table"), type="character", 
                        default="/r_objects/ForwardReads_DADA2.rds",
                        help="initial table, relative to project path, with filename, 
                        should be an r object"),
  optparse::make_option(c("-r", "--output_file_prefix"), type="character", 
                        default="",
                        help="prefix to identify different types of output files")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

print("Establishing directory layout and other constants.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')

source(file.path(home_dir, "lib", "table_manipulations.R"))

#### Set up constants ####
random_seed <- 36

#### Importing and preprocessing tables, starting with DADA2 ####
print("Importing and preprocessing tables, starting with DADA2")

initial_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

print("creating DADA2 lognorm, ALR and CLR")
if (!dir.exists(file.path(output_dir,"r_objects", "lognorm_asv.rds"))) {
  df <- lognorm(initial_table)
  saveRDS(df, file = file.path(output_dir,"r_objects", "lognorm_asv.rds"))
  write.csv(df, file = file.path(output_dir,"tables", "lognorm_dada2.csv"))
}
# if (!dir.exists(file.path(output_dir,"r_objects", "silva_lognorm.rds"))) {
  silva_ps_robj <- readRDS(file.path(output_dir, "r_objects","ref_tree_phyloseq_obj.rds"))
  df <- na.omit(silva_ps_robj@otu_table) 
  df <- na.omit(lognorm(df))
  saveRDS(df, file = file.path(output_dir,"r_objects", "lognorm_Silva.rds"))
  write.csv(df, file = file.path(output_dir,"tables", "lognorm_Silva.csv"))
# }
my_zeros <- apply(initial_table, 2, function(x) {
  return(sum(x == 0))
})
alr_col <- which(my_zeros == min(my_zeros))[1]
# alr_col_num <- grep(alr_col, colnames(initial_table))
print("creating DADA2 ALR")
df <- as.data.frame(compositions::alr(as.matrix(initial_table + 1), ivar = as.numeric(alr_col)))
saveRDS(df, file = file.path(output_dir,"r_objects", "alr_asv.rds"))
write.csv(df, file = file.path(output_dir,"tables", "alr_asv.csv"))
  
print("creating DADA2 CLR")
df <- as.data.frame(compositions::clr(as.matrix(initial_table + 1)))
saveRDS(df, file = file.path(output_dir,"r_objects", "clr_asv.rds"))
write.csv(df, file = file.path(output_dir,"tables", "clr_asv.csv"))

print("creating DADA2 proportions")
df <- proportions_transform(initial_table)
saveRDS(df, file = file.path(output_dir,"r_objects", "propotions_asv.rds"))
write.csv(df, file = file.path(output_dir,"tables", "propotions_asv.csv"))

print("creating DADA2 heilinger")
df <- heilinger_transform(initial_table)
saveRDS(df, file = file.path(output_dir,"r_objects", "heilinger_asv.rds"))
write.csv(df, file = file.path(output_dir,"tables", "heilinger_asv.csv"))

print("Creating rarefied dataset at 1000 read depth")
rd <- 1000
low_abund_removed <- filt_seq_dpth(rd, initial_table)
df <- vegan::rrarefy(low_abund_removed, rd)
saveRDS(df, file = file.path(output_dir,"r_objects", paste0(rd,"_rrarefy_asv.rds")))
write.csv(df, file = file.path(output_dir,"tables", paste0(rd,"_rrarefy_asv.csv")))

#### Creating the filtered dataset and repeating all transformations ####
initial_table <- initial_table[sapply(initial_table, function(x) mean(x == 0) <= 0.9)]

print("creating filtered DADA2 lognorm, ALR and CLR")
if (!dir.exists(file.path(output_dir,"r_objects", "lognorm_asv.rds"))) {
  df <- lognorm(initial_table)
  saveRDS(df, file = file.path(output_dir,"r_objects", "filtered_90prcnt_lognorm_asv.rds"))
  write.csv(df, file = file.path(output_dir,"tables", "filtered_90prcnt_lognorm_dada2.csv"))
}
#TODO: Add filtered silva lognorm
# # if (!dir.exists(file.path(output_dir,"r_objects", "silva_lognorm.rds"))) {
# silva_ps_robj <- readRDS(file.path(output_dir, "r_objects","ref_tree_phyloseq_obj.rds"))
# df <- na.omit(silva_ps_robj@otu_table) 
# df <- na.omit(lognorm(df))
# saveRDS(df, file = file.path(output_dir,"r_objects", "lognorm_Silva.rds"))
# write.csv(df, file = file.path(output_dir,"tables", "lognorm_Silva.csv"))
# # }
my_zeros <- apply(initial_table, 2, function(x) {
  return(sum(x == 0))
})
alr_col <- which(my_zeros == min(my_zeros))[1]
# alr_col_num <- grep(alr_col, colnames(initial_table))
print("creating DADA2 ALR")
df <- as.data.frame(compositions::alr(as.matrix(initial_table + 1), ivar = as.numeric(alr_col)))
saveRDS(df, file = file.path(output_dir,"r_objects", "filtered_90prcnt_alr_asv.rds"))
write.csv(df, file = file.path(output_dir,"tables", "filtered_90prcnt_alr_asv.csv"))

print("creating DADA2 CLR")
df <- as.data.frame(compositions::clr(as.matrix(initial_table + 1)))
saveRDS(df, file = file.path(output_dir,"r_objects", "filtered_90prcnt_clr_asv.rds"))
write.csv(df, file = file.path(output_dir,"tables", "filtered_90prcnt_clr_asv.csv"))

print("creating DADA2 proportions")
df <- proportions_transform(initial_table)
saveRDS(df, file = file.path(output_dir,"r_objects", "filtered_90prcnt_propotions_asv.rds"))
write.csv(df, file = file.path(output_dir,"tables", "filtered_90prcnt_propotions_asv.csv"))

print("creating DADA2 heilinger")
df <- heilinger_transform(initial_table)
saveRDS(df, file = file.path(output_dir,"r_objects", "filtered_90prcnt_heilinger_asv.rds"))
write.csv(df, file = file.path(output_dir,"tables", "filtered_90prcnt_heilinger_asv.csv"))

print("Creating rarefied dataset at 1000 read depth")
filtered_prefix <- "filtered_90prcnt_"
rd <- 1000
low_abund_removed <- filt_seq_dpth(rd, initial_table)
df <- vegan::rrarefy(low_abund_removed, rd)
saveRDS(df, file = file.path(output_dir,"r_objects", paste0(filtered_prefix,rd,"_rrarefy_asv.rds")))
write.csv(df, file = file.path(output_dir,"tables", paste0(filtered_prefix,rd,"_rrarefy_asv.csv")))

print("Reached end of script.")

