# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making lognorm and coda tables to run through my random forest
# Requires raw data2 output

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
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
library("rgr")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
if (!requireNamespace("vegan", quietly = TRUE)) BiocManager::install("vegan")
library("vegan")
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
  optparse::make_option(c("-m", "--metadata"), type="character", default=NULL,
                        help="metadata file path with filename", metavar="character"),
  optparse::make_option(c("-l", "--metadata_delim"), type="character", default="\t",
                        help="metadata file deliminator", metavar="character"),
  optparse::make_option(c("-r", "--metadata_rowname"), type="character", default=NULL,
                        help="metadata file row to use for row names", metavar="character")
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
df <- as.data.frame(rgr::alr(as.matrix(initial_table + 1), j = as.numeric(alr_col)))
saveRDS(df, file = file.path(output_dir,"r_objects", "alr_asv.rds"))
write.csv(df, file = file.path(output_dir,"tables", "alr_asv.csv"))
  
print("creating DADA2 CLR")
df <- as.data.frame(rgr::clr(as.matrix(initial_table + 1)))
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

# print("Creating hashseq lognorm, ALR and CLR.")
# hashseq <- data.frame(data.table::fread(file = file.path(output_dir,"hashseq", "hashseq.csv"),
#                                         header="auto", data.table=FALSE), row.names = 1)
# if (dir.exists(file.path(output_dir,"r_objects", "lognorm_hashseq.rds"))) {
#   ln_hs_tab <- readRDS(file.path(output_dir,"r_objects", "lognorm_HashSeq.rds"))
# }else{
# 	ln_hs_tab <- lognorm(hashseq)
#   saveRDS(ln_hs_tab, file = file.path(output_dir,"r_objects", "lognorm_hashseq.rds"))
#   write.csv(ln_hs_tab, file = file.path(output_dir,"tables", "lognorm_hashseq.csv"))
# }
# print("Making HashSeq clr.")
# if (dir.exists(file.path(output_dir,"r_objects", "r_objects", "clr_hashseq.rds"))) {
#   HashSeq_clr <- readRDS(file.path(output_dir,"r_objects", "clr_hashseq.rds"))
# }else{
#   HashSeq_clr <- as.data.frame(rgr::clr(as.matrix(hashseq + 1)))
#   saveRDS(ln_hs_tab, file = file.path(output_dir,"r_objects", "clr_hashseq.rds"))
#   write.csv(ln_hs_tab, file = file.path(output_dir,"tables", "clr_hashseq.csv"))
# }
# print("Making HashSeq alr.")
# my_zeros <- apply(hashseq, 2, function(x) {
#   return(sum(x == 0))
# })
# alr_col <- which(my_zeros == min(my_zeros))[1]
# if (dir.exists(file.path(output_dir,"r_objects", "r_objects", "alr_hashseq.rds"))) {
#   HashSeq_alr <- readRDS(file.path(output_dir,"r_objects", "alr_hashseq.rds"))
# }else{
#   HashSeq_alr <- as.data.frame(rgr::alr(as.matrix(hashseq + 1), j = as.numeric(alr_col)))
#   saveRDS(HashSeq_alr, file = file.path(output_dir,"r_objects", "alr_hashseq.rds"))
#   write.csv(HashSeq_alr, file = file.path(output_dir,"tables", "alr_hashseq.csv"))
# }

print("Reached end of script.")
