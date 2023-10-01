# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for determining what proportion of ASV's went into
# the Silva LTP

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
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
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
)

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

print("Establishing directory layout and other constants.")
home_dir <- opt$homedir
projects <- c("Jones", "Vangay", "Zeller", "Noguera-Julian")

source(file.path(home_dir, "lib", "table_manipulations.R"))

#### Importing and preprocessing tables, starting with DADA2 ####
print("Importing and preprocessing tables, starting with DADA2")

output_df <- data.frame("project" = character(length = 4),
                        "raw_counts" = numeric(length = 4),
                        "raw_counts_SilvaLTP" = numeric(length = 4),
                        "raw_proportion_mapped" = numeric(length = 4),
                        "prev_filt_counts" = numeric(length = 4),
                        "prev_filt_SilvaLTP" = numeric(length = 4),
                        "pf_proportion_mapped" = numeric(length = 4)
                        )

print(output_df)

for (pro in 1:length(projects)) {
	print(pro)
  project <- projects[pro]
	print(project)
  output_df$project[pro] <- project
  output_dir <- file.path(home_dir, project, "output")
  initial_table <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds"))
  output_df$raw_counts[pro] <- ncol(initial_table)
  silva <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
  output_df$raw_counts_SilvaLTP[pro] <- ncol(silva@otu_table)
  output_df$raw_proportion_mapped[pro] <- ncol(silva@otu_table) / ncol(initial_table)

  initial_table <- readRDS(file.path(output_dir, "r_objects", "filtered_90prcnt_dada2.rds"))
  output_df$prev_filt_counts[pro] <- ncol(initial_table)
  silva <- readRDS(file.path(output_dir, "r_objects", "ltp_90prcnt_filt.rds"))
  output_df$prev_filt_SilvaLTP[pro] <- ncol(silva@otu_table)
  output_df$pf_proportion_mapped[pro] <- ncol(silva@otu_table) / ncol(initial_table)
}
