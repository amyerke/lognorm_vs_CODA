# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for determining median read depth

rm(list = ls()) #clear workspace

####-Functions--------------------------------------------------------####
print("Loading custom functions.")

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
for (project in projects) {
  output_dir <- file.path(home_dir, project, 'output')
  initial_table <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds"))
  read_depth <- rowSums(initial_table)
  print(project)
  print(median(read_depth))
}
