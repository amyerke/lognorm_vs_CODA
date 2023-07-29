# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making building an alignment from an asv counts table

rm(list = ls()) #clear workspace

####-Functions--------------------------------------------------------####
print("Loading custom functions.")

####-Load dependencies------------------------------------------------####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
if (!requireNamespace("vegan", quietly = TRUE)) BiocManager::install("vegan")
library("vegan")
if (!requireNamespace("DECIPHER", quietly = TRUE)) BiocManager::install("DECIPHER")
library("DECIPHER")
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")}
library("optparse")
print("finished loading libraries")

####-Read command line arguements-------------------------------------####
option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git',"lognorm_vs_CODA"), 
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default="Vangay", 
                        help="project folder", metavar="character"),
  optparse::make_option(c("-i", "--input_table"), type="character", 
                        default="ForwardReads_DADA2_taxonomy.rds",
                        help="should be an r object, in project r_objects file"),
  optparse::make_option(c("-o", "--output_file"), type="character", 
                        default="ForwardReads_DADA2_taxonomy",
                        help="output file name - will be in output/tables and output/r_objects, don't add extension")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

####-Establish directory layout---------------------------------------####
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')

#asv tables are in each project in "project"/output/tables/ForwardReads_DADA2.txt
asv_table <- file.path("output", "r_objects", opt$input_table)
print("Finished establishing directory layout")

alignment <- DECIPHER::AlignSeqs(DNAStringSet(names(asv_table)), anchor=NA,verbose=FALSE)

saveRDS(alignment, file.path(output_dir, "r_objects", paste0(opt$output_file, ".rds")))
write.table(alignment, 
            file = file.path(output_dir, "tables", paste0(opt$output_file, ".aln")),
            sep = ",")
