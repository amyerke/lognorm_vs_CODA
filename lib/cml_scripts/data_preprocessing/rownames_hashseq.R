# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for running mundging hashseq output https://github.com/FarnazFouladi/HashSeq
print("Clearing workspace.")
rm(list = ls()) #clear workspace

print("Install dependencies.")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library("optparse")

print("Loading cml options")
option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git',"lognorm_vs_CODA"), 
              help="dataset dir path", metavar="home dir"),
  optparse::make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder", metavar="project")); 

opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

print("Establish directory layout.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')

print("Reading data.")
hashseq <- data.frame(data.table::fread(file = file.path(output_dir,"hashseq", "SvTable.txt"),
                                        header=TRUE, data.table=FALSE), row.names = 1)
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))
# hashseq <- hashseq[, colSums(hashseq != 0) > 0.01*nrow(hashseq)]#remove columns that don't have at least 10%
# print(paste("HashSeq has", ncol(hashseq), "columns after column reduction."))
print("Removing '_' and anything after it.")
row.names(hashseq) <- sapply(as.character(row.names(hashseq)), function(x) {
  strsplit(x,"_")[[1]][1]
  })
print("Fixed formating of row names.")
print("First 10 row names.")
print(paste(row.names(hashseq)[1:10], collapse = " "))
print("Last 10 row names.")
print(paste(tail(row.names(hashseq), 10), collapse = " "))
print(paste("number of elements:", nrow(hashseq), collapse = " "))

print("Attempting to reorder rownames based on ASV table from dada2\n")
hashseq <- hashseq[ order(row.names(asv_table)),]
print("First 10 row names hashseq:")
print(paste(row.names(hashseq)[1:10], collapse = " "))
print("First 10 row names DADA2:")
print(paste(row.names(hashseq)[1:10], collapse = " "))
print("Last 10 row names hashseq:")
print(paste(tail(row.names(hashseq), 10), collapse = " "))
print("Last 10 row names hashseq:")
print(paste(tail(row.names(hashseq), 10), collapse = " "))
print(paste("Hashseq and DADA2 are equal?",
 identical(row.names(hashseq), row.names(asv_table)), "."))
print(all(row.names(hashseq) == row.names(asv_table)))
print(all.equal(row.names(hashseq), row.names(asv_table), check.attributes=FALSE))
print(paste("The different elements are:", collapse = " "))
print(paste(setdiff(row.names(hashseq), row.names(asv_table))))

print("Writing munged hashseq table.\n")

write.table(hashseq, file = file.path(output_dir,"hashseq", "hashseq.csv"),
            sep = ",")

print("R script completed.")



