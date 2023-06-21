#!/usr/bin/env Rscript
# Author: Aaron Yerke
# This is a pipeline that was created to explore the philr transform.
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
  make_option(c("-r", "--r2_pattern"), type="character", default="_R2.fastq.gz", 
              help="pattern for reverse filenames for paired-end reads", 
              metavar="reverse read suffix"),
  make_option(c("-t", "--trunLen"), type="integer", default="240", 
              help="truncation length for paired-end reads", metavar="truncation length")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

##-Load Depencencies------------------------------------------------##
# ‘ape’, ‘dplyr’, ‘reshape2’, ‘plyr’
# .cran_packages <- c("ggplot2", "gridExtra")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DECIPHER", quietly = TRUE)) BiocManager::install("DECIPHER")
if (!requireNamespace("dada2", quietly = TRUE)) BiocManager::install("dada2",type = "source", checkBuilt = TRUE)
library("dada2")
library("DECIPHER")

print("Establish directory layout.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
f_path <- file.path(home_dir, project, "downloaded_seqs") # CHANGE ME to the directory containing the fastq files after unzipping.
# list.files(f_path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs <- sort(list.files(f_path, pattern=opt$r1_pattern, full.names = TRUE))
sampleNames <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filt_path <- file.path(f_path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sampleNames, "_R1_filt.fastq"))

#Only FAT the seqs that need it
completed_filtFs <- list.files(filt_path, full.names = TRUE)
# uncompleted_filtFs <- setdiff.Vector(filtFs, completed_filtFs)
# uncompleted_fnFs <- tail(fnFs, length(uncompleted_filtFs))


print("Attempting filter.")
out <- dada2::filterAndTrim(fnFs, filtFs,truncLen=opt$trunLen,
                            maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                            compress=FALSE, multithread=FALSE)
print("Completed filter and trim")

if (file.exists(file.path(output_dir,"r_objects", "dada2_dds.rds"))) {
  dds <- readRDS(file.path(output_dir,"r_objects", "dada2_dds.rds"))
  index <- length(dds) + 1
}else{
  index <- 1
  dds <- vector("list")
}

for (i in index:length(filtFs)){
  if (index <= length(filtFs)){
    print(paste0("On run", i, " of ", length(filtFs), "for dada2::dada"))
    f <- filtFs[i]
    errF <- dada2::learnErrors(f, multithread = FALSE)
    derepFs <- dada2::derepFastq(f, verbose=TRUE)
    dadaFs <- dada2::dada(derepFs, err=errF, multithread=FALSE)
    dds[[index]] <- dadaFs
    saveRDS(dds, file = file.path(output_dir,"r_objects", "dada2_dds.rds"))
    index <- index + 1
  }else{
    print("Already completed dada2.")
  }
}
names(dds) <- sampleNames

seqtab <- dada2::makeSequenceTable(dds)

#Removing chimeras
seqtab <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE)

# seqtabPath  <-file.path(output_dir,"dada2","ForwardReads_DADA2.rds")
# print(seqtabPath)
# fil <- tempfile(seqtabPath, fileext = ".rds")
setwd(file.path(output_dir, "tables"))
saveRDS(seqtab, file.path(output_dir, "r_objects","ForwardReads_DADA2.rds"))
write.table(seqtab,file.path(output_dir,"tables","ForwardReads_DADA2.txt"),sep="\t")
saveRDS(rowSums(seqtab), file.path(output_dir, "r_objects","raw_ASV_total_row_seqs.rds"))

fastaRef <- file.path(home_dir, "lib, "taxonomy", "rdp_train_set_16.fa.gz")
taxTab <- dada2::assignTaxonomy(seqtab, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))

seqs <- dada2::getSequences(seqtab)

names(seqs) <- seqs # This propagates to the tip labels of the tree

alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

saveRDS(taxTab, file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))
write.table(taxTab, 
            file = file.path(output_dir, "tables", "ForwardReads_DADA2_taxonomy.csv"),
            sep = ",")
saveRDS(alignment, file.path(output_dir, "r_objects","ForwardReads_DADA2_alignment.rds"))
write.table(alignment, 
            file = file.path(output_dir, "tables", "ForwardReads_DADA2_taxonomy.aln"),
            sep = ",")


print("Alignment completed")