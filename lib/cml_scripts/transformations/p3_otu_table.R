#!/usr/bin/env Rscript
# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making ratio table of OTUS
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

# rm(list = ls()) #clear workspace

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
  optparse::make_option(c("-t", "--taxonomic_level"), type="numeric", default="6",
                        help="taxonimic level for making otu table 1-6", metavar="numeric"),
  optparse::make_option(c("-o", "--output_prefix"), type="character", default="genus_",
                        help="the file name output prefix to differentiate it from other otu tables", metavar="character")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

##-load other dependencies------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("phyloseq")
}

##-Establish directory layout---------------------------------------##
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')

##-Source my own functions------------------------------------------##
source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Import R objects and data preprocessing--------------------------##
con <- gzfile(file.path( output_dir, "r_objects","ForwardReads_DADA2.rds"))
sv_table <- data.frame(readRDS(con))
close(con)

con <- gzfile(file.path( output_dir, "r_objects","ForwardReads_DADA2_taxonomy.rds"))
tax_table <- data.frame(readRDS(con))
close(con)

print(paste("dim sv_table:", dim(sv_table)))
print(paste("dim tax_table:", dim(tax_table)))

##-Create "OTU" table-----------------------------------------------##
otu_labels <- vector(length = ncol(sv_table), mode = "character")
#assign otus to the asvs
for (i in 1:ncol(sv_table)){
  asv <- names(sv_table)[i]
  otu_labels[i] <- tolower(as.character(tax_table[asv, as.numeric(opt$taxonomic_level)]))#this line sets the taxonomic level
}
uniq_otus <- unique(otu_labels)
uniq_otus <- uniq_otus[!is.na(uniq_otus)]

# set up main data.frame (read.table allows for prenaming columns)
# otu_tab <- read.table(text = "",
#                           colClasses = rep(x = "numeric", length(uniq_otus)),
#                           col.names = uniq_otus)
# otu_tab <- list()
otu_tab <- data.frame(row.names = row.names(sv_table))

for(otu1 in 1:length(uniq_otus)){
  my_otu_lab <- uniq_otus[otu1]
  my_otus <- otu_labels == my_otu_lab
  my_otus <- replace(my_otus, is.na(my_otus), F)
  #need to add vectors so that they are the same length as nrow(ps.philr)
  my_otus <- as.data.frame(sv_table[,my_otus])
  my_otus <- unlist(rowSums(my_otus, na.rm=T))
  otu_tab[,otu1] <- my_otus
}
colnames(otu_tab) <- uniq_otus

print("Saving OTU tables")

write.table(otu_tab, 
            file = file.path(output_dir, "tables", paste0(opt$output_prefix, "otu_table.csv")),
            sep = ",")

saveRDS(otu_tab, file.path(output_dir, "r_objects", paste0(opt$output_prefix, "otu_table.rds")))

##-Create "OTU" table-----------------------------------------------##
ln_otu_table <- lognorm(otu_tab)

print("Saving lognorm OTU tables")

write.table(otu_tab, 
            file = file.path(output_dir, "tables", paste0("lognorm_", opt$output_prefix, "otu_table.csv")),
            sep = ",")

saveRDS(otu_tab, file.path(output_dir, "r_objects", paste0("lognorm_", opt$output_prefix, "otu_table.rds")))

