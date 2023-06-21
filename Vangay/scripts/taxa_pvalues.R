# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for creating taxo anova pvalues.
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
makeTaxaTable <- function(no, tax, tax_lev_int){
  #no is otu table (dataframe)
  #tax is taxonomic table (dataframe)
  #tax_lev_int is the columns of the tax table usually 1-6
  
  
  otu_labels = vector(length = ncol(no), mode = "character")
  #assign otus to the asvs
  for (i in 1:ncol(no)){
    asv = names(no)[i]
    otu_labels[i] = tolower(as.character(tax[asv, tax_lev_int]))#this line sets the taxonomic level
  }
  uniq_otus = unique(otu_labels)
  
  otu_tab = data.frame(row.names = row.names(no))
  
  for(otu1 in 1:length(uniq_otus)){
    my_otu_lab = uniq_otus[otu1]
    my_otus = otu_labels == my_otu_lab
    my_otus = replace(my_otus, is.na(my_otus), F)
    #need to add vectors so that they are the same length as nrow(ps.philr)
    my_otus = as.data.frame(no[,my_otus])
    my_otus = unlist(rowSums(my_otus, na.rm=T))
    otu_tab[,otu1] =  my_otus
  }
  colnames(otu_tab) = uniq_otus
  return(otu_tab)
}


##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Vangay"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "r_libraries", "statistical_functions.R"))
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

##-Import tables and data preprocessing-----------------------------##
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
asv_tax <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds")))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(asv_table),]

##-Create pval tables-----------------------------------------------##
taxa_int <- c()
taxa_name <- c()
meta_name <- c()
meta_col <- c()
pval <- c()
taxa_lev <- c()
for( taxa_col in 1:ncol(asv_tax)){
  my_table <- makeTaxaTable(asv_table, asv_tax, taxa_col)
  for( tx in 1:ncol(my_table)){
    my_taxa <- names(my_table)[tx]
    for(meta in 1:ncol(metadata)){
      my_meta <- metadata[,meta]
      my_lm <- lm(unlist(asv_table[,tx]) ~ my_meta)
      my_pval <- anova(my_lm)$"Pr(>F)"[1]
      pval <- c(pval, my_pval)
      taxa_int <- c(taxa_int, taxa_col)
      taxa_lev <- c(taxa_lev, names(asv_tax)[taxa_col])
      taxa_name <- c(taxa_name, my_taxa)
      meta_name <- c(meta_name, names(metadata)[meta])
      meta_col <- c(meta_col, meta)
    }
  }
}  

adj_pval <- p.adjust(pval, method = "BH")

dFrame <- data.frame(taxa_lev, taxa_name, taxa_int, meta_name, meta_col, pval, adj_pval)
dFrame <- dFrame [order(dFrame$pval),]
dFrame$adj_pval <- p.adjust( dFrame$pval, method = "BH" )	
write.table(dFrame, file=file.path(output_dir, "tables", "Vangay_pValuesUnivariate_taxaVmetadata.tsv"), 
            sep="\t", 
            row.names=FALSE)



