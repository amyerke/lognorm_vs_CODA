# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making table of metadata to read depth pvalues or correlations

rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
print("finished loading libraries")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
projects <- c("Jones", "Zeller","Vangay", "Noguera-Julian")
proj_metaname <- c("patient_metadata.tsv", "patient_metadata.csv", "patient_metadata.tsv", "patient_metadata.tsv")
proj_meta_delim <- c("\t",",","\t","\t")
proj_meta_colname <- c("Run", "Run", "run_accession", "Run")
final_output_dir <- file.path(home_dir, "metastudies", "output")
##-Functions--------------------------------------------------------##

#Correlations for continuous data and pvalues for factors
cor_df <- data.frame(project <- c(), feature <- c(), kendal_corr <- c())
# corr_proj <- c()
# corr_meta <- c()
# corr_value <- c()

pval_df <- data.frame(project <- c(), feature <- c(), pvalue <- c(), "pval < 0.05" <- c())
# pval_proj <- c()
# pval_meta <- c()
# pval_val <- c()
# pval_adj_pval <- c()

for (proj in 1:length(projects)) {
  project <- projects[proj]
  print(project)
  proj_out <- file.path(home_dir, project, "output")
  metadata <- data.frame(data.table::fread(file = file.path(home_dir, project, proj_metaname[proj]),
                                           header=TRUE, data.table=FALSE), row.names = proj_meta_colname[proj])
  my_table <- readRDS(file.path(proj_out, "r_objects", "ForwardReads_DADA2.rds"))
  read_depth <- base::rowSums(my_table)
  metadata <- metadata[row.names(my_table),]
  
  meta_class <- sapply(metadata,class)
  for (colmn in 1:ncol(metadata)){
    if (meta_class[colmn] == "integer" || meta_class[colmn] == "numeric"){
      my_cor <- base::round(stats::cor(read_depth, metadata[,colmn], method = "kendal"),4)
      # corr_proj <- c(corr_proj, project)
      # corr_meta <- c(colnames(metadata)[colmn], corr_meta)
      # corr_value <- c(corr_value, my_cor)
      new_row <- c(project, colnames(metadata)[colmn], my_cor)
      cor_df <- rbind(cor_df, new_row)
    }else{
      my_lm <- lm(read_depth ~ metadata[,colmn])
      my_pval <- anova(my_lm)$"Pr(>F)"[1]
      if (my_pval < 0.0001){
        my_pval <- "<0.0001"
      }else{
        my_pval <- base::round(my_pval, 4)
      }
      sig <- my_pval < 0.05
      new_row <- c(project, colnames(metadata)[colmn], my_pval, sig)
      pval_df <- rbind(pval_df, new_row)
    }
  }
}
write.csv(cor_df, 
          file.path(home_dir, "metastudies", "output", "numeric_meta_vs_read_depth_cor.csv"),
          row.names = FALSE)
write.csv(pval_df, 
          file.path(home_dir, "metastudies", "output", "factor_meta_vs_read_depth_cor.csv"),
          row.names = FALSE)

