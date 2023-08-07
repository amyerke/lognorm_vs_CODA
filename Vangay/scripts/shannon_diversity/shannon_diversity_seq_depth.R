# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against different sequence depth to find the best one
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("reshape2", quietly = TRUE)) BiocManager::install("reshape2")
library("vegan")
library("ggplot2")
library("reshape2")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Vangay"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Import tables and data preprocessing-----------------------------##
asv_table <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds"))
metadata <- read.table(file.path(home_dir, project, "fullMetadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "run_accession",
                       stringsAsFactors=TRUE)
# wanted_cols <- "Recruitment.Date	Sample.Date	Recruitment.Location	Researcher	Sub.Study	Birth.Year	Age	Public.Housing	Medical.Assistance	Children.Free.Lunch	Highest.Education	Ethnicity	Religion	Birth.Location	Type.Birth.Location	Arrival.in.US	Years.in.US	Location.before.US	Type.location.before.US	Years.lived.in.Location.before.US	Tobacco.Use	Alcohol.Use	Height	Weight	Waist	BMI	BMI.Class	Medications	Breastfed	Years.Breastfed"
# wanted_cols <- "Researcher	Sub.Study	Public.Housing	Medical.Assistance	Children.Free.Lunch	Highest.Education	Ethnicity	Religion	Birth.Location	Type.Birth.Location	Location.before.US	Type.location.before.US	Tobacco.Use	Alcohol.Use	BMI.Class\tBreastfed"
# wanted_cols <- unlist(strsplit(wanted_cols, "\t"))
# metadata <- metadata[,wanted_cols]
metadata <- metadata[ lapply( metadata, function(x) sum(is.na(x)) / length(x)) < 0.05 ]
metadata <- metadata[ lapply( metadata, function(x) (length(unique(x)))) > 1 ]
metadata <- metadata[row.names(metadata) %in% row.names(asv_table), ]


metad_cols <- which(unlist(lapply(metadata, is.numeric)))

#Plot shannon diversity against log10(total_seqs)
#Plot other normalization methods: otu log 100 and alr and clr 
min_seq_depths <- c(0, 500, 1000, 5000, 10000, 20000, 25000, 30000)
mds_depth <- 5

total_seqs <- rowSums(asv_table)
total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))

kend <- vector(mode = "numeric", length = length(min_seq_depths) * mds_depth)
perma_r2 <- vector(mode = "numeric", length = length(min_seq_depths) * mds_depth)
mds_lev <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
seq_depth <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
var_exp <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
PCAi_shannon_div_spear_cor <- vector(mode = "integer", length =length(min_seq_depths) * mds_depth)
samples_left <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
taxa_left <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)
zero_count <- vector(mode = "integer", length = length(min_seq_depths) * mds_depth)

meta_cor <- data.frame(matrix(ncol = length(metad_cols), nrow = length(min_seq_depths) * mds_depth))
colnames(meta_cor) <- names(metad_cols)

meta_p <- data.frame(matrix(ncol = length(metad_cols), nrow = length(min_seq_depths) * mds_depth))
colnames(meta_p) <- names(metad_cols)

counter <- 1
for(s in 1:length(min_seq_depths)){
  seq_d <- min_seq_depths[s]#new sequencing depth
  sd_filt_asv <- asv_table[total_seqs$total_seqs >= seq_d,]#dataset 1
  cat(paste("sd_filtered dim:", dim(sd_filt_asv)))
  shan_div <- vegan::diversity(sd_filt_asv)
  print(paste("finished seq depth filter:", s))
  
  ##-Create a PCA-----------------------------------------------------##
  my_prcmp <- prcomp(sd_filt_asv, 
                     center = TRUE,
                     rank = mds_depth)#,
  # scale = TRUE)
  ##-Extract PCA matrix and convert to dataframe----------------------##
  myPCA <- data.frame(my_prcmp$x)
  my_var_exp <- my_prcmp$sdev^2/sum(my_prcmp$sdev^2)
  shan_div <- vegan::diversity(sd_filt_asv)
  if(nrow(sd_filt_asv) > 3){
    for (md in 1:mds_depth){
      for (m in 1:length(metad_cols)){
        my_cor <- cor.test(shan_div, metadata[row.names(sd_filt_asv),metad_cols[m]],
                           method = "spearman")
        meta_cor[counter, names(metad_cols)[m]] <- my_cor$estimate
        meta_p[counter, names(metad_cols)[m]] <- my_cor$p.value
      }
      mds_lev[counter] <- md
      seq_depth[counter] <- seq_d
      var_exp[counter] <- my_var_exp[md]
      PCAi_shannon_div_spear_cor[counter] <- cor(shan_div, myPCA[,md], method = "spearman")
      samples_left[counter] <- nrow(sd_filt_asv)
      taxa_left[counter] <- sum(colSums(sd_filt_asv) > 0)
      zero_count[counter] <- sum(sd_filt_asv == 0)
      counter <- counter + 1
    }
  }#if statement
}
result_df <- data.frame( mds_lev, seq_depth, var_exp, PCAi_shannon_div_spear_cor, samples_left, zero_count, taxa_left)
print("created resulting DF")

comb_result_df <- cbind(result_df, meta_cor)

write.table(comb_result_df, 
            file = file.path(output_dir, "tables", paste0(project, "_PCA_shandiv_filt_results_cor.csv")),
            sep = ",")

comb_result_df <- read.table(file = file.path(output_dir, "tables", paste0(project, "_PCA_shandiv_filt_results_cor.csv")),
                        sep = ",")

melted_results <- reshape2::melt(comb_result_df,
                       variable.name = "Metadata",
                       measure.vars = c(colnames(meta_cor), "PCAi_shannon_div_spear_cor"))

pdf(file = file.path(output_dir, "graphics", "shan_div_cor_PCA12345_line.pdf"))
for (i in 1:max(result_df$mds_lev)){
  g <- ggplot2::ggplot(melted_results[melted_results$mds_lev == i, ],
                       aes(x=seq_depth, y=value^2, group = Metadata)) +
    ggplot2::geom_point(aes(color = factor(Metadata))) +
    ggplot2::geom_line(aes(color = factor(Metadata))) +
    ggplot2::ggtitle(paste0("PCA",i," Seq depth vs Rsqd of Shannon Diversity vs metadata")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("Spearman corr") +
    ggplot2::labs(fill = "Metadata") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
}
dev.off()

comb_result_df <- cbind(result_df, meta_p)
write.table(comb_result_df, 
            file = file.path(output_dir, "tables", paste0(project, "_PCA_shandiv_filt_results_p.csv")),
            sep = ",")

comb_result_df <- read.table(file = file.path(output_dir, "tables", paste0(project, "_PCA_shandiv_filt_results_p.csv")),
                             sep = ",")

melted_results <- reshape2::melt(comb_result_df,
                                 variable.name = "Metadata",
                                 measure.vars = c(colnames(meta_cor)))


pdf(file = file.path(output_dir, "graphics", "shan_div_p_PCA12345_line.pdf"))
  g <- ggplot2::ggplot(melted_results,
                       aes(x=seq_depth, y=value, group = Metadata)) +
    ggplot2::geom_point(aes(color = factor(Metadata))) +
    ggplot2::geom_line(aes(color = factor(Metadata))) +
    ggplot2::ggtitle(paste0("PCA",i," Read depth vs P-val of Shannon Diversity vs metadata")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("Pvale") +
    ggplot2::labs(fill = "Metadata") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  print(g)
  
  g <- ggplot2::ggplot(melted_results,
                       aes(x=seq_depth, y=log10(value), group = Metadata)) +
    ggplot2::geom_point(aes(color = factor(Metadata))) +
    ggplot2::geom_line(aes(color = factor(Metadata))) +
    ggplot2::ggtitle(paste0("PCA",i," Read depth vs log20(P-val) of Shannon Diversity vs metadata")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("log10(pval)") +
    ggplot2::labs(fill = "Metadata") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  print(g)
  
  g <- ggplot2::ggplot(result_df,
                       aes(x=seq_depth, y=samples_left)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::ggtitle(paste0("Number samples vs read depth filt")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("Num samples") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  print(g)

dev.off()



