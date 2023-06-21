# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for ROC of a all the spreadsheets
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-------------------Load Depencencies------------------------------##
if (!requireNamespace("randomForest", quietly = TRUE)){
  install.packages("randomForest")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)
#set color palette
palette( brewer.pal(7,"Accent") )
library(randomForest)
library(compositions)
##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Vangay"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

meta_col <- "host_phenotype"

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "r_libraries", "statistical_functions.R"))
source(file.path(home_dir, "r_libraries", "table_manipulations.R"))

SNF_fun <- function(df1, df2) {
  library(SNFtool)
  ##-Set variables for SNFtool----------------------------------------##
  K = 20;		# number of neighbors, usually (10~30)
  alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
  SNF_T = 10; 	# Number of Iterations, usually (10~20)
  
  ##-The SNFtool process from the readme------------------------------##
  Dist1 <- dist2(as.matrix(df1),as.matrix(df1));
  Dist2 <- dist2(as.matrix(df2),as.matrix(df2));
  
  W1 <- affinityMatrix(Dist1, K, alpha)
  W2 <- affinityMatrix(Dist2, K, alpha)
  
  W <- SNF(list(W1,W2), K, SNF_T)
  
  C <- 2 					# number of clusters
  group <- spectralClustering(W, C); 	# the final subtypes information
  
  # displayClusters(W, group)
  # SNFNMI <- calNMI(group, truelabel)
  
  ConcordanceMatrix <- concordanceNetworkNMI(list(W, W1,W2),C)
  detach(package:SNFtool)
  
  # return(data.frame(ConcordanceMatrix))
  return(data.frame(W))
}

##-Import tables and data preprocessing-----------------------------##
philr_trans_ref <- read.table(file.path(output_dir, "tables", "ref_tree_ps_philr_transform.csv"),
                              sep = ",",
                              header = TRUE)

philr_trans_upgma <- read.table(file.path(output_dir, "tables", "philr_denovo_tree_UPGMA.csv"),
                                sep = ",",
                                header = TRUE)

otu_tab <- read.table(file.path(output_dir, "tables", "genus_otu_table.csv"),
                      sep = ",",
                      header = TRUE)

otu_ref_tree_combo <- SNF_fun(philr_trans_ref, otu_tab)

ln_otu_tab <- read.table(file.path(output_dir, "tables", "lognorm_genus_otu_table.csv"),
                      sep = ",",
                      header = TRUE)

ln_otu_philr_upgma <- SNF_fun(ln_otu_tab, philr_trans_upgma)

# otu_ratio <- simple_ratio_transform(otu_tab)

# asv_table <- read.table(file.path(output_dir, "tables", "ForwardReads_DADA2.txt"),
#                         sep = "\t",
#                         header = TRUE)

# my_alr <- as.data.frame(compositions::alr(as.matrix(asv_table)))

# my_clr <- as.data.frame(compositions::clr(as.matrix(asv_table)))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
# print(paste("metadata dim:", dim(metadata)))
metadata <- metadata[row.names(metadata) %in% row.names(philr_trans_ref), ]
# print(paste("metadata dim:", dim(metadata)))
print(sapply(metadata,class))
# metadata$Genotype <- as.factor(metadata$Genotype)
# metadata$Visit <- as.factor(metadata$Visit)
# metadata$Treatment <- as.factor(metadata$Treatment)
# metadata$sex <- as.factor(metadata$sex)

# my_datasets = list(asv_table, otu_tab, otu_ratio, my_alr, my_clr, philr_trans_ref, philr_trans_upgma, cbind(otu_ratio, philr_trans_ref), cbind(otu_tab, otu_ratio))
# my_ds_names = c("ASV", "OTU", "OTU ratio", "alr", "clr", "ref tree philR", "UPGMA philr trans", "OTU + ref tree philR", "OTU + OTU ratio")
# my_datasets <- list(asv_table, ln_otu_tab, my_alr, my_clr, philr_trans_ref, philr_trans_upgma )
# my_ds_names <- c("Pre-filtered ASV", "lognorm OTU", "Pre-filter alr", "Pre-filter clr", "ref tree philR", "UPGMA philr trans")
# my_datasets <- list(asv_table, ln_otu_tab, philr_trans_ref, philr_trans_upgma )
# my_ds_names <- c("Pre-filtered ASV", "lognorm OTU", "ref tree philR", "UPGMA philr trans")

my_datasets <- list(otu_ref_tree_combo, ln_otu_philr_upgma, philr_trans_ref, philr_trans_upgma )
my_ds_names <- c("OTU + ref Tree", "LN OTU + UPGMA", "ref tree philR", "UPGMA philr trans")


##-----------------Create training/testing sets---------------------##
set.seed(36)
train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]

# pdf(file = file.path(output_dir, "graphics","roc_all_metadata3.pdf"))
my_rocs <- list()
for(mta in 1:ncol(metadata)){
  resp_var_test <- metadata[,mta][test_index]
  resp_var_train <- metadata[,mta][train_index]
  for( ds in 1:length(my_datasets)){
    my_table <- as.data.frame(my_datasets[ds])
    my_table_train <- my_table[train_index,]
    my_table_test <- my_table[test_index,]
    
    rf <- randomForest(
      my_table_train, resp_var_train
    )
    my_roc = roc_axes(test_data = my_table_test,
                      true_resp = resp_var_test
    )
    print(paste(ds, mta))
    my_rocs[[ds]] = my_roc
  }#for ds
  par(bg = 'grey96')
  plot(true_pos ~ false_pos,
       data = as.data.frame(my_rocs[1]),
       type = "l",
       xlab = "False positives",
       ylab = "True positives",
       col = palette()[1],
       main = paste(project, colnames(metadata)[mta]),
       ylim = c(0,1))
  for(i in 2:length(my_ds_names)){
    lines(true_pos ~ false_pos,
          data = as.data.frame(my_rocs[i]),
          col = palette()[i],
          type = "l")
  }
  abline(
    a = 0,
    b = 1)
  legend('bottomright', 
         legend = my_ds_names, 
         col = palette(), 
         pch = 15,
         cex=0.5)
}#for mta
# dev.off()

