# Author: Aaron Yerke
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
project <- "Jones"
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))


library(SNFtool)
##-Set variables for SNFtool----------------------------------------##
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
SNF_T = 10; 	# Number of Iterations, usually (10~20)


##-Import tables and data preprocessing-----------------------------##
philr_trans_ref <- read.table(file.path(output_dir, "tables", "ref_tree_ps_philr_transform.csv"),
                              sep = ",",
                              header = TRUE)

philr_trans_upgma <- read.table(file.path(output_dir, "tables", "philr_denovo_tree_UPGMA.csv"),
                                sep = ",",
                                header = TRUE)

##-The SNFtool process from the readme------------------------------##
Dist1 <- dist2(as.matrix(philr_trans_ref),as.matrix(philr_trans_ref));
Dist2 <- dist2(as.matrix(philr_trans_upgma),as.matrix(philr_trans_upgma));

W1 <- affinityMatrix(Dist1, K, alpha)
W2 <- affinityMatrix(Dist2, K, alpha)

W <- SNF(list(W1,W2), K, SNF_T)

C = 2 					# number of clusters
group = spectralClustering(W, C); 	# the final subtypes information

displayClusters(W, group);
SNFNMI = calNMI(group, truelabel)

ConcordanceMatrix = concordanceNetworkNMI(list(W, W1,W2),C)

my_prcmp <- prcomp(W, 
                   center = TRUE,
                   scale = TRUE)
myPCA <- data.frame(my_prcmp$x)

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(philr_trans_ref), ]



for (m in 1:ncol(metadata)) {
  plot(myPCA$PC1, 
       myPCA$PC2,
       col = metadata$sex
  )
}


