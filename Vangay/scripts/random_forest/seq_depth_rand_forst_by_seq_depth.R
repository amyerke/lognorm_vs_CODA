# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing each transformation against different sequence depth to find the best one
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
munge_ref_ps <- function(ps){
  #ps must be a phyloseq object
  #requires phyloseq and ape packages to be loaded in the env
  # ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  # ps <- filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
  ps <- transform_sample_counts(ps, function(x) x+1)
  phy_tree(ps) <- makeNodeLabel(phy_tree(ps), method="number", prefix='n')
  return(ps)
}

roc_axes <- function(test_data, 
                     true_resp = resp_var_test,
                     ml_model = rf,
                     error_range = 0.10){
  # This function takes predictor variables (test_data) and tests it 
  # against the correct data (true_resp) and returns a dataframe of the 
  # true postives and false positives as ratios between 0 and 1
  # test_data is the testing variable
  # the true_resp is the correct variable
  # the error range provides an acceptable amount of error for when 
  # the numbers are close, but not exactly the same.
  groups = unique(true_resp)
  true_pos = c(0)
  false_pos = c(0)
  true_neg = c(0)
  # for categorical data
  for (grp in groups){
    if (class(true_resp) %in% c("double", "integer")){
      upper_lim_grp = grp + grp * error_range
      lower_lim_grp = grp - grp * error_range
    }
    for (rw in 1:nrow(test_data)){
      #for when data are not categoric
      if (class(true_resp) %in% c("double", "integer")){
        pred = predict(ml_model, newdata=test_data[rw,])
        upper_lim_tr = true_resp[rw] + true_resp[rw] * error_range
        lower_lim_tr = true_resp[rw] - true_resp[rw] * error_range
        decision = upper_lim_grp > pred & pred > lower_lim_grp
        truth = upper_lim_tr > pred & pred > lower_lim_tr
      }else{
        #for categoric data
        pred = predict(ml_model, newdata=test_data[rw,])
        decision = pred == grp
        truth = true_resp[rw] == grp
      }
      if (decision == truth & truth == T){
        true_pos = c(true_pos, tail(true_pos)+1)
        true_neg = c(true_neg, tail(true_neg))
        false_pos = c(false_pos, tail(false_pos))
      }else if (decision == truth & truth == F){
        true_neg = c(true_neg, tail(true_neg)+1)
        true_pos = c(true_pos, tail(true_pos))
        false_pos = c(false_pos, tail(false_pos))
      }else if (decision != truth & truth == T){
        true_neg = c(true_neg, tail(true_neg))
        true_pos = c(true_pos, tail(true_pos))
        false_pos = c(false_pos, tail(false_pos)+1)
      }
    }#for (rw in 1:nrow(test_data)){
  }#for (grp in groups){
  print(paste("max(false_pos):", max(false_pos), "max(true_pos):", max(true_pos), "any na:", any(is.na(true_pos))))
  if( max(true_neg) != 0){
    true_neg = true_neg/max(true_neg)
  }
  if( max(true_pos) != 0){
    true_pos = true_pos/max(true_pos)
  }
  if( max(false_pos) != 0){
    false_pos = false_pos/max(false_pos)
  }
  return(data.frame(true_pos, false_pos))
}#end roc_data

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library("compositions")
library("phyloseq")
library("vegan")
library("DESeq2")
library("philr")
library("ape")
library("ALDEx2")
library("ggplot2")
library("RColorBrewer")
#set color palette
palette( brewer.pal(7,"Accent") )
library("randomForest")

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
asv_table <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds"))

ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
asv_tax <- readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))

metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(asv_table),]
metadata$type <- droplevels(metadata$type)

##-----------------Create training/testing sets---------------------##
set.seed(36)
train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]


#Plot shannon diversity against log10(total_seqs)
#Plot other normalization methods: otu log 100 and alr and clr 
# my_ds_names <- c( "raw seqs", "clr(raw seqs)", "lognorm raw seqs", "philr ref", "DESeq2", "ALDEx2.clr")
my_ds_names <- c( "raw seqs", "clr(raw seqs)", "lognorm raw seqs", "DESeq2", "ALDEx2.clr")

# my_ds_names <- c( "raw seqs", "clr(raw seqs)", "lognorm raw seqs")

min_seq_depths <- c(0, 500, 1000, 5000, 10000, 20000, 40000)
mds_depth <- 5
mta = "type"

total_seqs <- rowSums(asv_table)
total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))

pdf(file = file.path(output_dir, "graphics","roc_seq_depth_by_seq_dep.pdf"))
my_rocs <- list()
for(s in 1:length(min_seq_depths)){
  seq_d <- min_seq_depths[s]#new sequencing depth
  sd_filt_asv <- asv_table[total_seqs$total_seqs >= seq_d,]#dataset 1
  
  print(paste("sd_filtered dim:", paste(dim(sd_filt_asv))))
  safe_rns <- intersect(row.names(ref_ps@otu_table), row.names(sd_filt_asv)) #rows for this iterate
  new_metadata <- metadata[safe_rns,]
  train_index <- sample(x = nrow(new_metadata), size = 0.75*nrow(new_metadata), replace=FALSE)
  test_index <- c(1:nrow(new_metadata))[!(1:nrow(new_metadata) %in% train_index)]
  
  
  my_clr <- compositions::clr(sd_filt_asv)#dataset 2
  new_tree <- phyloseq::prune_taxa(colnames(sd_filt_asv), ref_ps@phy_tree)#update tree for new phyloseq obj
  new_ref_ps <- phyloseq::prune_samples(safe_rns, ref_ps) #remove non-safe rows from ps
  new_ref_ps <- munge_ref_ps(new_ref_ps)
  print(paste("new dim ref ps:", dim(data.frame(new_ref_ps@otu_table))))
  #create DESeq2 dtaset from new ref ps
  new_DESeq2 <- phyloseq::phyloseq_to_deseq2(new_ref_ps, design= ~ 1)#dataset 5
  new_DESeq2 <- DESeq2::estimateSizeFactors(new_DESeq2)
  new_DESeq2 <- t(counts(new_DESeq2, normalized=T))

  print(paste("new DSeq:", paste(dim(new_DESeq2))))

  # ref_philr <- philr::philr(new_ref_ps@otu_table, new_ref_ps@phy_tree,
  #                           part.weights='enorm.x.gm.counts',
  #                           ilr.weights='blw.sqrt')#dataset 4
  # print(paste("made new philr", dim(as.data.frame(ref_philr))))

  ln_asv <- lognorm(sd_filt_asv)#dataset 6

  ald <- ALDEx2::aldex.clr(sd_filt_asv, mc.samples=12, denom="all", verbose=F)
  ald <-  data.frame(ald@analysisData)
  print(paste("size of ald:", object.size(ald)))
  print(paste("ald dim:", paste(dim(ald))))
  # my_datasets <- list(sd_filt_asv, my_clr,  ln_asv, new_DESeq2, ald)
  my_datasets <- list(sd_filt_asv, my_clr,  ln_asv, new_DESeq2, ald)
  # my_datasets <- list(sd_filt_asv, my_clr,  ln_asv, ref_philr, new_DESeq2, ald)
  print(paste("finished seq depth filter:", s))
  
  for( ds in 1:length(my_datasets)){
    print(my_ds_names[ds])
    my_table <- as.data.frame(my_datasets[ds])
    
    resp_var_test <- new_metadata[,mta][test_index]
    resp_var_train <- new_metadata[,mta][train_index]
    my_table_train <- my_table[train_index,]
    my_table_test <- my_table[test_index,]
    print(dim(my_table))
    rf <- randomForest::randomForest(
      my_table_train, as.factor(resp_var_train)
    )
    my_roc = roc_axes(test_data = my_table_test,
                      true_resp = resp_var_test,
                      ml_model = rf,
                      error_range = 0)
    my_rocs[[ds]] = my_roc

  }#for ds
  par(bg = 'grey96')
  plot(true_pos ~ false_pos,
       data = as.data.frame(my_rocs[1]),
       type = "l",
       xlab = "False positives",
       ylab = "True positives",
       col = palette()[1],
       main = paste(project, "stool/swab","seq_dept:", seq_d),
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
  plot.new()
  legend('bottomright', 
         legend = my_ds_names, 
         col = palette(), 
         pch = 15,
         cex=0.5)
  
}#for sd
dev.off()

