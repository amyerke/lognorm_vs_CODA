# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making distribution of AUC of ILR balances to compare to philr

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
make_ilr_taxa_auc_df <- function(ps_obj,
                                 metadata = metadata,
                                 metadata_cols,
                                 train_index = train_index,
                                 test_index = test_index,
                                 philr_ilr_weights = philr_ilr_weights,  
                                 philr_taxa_weights = philr_taxa_weights,
                                 just_otu = FALSE){
  #Function for making random forest AUC values
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("ROCR", quietly = TRUE)) BiocManager::install("ROCR")
  library("ROCR")
  all_auc <- c()
  metadata_col <- c()
  taxa_weight <- c()
  ilr_weight <- c()
  rf_imp_seq = c()
  rf_type = c()
  rf_ntree = c()
  rf_inbag = c()
  for (ilr_w in 1:length(philr_ilr_weights)){
    for (tax_w in 1:length(philr_taxa_weights)){
      tryCatch({
        if (just_otu == TRUE){
          my_table <- ps_obj
        }else{
          my_table <- philr::philr(ps_obj@otu_table, ps_obj@phy_tree, 
                                   part.weights = philr_taxa_weights[tax_w],
                                   ilr.weights = philr_ilr_weights[ilr_w])
        }
        my_table_train <- data.frame(my_table[train_index,])
        my_table_test <- data.frame(my_table[test_index,])
        },
        error=function(cond) {
          print('Opps, an error is thrown')
          message(cond)
          print(message(cond))
        },
        warning=function(cond) {
          print('Opps, a warning is thrown')
          message(cond)
          print(message(cond))
        }
      )
      for(mta in metadata_cols){
        tryCatch(
          { 
            print(paste("starting metadata col:", mta, colnames(metadata)[mta]))
            if (any(is.na(my_table))) break
            resp_var_test <- metadata[,mta][test_index]
            resp_var_train <- metadata[,mta][train_index]
            #rf requires rownames on resp var
            names(resp_var_test) <- row.names(my_table_test)
            rf <- randomForest::randomForest(my_table_train, resp_var_train)
            pred <- predict(rf, my_table_test)
            # print(paste("pred:", pred))
            # print(paste("num factors", nlevels(resp_var_test)))
            roc_data <- data.frame(pred = pred, resp_var_test = resp_var_test)
            
            if (nlevels(resp_var_test) > 2){
              print("multilevels")
            }else{
              preds <- ROCR::prediction(as.numeric(pred), as.numeric(resp_var_test))
              auc <- ROCR::performance(preds, "auc")@y.values[[1]]
            }
            print("ROC made")
            # auc <- pROC::auc(my_roc)
            # print(paste("auc: ")
            #update all output
            metadata_col <- append(metadata_col, colnames(metadata)[mta])
            all_auc <- append(all_auc, auc)
            taxa_weight <- c(taxa_weight, philr_taxa_weights[tax_w])
            ilr_weight <- c(ilr_weight, philr_ilr_weights[ilr_w])
            ##-Section for saving random forest parameters----------------------##
            my_df <- rf$importance
            maxImp <- max(rf$importance)
            maxRow <- which(rf$importance == maxImp)
            rf_imp_seq = c(rf_imp_seq, row.names(my_df)[maxRow])
            rf_type = c(rf_type, rf$type)
            rf_ntree = c(rf_ntree, rf$ntree)
            ##------------------------------------------------------------------##
          },
          error=function(cond) {
            print('Opps, an error is thrown')
            message(cond)
            print(message(cond))
          },
          warning=function(cond) {
            print('Opps, a warning is thrown')
            message(cond)
            print(message(cond))
            }
        )
      }#for mta
      if (just_otu == TRUE) break
    }#taxa
    if (just_otu == TRUE) break
  }#ilr
  return(data.frame(all_auc,
                    metadata_col,
                    taxa_weight,
                    ilr_weight,
                    rf_imp_seq,
                    rf_type,
                    rf_ntree))
}#end function

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
library("ape")
if (!requireNamespace("randomForest", quietly = TRUE)) BiocManager::install("randomForest")
library("randomForest")
if (!requireNamespace("ROCR", quietly = TRUE)) BiocManager::install("ROCR")
library("ROCR")
if (!requireNamespace("ggpubr", quietly = TRUE)) BiocManager::install("ggpubr")
library("ggpubr")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("ggplot2")
if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
library("rgr")

print("finished loading libraries")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Jones"
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##

source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Set up constants-------------------------------------------------##
rf_cols <- 7
num_cycles <- 20
if(num_cycles < 3) stop("num_cycles should be 3 or more")
main_output_label <- paste0("alr_calibration", num_cycles)
philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")
random_seed <- 36

##-Import tables and data preprocessing-----------------------------##
asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

print("loading and munging metadata")
metadata <- read.table(file.path(home_dir, project, "patient_metadata.tsv"), 
                       sep="\t", 
                       header=TRUE, 
                       row.names = "Run", 
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
metadata <- metadata[row.names(metadata) %in% row.names(asv_table), ]
metadata$type <- droplevels(metadata$type)
metadata$type <- factor(metadata$type)

limit_num_0 <- 50
print(paste("Droping columns with only", limit_num_0, "nonzero value."))
asv_table <- asv_table[, colSums(asv_table != 0) > limit_num_0] 
dim(asv_table)
print(ncol(asv_table)/100*6/60)

##-Random num seed--------------------------------------------------##
print("Setting random seed to:", random_seed)
set.seed(random_seed)

##-Create plot data-------------------------------------------------##
all_plot_data <- data.frame(all_auc = c(),
                            metadata_col = c(),
                            taxa_weight = c(),
                            ilr_weight = c(),
                            random_batch = c(),
                            trans_group = c(),
                            rf_imp_seq = c(),
                            rf_type = c(),
                            rf_ntree = c(),
                            j_num_zeros = c())
skips <- 0
counter <- 0

train_index <- sample(x = nrow(metadata), size = 0.75*nrow(metadata), replace=FALSE)
test_index <- c(1:nrow(metadata))[!(1:nrow(metadata) %in% train_index)]
print("creating multiple iterations of the ALR")
start.time <- Sys.time()
for (alr_col in 1:ncol(asv_table)) {
  my_alr <- as.data.frame(rgr::alr(as.matrix(asv_table + 1), j = as.numeric(alr_col)))
  
  j_num_zero <- sum(asv_table[,alr_col] == 0)
  
  print(paste("counter:", counter, " generate alr data"))
  my_plot_data <- make_ilr_taxa_auc_df(ps_obj = my_alr,
                                       metadata_cols = rf_cols,
                                       metadata = metadata,
                                       train_index = train_index,
                                       test_index = test_index,
                                       philr_ilr_weights = philr_ilr_weights,
                                       philr_taxa_weights = philr_taxa_weights,
                                       just_otu = TRUE )
  my_plot_data$random_batch <- rep("None", nrow(my_plot_data))
  my_plot_data$trans_group <- rep("alr", nrow(my_plot_data))
  my_plot_data$j_num_zeros <- rep(j_num_zero, nrow(my_plot_data))
  all_plot_data <- rbind(all_plot_data, my_plot_data)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

cor(all_plot_data$all_auc, all_plot_data$j_num_zeros)

print("Make all the boxplots")
##-Make all the boxplots--------------------------------------------##
pdf(file = file.path(output_dir, "graphics", paste0(main_output_label, ".pdf")))
g <- ggplot2::ggplot(all_plot_data, aes(y = all_auc, x= j_num_zeros)) + 
  ggplot2::ggtitle(paste(project, "Type, cor:", cor(all_plot_data$all_auc, all_plot_data$j_num_zeros))) +
  ggplot2::geom_point() +
  # ggplot2::geom_hline(yintercept = 0) +
  ggplot2::theme(axis.text.x = element_text(angle = 45),
                 axis.text = element_text(size = 20)) +
  ggplot2::theme_classic() +
  # ggplot2::scale_y_discrete(labels = seq(0, 1, by = 0.2)) +
  ggplot2::ylab("AUC") +
  ggplot2::xlab("Number of zeros in denumerator column") +
  ggplot2::labs(color = "ilr weight")
print(g)
dev.off()
