# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making distribution of AUC of ILR balances to compare to philr
# Requires UPGMA_Tree, Silva_tree, Iqtree, and metadata

rm(list = ls()) #clear workspace

# --------------------------------------------------------------------------
print("Defining functions")
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
print("Loading dependencies")
# --------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
library("ape")
if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
library("philr")
if (!requireNamespace("randomForest", quietly = TRUE)) BiocManager::install("randomForest")
library("randomForest")
if (!requireNamespace("ggpubr", quietly = TRUE)) BiocManager::install("ggpubr")
library("ggpubr")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
library("phyloseq")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("ggplot2")
if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
library("rgr")
if (!requireNamespace("pROC", quietly = TRUE)) BiocManager::install("pROC")
library("pROC")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
if (!requireNamespace("MLmetrics", quietly = TRUE)) BiocManager::install("MLmetrics")
library("MLmetrics")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
library("optparse")
# --------------------------------------------------------------------------
print("Reading cml arguments")
# --------------------------------------------------------------------------
option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git',"lognorm_vs_CODA"), 
                        help="dataset dir path"),
  optparse::make_option(c("-p", "--project"), type="character", default="Jones", 
                        help="project folder")
);

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

# --------------------------------------------------------------------------
print("Establishing directory layout and other constants.")
# --------------------------------------------------------------------------
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
input_dir <- file.path(output_dir, "tables", "r_v_py_train_test_tables")
base::setwd(input_dir)

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Set up constants-------------------------------------------------##
main_output_text <- "from_file_random_forest_score_R"
main_output_fn <- paste0(main_output_text, ".csv")
main_output_fpath <- file.path(output_dir, "tables", main_output_fn)

random_seed <- 36
main_header <- "all_score, metadata_col, rf_imp_se, rf_type, rf_ntree, trans_group, cycle"

##-Import tables and data preprocessing-----------------------------##
print(paste("Reading files at ", input_dir))
##-Import tables and data preprocessing-----------------------------##
my_files <- list.files(path = input_dir)
metadata <- c()
iteration_min <- NaN
iteration_max <- NaN
for (filen in my_files){
  # print(filen)
  my_splits = unlist(strsplit(filen, split="(_)", fixed=TRUE))
  meta <- my_splits[1]
  iteration <- base::as.integer( my_splits[2] )
  if (! meta %in% metadata){
    metadata <- c(metadata, meta)
  }
  # if iteration_min is NaN, then make iteration, otherwise leave it iteration_min
  iteration_min <- base::ifelse(base::is.nan(iteration_min), iteration, iteration_min)
  iteration_max <- base::ifelse(base::is.nan(iteration_max), iteration, iteration_max)
  iteration_min <- base::ifelse(iteration < iteration_min, iteration, iteration_min)
  iteration_max <- base::ifelse(iteration > iteration_max, iteration, iteration_max)
}

# ##-Import tables and data preprocessing-----------------------------##
# print(paste("Iterating through files at ", input_dir))
# ##-Import tables and data preprocessing-----------------------------##
# print(paste("Initilizing", main_output_fpath, "."))
# print(paste("File header: ", main_header))
# cat(main_header,
#     file = main_output_fpath,
#     append=FALSE)

# for (feat in metadata){
#   for (i in iteration_min:iteration_max){
#     print(paste(feat, i))
#     pred_train <- data.frame(data.table::fread(file = paste(feat, i,"pred", "train.csv", sep = "(_)"),
#                                             header=TRUE, data.table=FALSE), row.names = 1)
#     print(paste("pred_train:", nrow(pred_train)))
#     pred_test <- data.frame(data.table::fread(file = paste(feat, i,"pred", "test.csv", sep = "(_)"),
#                                                header=TRUE, data.table=FALSE), row.names = 1)
#     print(paste("pred_test:", nrow(pred_test)))
#     resp_train <- read.csv(paste(feat, i,"resp", "train.csv", sep = "(_)"), 
#                            header = TRUE,
#                            row.names = 1,
#                            stringsAsFactors = TRUE,
#                            check.names = FALSE)
#     print(paste("resp_train:", nrow(resp_train)))
#     resp_test <- read.csv(paste(feat, i,"resp", "test.csv", sep = "(_)"), 
#                           header = TRUE,
#                           row.names = 1,
#                           stringsAsFactors = TRUE,
#                           check.names = FALSE)
#     feat_sub = gsub("∕","/",feat)
#     if (class(resp_train[,feat_sub]) == "character"){
#       resp_train[,feat_sub] <- factor(resp_train[,feat_sub])
#       resp_test[,feat_sub] <- factor(resp_test[,feat_sub])
#     }
#     if (is.factor(resp_train[,feat_sub])){
#       # print(paste("nlevl resp_test:" length(levels(resp_var_test)))
#       all_levels <- base::union(levels(pred_train), levels(resp_train[,feat_sub]))
#       levels(pred_train) <- all_levels#hack for when the levels are different
#       levels(resp_train[,feat_sub]) <- all_levels
#     }
#     print(paste("resp_test:", nrow(resp_test)))
#     rf <- randomForest::randomForest(pred_train, resp_train[row.names(pred_train),feat_sub])
#     print("made rf")
#     pred <- predict(rf, pred_test)
#     print("pred")
#     if (is.factor(resp_train[,feat_sub])){
#       # print(paste("nlevl resp_test:" length(levels(resp_var_test)))
#       all_levels <- base::union(levels(pred), levels(resp_test[,feat_sub]))
#       levels(pred) <- all_levels#hack for when the levels are different
#       levels(resp_test[,feat_sub]) <- all_levels
#     }
#     my_test <- resp_test[row.names(pred_test),feat_sub]
#     names(my_test) <- row.names(pred_test)
#     score <- MLmetrics::Accuracy(pred, my_test)
#     print(paste("score:", score))
#     my_df <- rf$importance
#     maxImp <- max(rf$importance)
#     maxRow <- which(rf$importance == maxImp)
    
#     # Check its existence
#     if (file.exists(main_output_fpath)) {
#       print(paste0("Writing output to ", main_output_fpath, " ."))
#       # main_header <- "all_score,	metadata_col,	rf_imp_se, rf_type, rf_ntree, trans_group, random_batch, cycle"
#       cat(paste(paste0("\n", score), feat, row.names(my_df)[maxRow], rf$type, #ilr_weight,	rf_imp_se, rf_type,
#                 rf$ntree, "raw_DADA2", i, #rf_ntree, trans_group, cycle
#                 sep = ","),
#           file = main_output_fpath,
#           append=TRUE)
#     }
#   }
# }



# --------------------------------------------------------------------------
print("Reshaping data from long to wide format")
# --------------------------------------------------------------------------
print("Reading in data from file.")
all_plot_data <- data.frame(read.table(file = main_output_fpath,
                                       sep = ",", header = TRUE))
print("Building wide df.")
data_cols <- unlist(lapply(iteration_min:iteration_max, function(x) paste0("split",x)))
wide_header <- c("metadata_col", "trans_group", "color", data_cols)
wide_result_df <- data.frame(matrix(nrow = length(unique(all_plot_data$metadata_col)) * length(unique(all_plot_data$trans_group)),
                                    ncol = length(wide_header)))
colnames(wide_result_df) <- wide_header
wide_result_df$metadata_col <- rep(unique(all_plot_data$metadata_col), length(unique(all_plot_data$trans_group)))
wide_result_df$trans_group <- rep(unique(all_plot_data$trans_group), each = length(unique(all_plot_data$metadata_col)))

for (ro in 1:nrow(all_plot_data)){
  my_split <- paste0("split", all_plot_data[ro,"cycle"])
  wide_result_df[wide_result_df$metadata_col == all_plot_data[ro,"metadata_col"] & wide_result_df$trans_group == all_plot_data[ro,"trans_group"], my_split] = all_plot_data[ro,"all_score"]
}
print(wide_result_df)

write.table(wide_result_df, row.names = FALSE, sep = ",",
            file = file.path(output_dir, "tables", paste0("wide_", main_output_text, ".csv")))

print("Completed script!")

