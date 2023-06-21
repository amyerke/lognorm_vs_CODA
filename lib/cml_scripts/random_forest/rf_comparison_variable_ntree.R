# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making distribution of AUC of ILR balances to compare to philr
# Requires UPGMA_Tree, Silva_tree, Iqtree, and metadata

rm(list = ls()) #clear workspace

# --------------------------------------------------------------------------
print("Defining functions")
# --------------------------------------------------------------------------
add_PhILR_dfs_to_table <- function(lst,
                                   root_folder,
                                   base_fn,
                                   # philr_part_weights = c("anorm","enorm"),
                                   # philr_ilr_weights = c("blw.sqrt","mean.descendants"),
                                   philr_part_weights = c("anorm"),
                                   philr_ilr_weights = c("blw.sqrt"),
                                   color = "w"){
  if (!dir.exists(root_folder)){
    print(paste0(root_folder,"does not exist. Use PhILR_random_trees_and_counts_tables.R to create it."))
    # break
  }
  for (ilr_w in philr_ilr_weights){
    for (tax_w in philr_part_weights){
      my_label <- paste(base_fn, ilr_w, tax_w, sep = "_")
      table_fn <- paste0(my_label, ".csv")
      lst[[length(lst) + 1]] <- c(my_label, file.path(root_folder, table_fn), ',', color)
    }
  }
  return(lst)
}

add_random_tree_PhILRs_to_table  <- function(lst, 
                                             root_folder, 
                                             base_fn, 
                                             # philr_part_weights = c("anorm","enorm"),
                                             # philr_ilr_weights = c("blw.sqrt","mean.descendants"),
                                             philr_part_weights = c("anorm"),
                                             philr_ilr_weights = c("blw.sqrt"),
                                             color = "w", 
                                             num_rand_trees = 10){
  print(paste0("Adding random trees from ", base_fn, "."))
  if (!dir.exists(root_folder)){
    print(paste0(root_folder,"does not exist. Use PhILR_random_trees_and_counts_tables.R to create it."))
    # break
  }
  for (rand in 1:num_rand_trees){
    for (ilr_w in philr_ilr_weights){
      for (tax_w in philr_part_weights){
        my_label <- paste(paste0("Shuffle", rand, "_PhILR"), base_fn, ilr_w, tax_w, sep = "_")
        table_fn <- paste0(my_label, ".csv")
        lst[[length(lst) + 1]] <- c(my_label, file.path(root_folder, table_fn), ',', color)
      }
    }
  }
  return(lst)
}

df_factory <- function(my_path, my_sep){
  # Use for building df from "tables" list in random forest loop
  df <- read.table(my_path, sep=my_sep, 
                   header=0, row.names=0,
                   check.names = FALSE,
                   stringsAsFactors=TRUE)
  return(df)
}

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
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default=NULL, 
                        help="project folder", metavar="character"),
  optparse::make_option(c("-m", "--metadata"), type="character", default=NULL,
                        help="metadata file path with filename", metavar="character"),
  optparse::make_option(c("-l", "--metadata_delim"), type="character", default="\t",
                        help="metadata file deliminator", metavar="character"),
  optparse::make_option(c("-r", "--metadata_rowname"), type="character", default=NULL,
                        help="metadata file row to use for row names", metavar="character"),
  optparse::make_option(c("-n", "--num_cycles"), type="numeric", default=20,
                        help="Number of times to shuffle data and run loop again"),
  optparse::make_option(c("-s", "--start_ntree"), type="numeric", default=100,
                        help="Starting number of ntree"),
  optparse::make_option(c("-e", "--end_ntree"), type="numeric", default=2000,
                        help="Final number of ntree"),
  optparse::make_option(c("-i", "--ntree_increment"), type="numeric", default=200,
                        help="Number of times to shuffle data and run loop again")
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

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Set up constants-------------------------------------------------##
num_cycles <- opt$num_cycles
if(num_cycles < 3) stop("num_cycles should be 3 or more")
main_output_text <- "random_forest_score_R_var_ntree"
main_output_label <- paste0(main_output_text, num_cycles)
main_output_fn <- paste0(main_output_label, ".csv")
main_output_fpath <- file.path(output_dir, "tables", main_output_fn)
# philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
# philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")
philr_taxa_weights <- c("enorm")
philr_ilr_weights <- c("mean.descendants")
random_seed <- 36
main_header <- "all_score, metadata_col, rf_imp_se, rf_type, rf_ntree, trans_group, cycle, ntree"

##-Import tables and data preprocessing-----------------------------##
print("loading and munging metadata")
##-Import tables and data preprocessing-----------------------------##
metadata <- read.table(opt$metadata,
                       sep=opt$metadata_delim,
                       header=TRUE,
                       row.names = opt$metadata_rowname,
                       check.names = FALSE,
                       stringsAsFactors=TRUE)
print("Metadata columns:")
print(sapply(metadata, class))
print("Removing non-factor columns from metadata to decrease run time.")
metadata <- metadata[,sapply(metadata, is.factor)]
metadata <- Filter(function(x) length(unique(na.omit(x))) < 6, metadata)#hack to remove date columns if they get counted as factors
print("Metadata columns:")
print(colnames(metadata))
# df[,-which(sapply(df, class) == "factor")]
rf_cols <- 1:ncol(metadata)#hack so I don't have to fix this in the function


tables <- list()
tables[[length(tables) + 1]] <- c("Raw_DADA2",file.path(output_dir, "tables", "ForwardReads_DADA2.txt"),"\t","r")
# tables[[length(tables) + 1]] <- c("lognorm_DADA2", file.path(output_dir, "tables", "lognorm_dada2.csv"), ",", "y")
# # tables[[length(tables) + 1]] <- c("lognorm_HashSeq", file.path(output_dir,"tables", "lognorm_hashseq.csv"), ",", "y")
# tables[[length(tables) + 1]] <- c("alr_DADA2", file.path(output_dir, "tables", "alr_asv.csv"), ",", "g")
# # tables[[length(tables) + 1]] <- c("alr_HashSeq", file.path(output_dir,"tables", "alr_hashseq.csv"), ",", "g")
# tables[[length(tables) + 1]] <- c("clr_DADA2", file.path(output_dir, "tables", "clr_asv.csv"), ",", "g")
# # tables[[length(tables) + 1]] <- c("clr_HashSeq", file.path(output_dir,"tables", "clr_hashseq.csv"), ",", "m")
# tables[[length(tables) + 1]] <- c("Silva_DADA2", file.path(output_dir,"tables", "Silva_DADA2", "Silva_DADA2.csv"), ",", "#64baeb")

print(tables)

print(paste("Initilizing", main_output_fpath, "."))
print(paste("File header: ", main_header))
cat(main_header,
    file = file.path(output_dir, "tables", paste0("long_",main_output_label, ".csv")),
    append=FALSE)

rowname_table <- data.frame(data.table::fread(file = tables[[1]][2],#this is a hack to make the train and test index work
                                              header=TRUE, data.table=FALSE),#the row names of this table should be
                            row.names = 1)#available in all of the other tables or there will be an error

print("Initializing wide df.")
len_ntree <- length(seq(opt$start_ntree, opt$end_ntree, by=opt$ntree_increment))
len_meta <- length(rf_cols)
len_trans <- length(tables)
data_cols <- unlist(lapply(1:num_cycles, function(x) paste0("split",x)))
wide_header <- c("metadata", "dataset", "color", "ntrees", data_cols)
wide_result_df <- data.frame(matrix(nrow = len_meta * len_trans * len_ntree, ncol = length(wide_header)))
colnames(wide_result_df) <- wide_header
counter <- 0
print(paste("Counter:", counter, " entering main loop"))

for (counter in 1:num_cycles) {
  row_count <- 1
#Create training/testing sets-------------------------------------##
  train_index <- row.names(rowname_table)[sample(x = nrow(rowname_table), size = 0.75*nrow(rowname_table), replace=FALSE)]
  test_index <- row.names(rowname_table)[c(1:nrow(rowname_table))[!(1:nrow(rowname_table) %in% train_index)]]

  for(mta in rf_cols){
    print(paste("starting metadata col:", mta, colnames(metadata)[mta]))
    print(paste("Finding rows where", colnames(metadata)[mta], "has NA or Empty values from training and testing selectors."))
    na_or_empty_index <- which(is.na(metadata[,mta]) | metadata[,mta] == "")
    na_or_empty_rows <- row.names(metadata)[na_or_empty_index]
    print(length(train_index))
    train_index <- setdiff(train_index, na_or_empty_rows)
    print(train_index)
    test_index <- setdiff(test_index, na_or_empty_rows)
    print(length(train_index))
    # tryCatch({
      for (tabl in tables){
        print(tabl)
        transf_label <- tabl[1]
        my_table <- data.frame(data.table::fread(file = tabl[2],
                                                header=TRUE, data.table=FALSE),
                               row.names = 1)
        if(setequal(row.names(metadata), row.names(my_table)) == FALSE){
          print(paste0("These samples are in the metdata, but not ", transf_label,": ", setdiff(row.names(metadata), row.names(my_table))))
          print(paste0("These samples are in the ", transf_label, " but not metadata: ", setdiff(row.names(my_table), row.names(metadata))))
          # warning(paste("Metadata dataframe and", transf_label, "dataframe must have the same rows (order is not important)."))
        }
        my_table_train <- my_table[row.names(my_table) %in% train_index,]
        my_table_test <- my_table[row.names(my_table) %in% test_index,]
        resp_var_train <- metadata[row.names(metadata) %in% train_index,mta]
        resp_var_test <- metadata[row.names(metadata) %in% test_index,mta]
        resp_var_train <- droplevels(resp_var_train, exclude= NA, "")
        resp_var_test <- droplevels(resp_var_test, exclude= NA, "")
        print(levels(resp_var_test))
        print(length(resp_var_train))
        names(resp_var_test) <- row.names(my_table_test)
        if (length(levels(resp_var_test)) > 1 & length(resp_var_train) > 1){
          print("There is at least 2 groups and more than one sample in the resp var.")
          for (ntr in seq(opt$start_ntree, opt$end_ntree, by=opt$ntree_increment)){
            rf <- randomForest::randomForest(my_table_train, resp_var_train, ntree=ntr)
            print("made rf")
            pred <- predict(rf, my_table_test)
            # roc_data <- data.frame(pred = pred, resp_var_test = resp_var_test)
            score <- MLmetrics::Accuracy(pred, resp_var_test)
            print(paste("score:", score, "ntree:", ntr, "split:", counter))
            my_df <- rf$importance
            maxImp <- max(rf$importance)
            maxRow <- which(rf$importance == maxImp)
            
            if (counter == 1){
              wide_result_df$metadata[row_count] <- colnames(metadata)[mta]
              wide_result_df$dataset[row_count] <- transf_label
              wide_result_df$color[row_count] <- tabl[4]
              wide_result_df$ntrees[row_count] <- ntr
              wide_result_df[row_count, paste0("split",counter)] <- score
            }else{
              if (wide_result_df$metadata[row_count] == colnames(metadata)[mta] &&
                wide_result_df$dataset[row_count] == transf_label &&
                wide_result_df$color[row_count] == tabl[4] &&
                wide_result_df$ntrees[row_count] == ntr
              ){
                wide_result_df[row_count, paste0("split",counter)] <- score
              }else{
                print("There is a problem in the wide df logic.")
                this_is_going_to_exit()
              }
            }
            print(head(wide_result_df))
            row_count <- row_count + 1
            
            #Check its existence
            if (file.exists(main_output_fpath)) {
              print(paste0("Writing output to ", main_output_fpath, " ."))
              # main_header <- "all_score,	metadata_col,	rf_imp_se, rf_type, rf_ntree, trans_group, random_batch, cycle"
              cat(paste(paste0("\n", score), colnames(metadata)[mta], row.names(my_df)[maxRow], rf$type, #ilr_weight,	rf_imp_se, rf_type,
                        rf$ntree, transf_label, counter, ntr,#rf_ntree, trans_group, cycle
                        sep = ","),
                  file = main_output_fpath,
                  append=TRUE)
              }#for (ntr in 
            }#second if
          }#end if (length(levels(resp_var_test))...
        }#end for loop
  #     },#try catch
  #     error=function(cond) {
  #       print(paste("Opps, an error2 is thrown with", transf_label))
  #       message(paste(transf_label, cond))
  #     },
  #     warning=function(cond) {
  #       print(paste("Opps, a warning2 is thrown with", transf_label))
  #       message(paste(transf_label, cond))
  #     },
  #     finally=.Call(CfreadCleanup)
  #   )
  }#for mta
  print(paste("completed loop:", counter))
}

#Recording wide df----------------------------------------------------------
print("Recording wide df.")
# --------------------------------------------------------------------------
write.table(wide_result_df, row.names = FALSE, sep = ",",
            file = file.path(output_dir, "tables", paste0("wide_",main_output_label, ".csv")))


print("Reading in data from file.")
all_plot_data <- data.frame(read.table(file = main_output_fpath,
                                       sep = ",", header = TRUE))

print("finished R script")

