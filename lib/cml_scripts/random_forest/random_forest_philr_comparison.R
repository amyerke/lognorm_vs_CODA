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
    for (tax_w in philr_taxa_weights){
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
      for (tax_w in philr_taxa_weights){
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
                        help="Number of times to shuffle data and run loop again", 
                        metavar="character")
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

source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Set up constants-------------------------------------------------##
num_cycles <- opt$num_cycles
if(num_cycles < 3) stop("num_cycles should be 3 or more")
main_output_text <- "random_forest_score_R"
main_output_label <- paste0(main_output_text, num_cycles)
main_output_fn <- paste0(main_output_label, ".csv")
main_output_fpath <- file.path(output_dir, "tables", main_output_fn)
# philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
# philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")
philr_taxa_weights <- c("enorm")
philr_ilr_weights <- c("mean.descendants")
random_seed <- 36
main_header <- "all_score, metadata_col, rf_imp_se, rf_type, rf_ntree, trans_group, cycle"

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
print(paste("ncol(metadata)",ncol(metadata)))
# df[,-which(sapply(df, class) == "factor")]
rf_cols <- 1:ncol(metadata)#hack so I don't have to fix this in the function

tables <- list()
tables[[length(tables) + 1]] <- c("DADA2",file.path(output_dir, "tables", "ForwardReads_DADA2.txt"),"\t","r")
# tables[[length(tables) + 1]] <- c("HashSeq", file.path(output_dir,  "hashseq", "hashseq.csv"),",", "r")
tables[[length(tables) + 1]] <- c("lognorm_DADA2", file.path(output_dir, "tables", "lognorm_dada2.csv"), ",", "y")
# tables[[length(tables) + 1]] <- c("lognorm_HashSeq", file.path(output_dir,"tables", "lognorm_hashseq.csv"), ",", "y")
tables[[length(tables) + 1]] <- c("alr_DADA2", file.path(output_dir, "tables", "alr_asv.csv"), ",", "g")
# tables[[length(tables) + 1]] <- c("alr_HashSeq", file.path(output_dir,"tables", "alr_hashseq.csv"), ",", "g")
tables[[length(tables) + 1]] <- c("clr_DADA2", file.path(output_dir, "tables", "clr_asv.csv"), ",", "g")
# tables[[length(tables) + 1]] <- c("clr_HashSeq", file.path(output_dir,"tables", "clr_hashseq.csv"), ",", "m")
tables[[length(tables) + 1]] <- c("Silva_DADA2", file.path(output_dir,"tables", "Silva_DADA2", "Silva_DADA2.csv"), ",", "#64baeb")
tables <- add_PhILR_dfs_to_table(tables, file.path(output_dir, "tables", "Silva_DADA2"), "Silva_DADA2", color = "#c8e5f5")
tables <- add_random_tree_PhILRs_to_table(tables, file.path(output_dir, "tables", "Silva_DADA2"), "Silva_DADA2", color = "#1474aa", num_rand_trees=5)
tables[[length(tables) + 1]] <- c("Filtered_Silva_DADA2", file.path(output_dir,"tables", "Filtered_Silva_DADA2", "Filtered_Silva_DADA2.csv"), ",", "#f78646")
tables <- add_PhILR_dfs_to_table(tables, file.path(output_dir, "tables", "Filtered_Silva_DADA2"), "Filtered_Silva_DADA2", color = "#f5cbb3")
tables <- add_random_tree_PhILRs_to_table(tables, file.path(output_dir, "tables", "Filtered_Silva_DADA2"), "Filtered_Silva_DADA2", color = "#8c390b", num_rand_trees=5)
tables[[length(tables) + 1]] <- c("Filtered_UPGMA_DADA2", file.path(output_dir,"tables", "Filtered_UPGMA_DADA2", "Filtered_UPGMA_DADA2.csv"), ",", "#0000ff")
tables <- add_PhILR_dfs_to_table(tables, file.path(output_dir, "tables", "Filtered_UPGMA_DADA2"), "Filtered_UPGMA_DADA2", color = "#b7b7f3")
tables <- add_random_tree_PhILRs_to_table(tables, file.path(output_dir, "tables", "Filtered_UPGMA_DADA2"), "Filtered_UPGMA_DADA2", color = "#050598", num_rand_trees=5)
tables[[length(tables) + 1]] <- c("Filtered_IQtree", file.path(output_dir,"tables", "Filtered_IQtree", "Filtered_IQtree.csv"), ",", "orange")
tables <- add_PhILR_dfs_to_table(tables, file.path(output_dir, "tables", "Filtered_IQtree"), "Filtered_IQtree", color = "#f7d8a0")
tables <- add_random_tree_PhILRs_to_table(tables, file.path(output_dir, "tables", "Filtered_IQtree"), "Filtered_IQtree", color = "#e29302", num_rand_trees=5)

print(tables)

print(paste("Initilizing", main_output_fpath, "."))
print(paste("File header: ", main_header))
cat(main_header,
    file = main_output_fpath,
    append=FALSE)

counter <- 0

rowname_table <- data.frame(data.table::fread(file = tables[[1]][2],#this is a hack to make the train and test index work
                                              header=TRUE, data.table=FALSE),#the row names of this table should be
                            row.names = 1)#available in all of the other tables or there will be an error

print(paste("Counter:", counter, " entering main loop"))
for (counter in 1:num_cycles) {
  ##-Create training/testing sets-------------------------------------##
  train_index <- row.names(rowname_table)[sample(x = nrow(rowname_table), size = 0.75*nrow(rowname_table), replace=FALSE)]
  test_index <- row.names(rowname_table)[c(1:nrow(rowname_table))[!(1:nrow(rowname_table) %in% train_index)]]

  for(mta in rf_cols){
    print(paste("starting metadata col:", mta, colnames(metadata)[mta]))
    print(paste("Finding rows where", colnames(metadata)[mta], "has NA or Empty values from training and testing selectors."))
    na_or_empty_index <- which(is.na(metadata[,mta]) | metadata[,mta] == "")
    na_or_empty_rows <- row.names(metadata)[na_or_empty_index]
    print(length(train_index))
    train_index <- setdiff(train_index, na_or_empty_rows)
    # print(train_index)
    test_index <- setdiff(test_index, na_or_empty_rows)
    print(length(train_index))
    tryCatch({
      for (tabl in tables){
        print(tabl)
        transf_label <- tabl[1]
        my_table <- data.frame(data.table::fread(file = tabl[2],
                                                header=TRUE, data.table=FALSE),
                               row.names = 1)
        my_table_train <- my_table[row.names(my_table) %in% train_index,]
        my_table_test <- my_table[row.names(my_table) %in% test_index,]
        resp_var_train <- metadata[row.names(metadata) %in% train_index,mta]
        resp_var_test <- metadata[row.names(metadata) %in% test_index,mta]
        if (is.factor(resp_var_train)){
          resp_var_train <- droplevels(resp_var_train, exclude=c(NA, ""))
          resp_var_test <- droplevels(resp_var_test, exclude=c( NA, ""))
        }
        names(resp_var_test) <- row.names(my_table_test)
        if (length(resp_var_train) > 1){
          print("There is at least 2 groups and more than one sample in the resp var.")
          rf <- randomForest::randomForest(my_table_train, resp_var_train)
          print("made rf")
          pred <- predict(rf, my_table_test)
          roc_data <- data.frame(pred = pred, resp_var_test = resp_var_test)
          if (is.factor(resp_var_train)){
            # print(paste("nlevl resp_test:" length(levels(resp_var_test)))
            my_union <- base::union(levels(pred), levels(resp_var_test))
            levels(pred) <- my_union#hack for when the levels are different
            levels(resp_var_test) <- my_union
          }
          score <- MLmetrics::Accuracy(pred, resp_var_test)
          print(paste("score:", score))
          my_df <- rf$importance
          maxImp <- max(rf$importance)
          maxRow <- which(rf$importance == maxImp)
          
          #Check its existence
          if (file.exists(main_output_fpath)) {
            print(paste0("Writing output to ", main_output_fpath, " ."))
            # main_header <- "all_score,	metadata_col,	rf_imp_se, rf_type, rf_ntree, trans_group, random_batch, cycle"
            cat(paste(paste0("\n", score), colnames(metadata)[mta], row.names(my_df)[maxRow], rf$type, #ilr_weight,	rf_imp_se, rf_type,
                      rf$ntree, transf_label, counter, #rf_ntree, trans_group, cycle
                      sep = ","),
                file = main_output_fpath,
                append=TRUE)
            }#second if
          }#end if (length(levels(resp_var_test))...
        }#for loop
      },#try catch
      error=function(cond) {
        print(paste("Opps, an error2 is thrown with", transf_label))
        message(paste(transf_label, cond))
      },
      warning=function(cond) {
        print(paste("Opps, a warning2 is thrown with", transf_label))
        message(paste(transf_label, cond))
      }
      # finally=.Call(CfreadCleanup)
    )
  }#for mta
  print(paste("completed loop:", counter))
}


# --------------------------------------------------------------------------
print("Reshaping data from long to wide format")
# --------------------------------------------------------------------------
print("Reading in data from file.")
all_plot_data <- data.frame(read.table(file = main_output_fpath,
                                       sep = ",", header = TRUE))
print("Building wide df.")
data_cols <- unlist(lapply(1:num_cycles, function(x) paste0("split",x)))
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
            file = file.path(output_dir, "tables", paste0("wide_",main_output_label, ".csv")))

# --------------------------------------------------------------------------
print("Recording most important sequences.")
# --------------------------------------------------------------------------
unique(all_plot_data$rf_type)

best_seqs <- table(all_plot_data$rf_imp_seq)

write.table(best_seqs, row.names = FALSE,
            file = file.path(output_dir, "tables", paste0("best_seqs_",main_output_label, ".csv")))

weight_table <- data.frame(tree_type = c(F),
                           metadata = c(F),
                           taxa_pval = c(F),
                           ilr_pval = c(F))
weight_counter <- 1

print("Make all the boxplots")
##-Make all the boxplots--------------------------------------------##
pdf(file = file.path(output_dir, "graphics", paste0("boxplot_", main_output_label, ".pdf")))
for (mta in 1:length(unique(all_plot_data$metadata_col))){
  my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
  message(my_meta)
  #select plot data for each metadata cat
  plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]
  my_order <- lapply(tables, `[[`, 1)#takes first ele of each nested list
  my_colors  <- lapply(tables, `[[`, 4)
  plot_data$trans_group <- factor(plot_data$trans_group, levels = my_order)
  
  meta_mean <- mean(plot_data$all_score)
  g <- ggplot2::ggplot(plot_data, aes(y = all_score, x = trans_group, group=trans_group)) +
    ggplot2::geom_boxplot() +
    ggplot2::ggtitle(label = paste("R, 4-fold cv", project, my_meta)) +
    ggplot2::geom_hline(yintercept = meta_mean, color="red") +
    # scale_fill_discrete(labels=new_labels) +
    # ggplot2::ggtitle( label = paste("num_tg:", length(unique(new_pd$trans_group)))) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
    ggplot2::ylab("Score") +
    ggplot2::xlab("Transformation")
  print(g)
}

dev.off()
# 
# 
# print(paste("Making empty vectors to fill during plot building"))
# pval <- c()
# pw_name <- c()
# iw_name <- c()
# metadata_col <- c()
# transformation <- c()
# distance_metric <- c()
# pg_num <- c()
# 
# pdf(file = file.path(output_dir, "graphics", paste0("new_bp_", main_output_label, ".pdf")))
# index <- 1
# for (mta in 1:length(unique(all_plot_data$metadata_col))){
#   my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
#   message(my_meta)
#   
#   plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]
#   
#   not_uniform_tw <- which(plot_data$taxa_weight != "uniform" )
#   not_uniform_iw <- which(plot_data$ilr_weight != "uniform" )
#   philr_ds <- unique(plot_data$trans_group[c(not_uniform_tw,not_uniform_iw)] ) #pulls out philr only data
#   non_philr_ds <- unique(plot_data$trans_group[ !(plot_data$trans_group %in% philr_ds)])
#   non_philr_ds_pd <- subset(plot_data, trans_group %in% non_philr_ds)
#   philr_ds_pd <- data.frame(subset(plot_data, trans_group %in% philr_ds))
#   rownames(philr_ds_pd) <- seq(length=nrow(philr_ds_pd))
#   # my_means <- c()
#   # for (tg in 1:length(unique(new_pd$trans_group))){
#   #   trans_g <- new_pd$trans_group[tg]
#   #   my_vals <- which(new_pd$trans_group == trans_g)
#   #   my_means <- c(my_means, mean(new_pd$all_score[my_vals]))
#   #   names(my_means)[tg] <- trans_g
#   # }
#   new_pd <- rbind(non_philr_ds_pd, philr_ds_pd)
#   new_pd$trans_group <- factor(new_pd$trans_group, levels = c(non_philr_ds, philr_ds))
#   
#   for(tw in unique(plot_data$taxa_weight)){
#     for(iw in unique(plot_data$ilr_weight)){
#       philr_pd_tw_iw <- subset(philr_ds_pd, taxa_weight == tw & ilr_weight == iw)
#       jitter_pd <- rbind(non_philr_ds_pd, philr_pd_tw_iw)
#       jitter_pd$trans_group <- factor(jitter_pd$trans_group, levels = c(non_philr_ds, philr_ds))
#       #need to show means from new_pd, but show jitter of tw and iw
#       #or could just show selected points but show overal mean for each tw and iw
#       back_ground_points <- subset(philr_ds_pd, taxa_weight != tw & ilr_weight != iw)
#       bg_jitter <- rbind(non_philr_ds_pd, back_ground_points)
#       bg_jitter$trans_group <- factor(bg_jitter$trans_group, levels = c(non_philr_ds, philr_ds))
#       g <- ggplot2::ggplot(new_pd, aes(y = all_score, x = trans_group)) + 
#         ggplot2::geom_boxplot(data = jitter_pd, color = "blue", alpha = 0.5) +
#         ggplot2::geom_boxplot(data = bg_jitter, color = "red", alpha = 0.5) +
#         ggplot2::ggtitle(label = paste("Jones", my_meta, ", part_weight:", tw, ", ilr_weight:", iw)) +
#         # ggplot2::ggtitle( label = paste("num_tg:", length(unique(new_pd$trans_group)))) +
#         ggplot2::theme_classic() +
#         ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
#         ggplot2::ylab("AUC") +
#         ggplot2::xlab("Tree type")
#       print(g)
#       #build vectors for table
#       for (grp in unique(philr_pd_tw_iw$trans_group)){
#         # print(grp)
#         my_case <- philr_pd_tw_iw[philr_pd_tw_iw$trans_group == grp, ]$all_score
#         my_control <- back_ground_points[back_ground_points$trans_group == grp, ]$all_score
#         my_test <- t.test(my_case, my_control)
#         my_pval <- my_test$p.value
#         pval <- c(pval, my_pval)
#         pw_name <- c(pw_name, tw)
#         iw_name <- c(iw_name, iw)
#         metadata_col <- c(metadata_col, my_meta)
#         transformation <- c(transformation, grp)
#         pg_num <- c(pg_num, index)
#       }
#       index <- index + 1
#     }#end for iw
#   }#end for tw
# }
# 
# dev.off()
# 
# dFrame <- data.frame( pval, pw_name, iw_name, metadata_col, transformation, pg_num)
# dFrame$adj_pval <- p.adjust(dFrame$pval, method = "BH" )	
# dFrame <- dFrame [order(dFrame$adj_pval),]
# 
# write.table(dFrame, file=file.path(output_dir, "tables", paste0("new_bp_", main_output_label, ".csv")), 
#             sep=",", 
#             row.names=FALSE)
# 
# print(paste("completed"))
# 
# pdf(file = file.path(output_dir, "graphics", paste0("new_bp_NO_WEIGHT", main_output_label, ".pdf")))
# for (mta in 1:length(unique(all_plot_data$metadata_col))){
#   my_meta <- as.character(unique(all_plot_data$metadata_col)[mta])
#   message(my_meta)
#   
#   plot_data <- all_plot_data[all_plot_data$metadata_col == my_meta,]
#   
#   not_uniform_tw <- which(plot_data$taxa_weight != "uniform" )
#   not_uniform_iw <- which(plot_data$ilr_weight != "uniform" )
#   philr_ds <- unique(plot_data$trans_group[c(not_uniform_tw,not_uniform_iw)] ) #pulls out philr only data
#   non_philr_ds <- unique(plot_data$trans_group[ !(plot_data$trans_group %in% philr_ds)])
#   non_philr_ds_pd <- subset(plot_data, trans_group %in% non_philr_ds)
#   philr_ds_pd <- data.frame(subset(plot_data, trans_group %in% philr_ds))
#   rownames(philr_ds_pd) <- seq(length=nrow(philr_ds_pd))
# 
#   new_pd <- rbind(non_philr_ds_pd, philr_ds_pd)
#   new_pd$trans_group <- factor(new_pd$trans_group, levels = c(non_philr_ds, philr_ds))
#   
#   for(tw in unique(plot_data$taxa_weight)){
#     for(iw in unique(plot_data$ilr_weight)){
#       g <- ggplot2::ggplot(new_pd, aes(y = all_score, x = trans_group)) + 
#         ggplot2::geom_boxplot( color = "black",) +
#         ggplot2::ggtitle(label = paste(project, my_meta)) +
#         # ggplot2::ggtitle( label = paste("num_tg:", length(unique(new_pd$trans_group)))) +
#         ggplot2::theme_classic() +
#         ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
#         ggplot2::ylab("AUC") +
#         ggplot2::xlab("Tree type")
#       print(g)
#       #build vectors for table
#       for (grp in unique(philr_pd_tw_iw$trans_group)){
#         # print(grp)
#         my_case <- philr_pd_tw_iw[philr_pd_tw_iw$trans_group == grp, ]$all_score
#         my_control <- back_ground_points[back_ground_points$trans_group == grp, ]$all_score
#         my_test <- t.test(my_case, my_control)
#         my_pval <- my_test$p.value
#         pval <- c(pval, my_pval)
#         pw_name <- c(pw_name, tw)
#         iw_name <- c(iw_name, iw)
#         metadata_col <- c(metadata_col, my_meta)
#         transformation <- c(transformation, grp)
#         pg_num <- c(pg_num, index)
#       }
#       index <- index + 1
#     }#end for iw
#   }#end for tw
# }
# 
# dev.off()

print("finished R script")



