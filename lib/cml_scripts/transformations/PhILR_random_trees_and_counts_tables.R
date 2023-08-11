# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making all weightschemes in PhILR 

rm(list = ls()) #clear workspace

####-functions--------------------------------------------------------####
raw_ps_to_clean_ps <- function(ps) {
  #requires ape, phyloseq and philr_tutorial_normalization 
  clean_otu = data.frame(ps@otu_table)
  clean_otu = philr_tutorial_normalization(clean_otu)
  ps_clean = phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                 phy_tree(ps@phy_tree),
                                 tax_table(ps@tax_table), 
                                 sample_data(ps@sam_data))
  return(ps_clean)
}
make_PhILR_transform_tables <- function(counts_table,
                                        tree,
                                        table_name,
                                        output_folder,
                                        philr_taxa_weights = c("uniform","gm.counts","enorm","enorm.x.gm.counts"),
                                        # philr_taxa_weights = c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts"),
                                        philr_ilr_weights = c("uniform","blw","blw.sqrt","mean.descendants"),
                                        save_counts_table = FALSE
                                        ){
  print(paste0("Making all possible PhILR Weighted transformations of ", table_name, "." ))
  print(paste0("Using these taxa weights ", paste(philr_taxa_weights, collapse = " "), "." ))
  print(paste0("And these ilr weights ", paste(philr_ilr_weights, collapse = " "), "." ))
  if(!file.exists(output_folder)){
    print(paste("Creating folder at ", output_folder))
    dir.create(output_folder)
  }
  if (save_counts_table == TRUE){
    print("Saving counts table.")
    write.table(counts_table, 
              file = file.path(output_folder, paste0(table_name,".csv")),
              sep = ",", row.names = TRUE)
  }
  tree = ape::makeNodeLabel(tree, method="number", prefix='n')
	if (any(counts_table == 0)){
				    print("adding pseudocount of 1 before PhILR transform")
				    counts_table = counts_table + 1
				  }
  print(paste0("Building output files of ", table_name, "." ))
  for (ilr_w in 1:length(philr_ilr_weights)){
    iw <- philr_ilr_weights[ilr_w]
    for (tax_w in 1:length(philr_taxa_weights)){
      pw <- philr_taxa_weights[tax_w]
      table_name_full <- paste0(paste(table_name, iw, pw, sep = "_"),".csv")
      if(!file.exists(file.path(output_folder, table_name_full))){
        if(any(tree$node.label == "")){
          my_nodes <- paste0("n", 1:tree$Nnode)
          tree$node.label <- my_nodes
        }
        print(philr_taxa_weights[tax_w])
        print(philr_ilr_weights[ilr_w])
        my_table <- philr::philr(counts_table,
                                 tree,
                                 part.weights = philr_taxa_weights[tax_w],
                                 ilr.weights = philr_ilr_weights[ilr_w])
        print(paste0("Saving ", table_name_full, " to ", output_folder))
        write.table(my_table, 
                  file = file.path(output_folder, table_name_full),
                  sep = ",", row.names = TRUE)
      }else{
        print(paste(table_name_full, "already exists at", output_folder))
      }
    }
  }
}
make_random_tree_philrs <- function(counts_table,
                                    tree,
                                    table_name,
                                    output_folder,#combines all output into same folder
                                    num_random_trees,
                                    # philr_taxa_weights = c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts"),
                                    philr_taxa_weights = c("uniform","gm.counts","enorm","enorm.x.gm.counts"),
                                    philr_ilr_weights = c("uniform","blw","blw.sqrt","mean.descendants")
                                    ){
  random_seed <- 36
  print(paste("Setting random seed to:", random_seed))
  set.seed(random_seed)
  print("Making random trees")
  for (rand in 1:num_random_trees){
    print(paste("random tree", rand))
    my_rand_tree <- ape::rtree(n = length(tree$tip.label), tip.label = tree$tip.label)
    my_rand_tree <- ape::makeNodeLabel(my_rand_tree, method="number", prefix='n')
    table_name_full <- paste0("Shuffle", rand, "_PhILR_", table_name)
    make_PhILR_transform_tables(counts_table = counts_table,
                                tree = my_rand_tree,
                                table_name = table_name_full,
                                output_folder = output_folder,
                                philr_taxa_weights = philr_taxa_weights,
                                philr_ilr_weights = philr_ilr_weights)
    print(paste("Ending make_random_tree_philrs of", table_name))
  }
}
####-Load Depencencies------------------------------------------------####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
library("ape")
if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
library("philr")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
library("phyloseq")
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

####-reading cml input------------------------------------------------####
print("Finished loading libraries, now reading cml input.")


option_list <- list(
  make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git',"lognorm_vs_CODA"), 
              help="dataset dir path", metavar="home dir"),
  make_option(c("-p", "--project"), type="character", default="Jones", 
              help="project folder", metavar="project")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)
####-Establish directory constants------------------------------------####
print("Establishing directory constants.")

home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

source(file.path(home_dir, "lib", "table_manipulations.R"))

####-Set up other constants-------------------------------------------####
print("Setting up other constants.")

main_output_label <- paste0("philr_weights_shuffled")
# philr_taxa_weights <- c("uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts")
philr_taxa_weights <- c("uniform","gm.counts", "enorm","enorm.x.gm.counts")
philr_ilr_weights <- c("uniform","blw","blw.sqrt","mean.descendants")

####-Import and prepping Silva_DADA2----------------------------------####
print("Transforming non-filtered phyloseq objects first.")
print("Importing and prepping Silva_DADA2.")

# pdf(file = file.path(output_dir, "graphics", paste0("trees_", main_output_label, ".pdf")))
# phylo_obj <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
# phylo_label <- "Silva_DADA2"
# print(paste0("Counts table dimensions of ", phylo_label, ": ", dim(phylo_obj@otu_table), collapse = ""))
# make_PhILR_transform_tables(phylo_obj@otu_table,
#                             phylo_obj@phy_tree,
#                             phylo_label,
#                             file.path(output_dir, "tables", phylo_label),
#                             save_counts_table = TRUE)
# make_random_tree_philrs(phylo_obj@otu_table,
#                         phylo_obj@phy_tree,
#                         phylo_label,
#                         file.path(output_dir, "tables", phylo_label),
#                         3)
# phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_ref"))
# 
# phy_tree(phylo_obj) <- ape::makeNodeLabel(phy_tree(phylo_obj), method="number", prefix='n')
# print("Cleaning Silva tree otu with philr tutorial normalization")
# phylo_obj <- raw_ps_to_clean_ps(phylo_obj)
# phylo_label <- "Filtered_Silva_DADA2"
# phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_ref"))
# print(paste0("Counts table dimensions of ", phylo_label, ": ", dim(phylo_obj@otu_table), collapse = ""))
# make_PhILR_transform_tables(phylo_obj@otu_table,
#                             phylo_obj@phy_tree,
#                             phylo_label,
#                             file.path(output_dir, "tables", phylo_label),
#                             save_counts_table = TRUE)
# make_random_tree_philrs(phylo_obj@otu_table,
#                         phylo_obj@phy_tree,
#                         phylo_label,
#                         file.path(output_dir, "tables", phylo_label),
#                         3)
# 
# ####-Import and prepping UPGMA_DADA2----------------------------------####
# print("Importing UPGMA phyloseq")
# 
# phylo_obj <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_UPGMA_phyloseq_obj.rds"))
# phy_tree(phylo_obj) <- ape::makeNodeLabel(phy_tree(phylo_obj), method="number", prefix='n')
# phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_upgma"))
# print("Cleaning UPGMA with philr tutorial normalization")
# phylo_obj <- raw_ps_to_clean_ps(phylo_obj)
# phylo_obj <- phyloseq::transform_sample_counts(phylo_obj, function(x) x + 1 )
# phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_upgma"))
# phylo_label <- "Filtered_UPGMA_DADA2"
# print(paste0("Counts table dimensions of ", phylo_label, ": ", dim(phylo_obj@otu_table), collapse = ""))
# make_PhILR_transform_tables(phylo_obj@otu_table,
#                             phylo_obj@phy_tree,
#                             phylo_label,
#                             file.path(output_dir, "tables", phylo_label),
#                             save_counts_table = TRUE)
# make_random_tree_philrs(phylo_obj@otu_table,
#                         phylo_obj@phy_tree,
#                         phylo_label,
#                         file.path(output_dir, "tables", phylo_label),
#                         3)
# 
# ####-Import and prepping IQTREE2_DADA2--------------------------------####
# print("Importing IQTree phyloseq")
# 
# phylo_obj <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_iqtree_phyloseq_obj.rds"))
# phy_tree(phylo_obj) <- ape::makeNodeLabel(phy_tree(phylo_obj), method="number", prefix='n')
# phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_iqtree"))
# print("Cleaning IQ-tree with philr tutorial normalization")
# phylo_obj <- raw_ps_to_clean_ps(phylo_obj)
# phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_iqtree"))
# phylo_label <- "Filtered_IQtree"
# print(paste0("Counts table dimensions of ", phylo_label, ": ", dim(phylo_obj@otu_table), collapse = ""))
# make_PhILR_transform_tables(phylo_obj@otu_table,
#                             phylo_obj@phy_tree,
#                             phylo_label,
#                             file.path(output_dir, "tables", phylo_label),
#                             save_counts_table = TRUE)
# make_random_tree_philrs(phylo_obj@otu_table,
#                         phylo_obj@phy_tree,
#                         phylo_label,
#                         file.path(output_dir, "tables", phylo_label),
#                         3)
# 
# dev.off()

pdf(file = file.path(output_dir, "graphics", paste0("trees_prev_filtered_90", main_output_label, ".pdf")))
####-Import and prepping prevelence filtered Silva_DA-----------------####
print("Transforming 90% prevalence filtered phyloseq objects.")
print("Importing and prepping prevalence filtered filtered Silva_DADA2.")

phylo_obj <- readRDS(file.path(output_dir, "r_objects", "ltp_90prcnt_filt.rds"))
phylo_label <- "prev_filt90_Silva_DADA2"
print(paste0("Counts table dimensions of ", phylo_label, ": ", dim(phylo_obj@otu_table), collapse = ""))
make_PhILR_transform_tables(phylo_obj@otu_table,
                            phylo_obj@phy_tree,
                            phylo_label, 
                            file.path(output_dir, "tables", phylo_label),
                            save_counts_table = TRUE)
make_random_tree_philrs(phylo_obj@otu_table,
                        phylo_obj@phy_tree,
                        phylo_label, 
                        file.path(output_dir, "tables", phylo_label),
                        3)
phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_ref_90prcnt_filt"))

# phy_tree(phylo_obj) <- ape::makeNodeLabel(phy_tree(phylo_obj), method="number", prefix='n')
# print("Cleaning prevalence filtered Silva tree otu with philr tutorial normalization")
# phylo_obj <- raw_ps_to_clean_ps(phylo_obj)
# phylo_label <- "prev_filt90_Filtered_Silva_DADA2"
# phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_ref_90prcnt_filt"))
# print(paste0("Counts table dimensions of ", phylo_label, ": ", dim(phylo_obj@otu_table), collapse = ""))
# make_PhILR_transform_tables(phylo_obj@otu_table,
#                             phylo_obj@phy_tree,
#                             phylo_label, 
#                             file.path(output_dir, "tables", phylo_label),
#                             save_counts_table = TRUE)
# make_random_tree_philrs(phylo_obj@otu_table,
#                         phylo_obj@phy_tree,
#                         phylo_label, 
#                         file.path(output_dir, "tables", phylo_label),
#                         3)

####-Import and prepping UPGMA_DADA2----------------------------------####
print("Importing UPGMA phyloseq")

phylo_obj <- readRDS(file.path(output_dir, "r_objects", "filtered_90prcntdenovo_tree_UPGMA_phyloseq_obj.rds"))
phy_tree(phylo_obj) <- ape::makeNodeLabel(phy_tree(phylo_obj), method="number", prefix='n')
phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_upgma_90prcnt_filt"))
print("Cleaning prevalence filtered UPGMA with philr tutorial normalization")
phylo_obj <- raw_ps_to_clean_ps(phylo_obj)
phylo_obj <- phyloseq::transform_sample_counts(phylo_obj, function(x) x + 1 )
phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_upgma_90prcnt_filt"))
phylo_label <- "prev_filt90_Filtered_UPGMA_DADA2"
print(paste0("Counts table dimensions of ", phylo_label, ": ", dim(phylo_obj@otu_table), collapse = ""))
make_PhILR_transform_tables(phylo_obj@otu_table,
                            phylo_obj@phy_tree,
                            phylo_label, 
                            file.path(output_dir, "tables", phylo_label),
                            save_counts_table = TRUE)
make_random_tree_philrs(phylo_obj@otu_table,
                        phylo_obj@phy_tree,
                        phylo_label, 
                        file.path(output_dir, "tables", phylo_label),
                        3)
####-Import and prepping UPGMA_DADA2----------------------------------####
print("Importing prevalence filtered IQTree phyloseq")

phylo_obj <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_iqtree_phyloseq_obj.rds"))
phy_tree(phylo_obj) <- ape::makeNodeLabel(phy_tree(phylo_obj), method="number", prefix='n')
phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("orig_iqtree_90prcnt_filt"))
print("Cleaning IQ-tree with philr tutorial normalization")
phylo_obj <- raw_ps_to_clean_ps(phylo_obj)
phyloseq::plot_tree(phylo_obj, method = "treeonly", nodelabf=nodeplotblank, title = paste0("cln_iqtree_90prcnt_filt"))
phylo_label <- "prev_filt90_Filtered_IQtree"
print(paste0("Counts table dimensions of ", phylo_label, ": ", dim(phylo_obj@otu_table), collapse = ""))
make_PhILR_transform_tables(phylo_obj@otu_table,
                            phylo_obj@phy_tree,
                            phylo_label, 
                            file.path(output_dir, "tables", phylo_label),
                            save_counts_table = TRUE)
make_random_tree_philrs(phylo_obj@otu_table,
                        phylo_obj@phy_tree,
                        phylo_label, 
                        file.path(output_dir, "tables", phylo_label),
                        3)

dev.off()


print("script complete")

