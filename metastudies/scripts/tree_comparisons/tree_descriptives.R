# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making distribution of AUC of ILR balances to compare to philr

rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
if (!requireNamespace("randomForest", quietly = TRUE)) BiocManager::install("randomForest")
if (!requireNamespace("ROCR", quietly = TRUE)) BiocManager::install("ROCR")
if (!requireNamespace("ggpubr", quietly = TRUE)) BiocManager::install("ggpubr")
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
library("rgr")
library("phyloseq")
library("ggpubr")
library("ROCR")
library("philr")
library("ggplot2")
library("randomForest")
library("ape")
library("compositions")
print("finished loading libraries")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
projects <- c("Jones", "Vangay", "Zeller", "Noguera-Julian")

##-Functions--------------------------------------------------------##
raw_ps_to_clean_ps <- function(ps) {
  #requires ape, phyloseq and philr_tutorial_normalization 
  clean_otu = data.frame(ps@otu_table@.Data)
  clean_otu = philr_tutorial_normalization(clean_otu)
  ps_clean = phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                 phy_tree(ps@phy_tree),
                                 tax_table(ps@tax_table), 
                                 sample_data(ps@sam_data))
  return(ps_clean)
}


source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Create tree attribute vectors------------------------------------##
project_name <-c()
tree_name <- c()
num_nodes <- c()
num_tips <- c()
ave_branch_length <- c()
is_ultrametric <- c()
var_branch_length <- c()
philr_ncol <- c()
# compute.brlen(phy, method = "Grafen", power = 1, ...)
# base.freq()

for (project in projects) {
  ##-Set up constants-------------------------------------------------##
  print("making empty lists for trees and data.frames")
  phylo_objects <- list()
  phylo_objects_names <- c()
  output_dir <- file.path(home_dir, project, 'output')
  
  ##-Import tables and data preprocessing-----------------------------##
  asv_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))
  
  print("Cleaning Ref tree otu with philr tutorial normalization")
  ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))
  phylo_objects <- append(phylo_objects, ref_ps)
  phylo_objects_names <- append(phylo_objects_names, "Silva")
  
  
  #TODO: Fix this spaghetti code
  clean_otu <- data.frame(ref_ps@otu_table@.Data)
  clean_otu <- philr_tutorial_normalization(clean_otu)
  print(paste("nrow orginal ref:", nrow(ref_ps@otu_table), "nrow clean ref: ", nrow(clean_otu)))
  
  phy_tree(ref_ps) <- ape::makeNodeLabel(phy_tree(ref_ps), method="number", prefix='n')
  ref_ps_clean <- phyloseq::phyloseq( otu_table(clean_otu, taxa_are_rows = F), 
                                      phy_tree(ref_ps@phy_tree),
                                      tax_table(ref_ps@tax_table), 
                                      sample_data(ref_ps@sam_data))
  phylo_objects <- append(phylo_objects, ref_ps_clean)
  phylo_objects_names <- append(phylo_objects_names, "Filtered_Silva")
  
  print("Cleaning UPGMA tree otu with philr tutorial normalization")
  denovo_tree_ps <- readRDS(file.path(output_dir, "r_objects", "denovo_tree_UPGMA_phyloseq_obj.rds"))
  clean_den_otu <- philr_tutorial_normalization(data.frame(denovo_tree_ps@otu_table@.Data))
  print(paste("nrow orginal denovo:", nrow(denovo_tree_ps@otu_table), "nrow clean denovo otu: ", nrow(clean_den_otu)))
  cln_denovo_tree_ps <- phyloseq::phyloseq( otu_table(clean_den_otu, taxa_are_rows = F),
                                            phy_tree(ape::makeNodeLabel(phy_tree(denovo_tree_ps@phy_tree))),
                                            tax_table(denovo_tree_ps@tax_table), 
                                            sample_data(denovo_tree_ps@sam_data))
  denovo_tree_ps <- transform_sample_counts(denovo_tree_ps, function(x) x + 1 )
  phy_tree(ref_ps_clean) <- makeNodeLabel(phy_tree(ref_ps_clean), method="number", prefix='n')
  phy_tree(cln_denovo_tree_ps) <- makeNodeLabel(phy_tree(cln_denovo_tree_ps), method="number", prefix='n')
  phylo_objects <- append(phylo_objects, cln_denovo_tree_ps)
  phylo_objects_names <- append(phylo_objects_names, "Filtered_UPGMA")
  phylo_objects <- append(phylo_objects, denovo_tree_ps)
  phylo_objects_names <- append(phylo_objects_names, "UPGMA")
  
  # Adding iqtree to list----------------------------------------------#
  print("Adding iqtree to list")
  # Adding iqtree to list----------------------------------------------#
  iqtree <- readRDS(file.path(output_dir, "r_objects","denovo_tree_iqtree_phyloseq_obj.rds"))
  phylo_objects <- append(phylo_objects, iqtree)
  phylo_objects_names <- append(phylo_objects_names, "IQTREE")
  phylo_obj <- raw_ps_to_clean_ps(iqtree)
  phylo_objects <- append(phylo_objects, phylo_obj)
  phylo_objects_names <- append(phylo_objects_names, "Filtered_IQTREE")
  
  print("updating tree and table descriptives from phylo objects")
  for (phy in 1:length(phylo_objects)){
    print(phy)
    #tree
    my_phy <- phylo_objects[[phy]]
    my_tree <- my_phy@phy_tree
    project_name <- c(project_name, project)
    tree_name <- c(tree_name, phylo_objects_names[phy])
    ave_branch_length <- c(ave_branch_length, mean(my_tree$edge.length))
    var_branch_length <- c(var_branch_length, var(my_tree$edge.length))
    num_nodes <- c(num_nodes, my_tree$Nnode)
    num_tips <- c(num_tips, length(my_tree$tip.label))
    is_ultrametric <- c(is_ultrametric, ape::is.ultrametric(my_tree))
  }
  
}

print("forming dataframes from tree data")
my_trees <- data.frame(project_name,
                       tree_name,
                       num_nodes,
                       num_tips,
                       ave_branch_length,
                       var_branch_length,
                       is_ultrametric)

write.table(my_trees,
            sep = ",",
            row.names = FALSE,
            file = file.path(home_dir,"metastudies","output", "tree_descriptives.csv"))

