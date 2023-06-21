# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for writing phylo trees in phyloxml format
# This was helepful:https://uscbiostats.github.io/rphyloxml/

rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##
library("phyloseq")
library("ape")
BiocManager::install("rphyloxml")
devtools::install_github("USCBiostats/rphyloxml")
library(rphyloxml)

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Jones"
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Import tables and data preprocessing-----------------------------##
ref_ps <- readRDS(file.path(output_dir, "r_objects", "ref_tree_phyloseq_obj.rds"))

#view tree

plot_tree(ref_ps, label.tips="Genus", color="type", sizebase = 10)

plot_tree(ref_ps, label.tips="Genus", color="type", sizebase = 0.5, text.size = 1.5)
plot_tree(ref_ps, label.tips="Genus",shape = ".", sizebase = 0.5, text.size = 1.5)

# plot_tree(prune_tree(ref_ps))




##-Create phyloxml text---------------------------------------------##
z <- write_phyloxml(phy_tree(ref_ps))

##-write phyloxml text to file--------------------------------------##
xml2::write_xml(z, 
                file.path(output_dir, "trees", paste0(project, "_ref_tree.xml")))

write.table(ref_ps@tax_table@.Data,
            file.path(output_dir, "taxonomy", paste0(project, "_ref_tree_tax.csv")),
            sep = ",")

write.tree(phy_tree(ref_ps), 
file = file.path(output_dir, "trees", paste0(project, "_ref_tree.newick")))




