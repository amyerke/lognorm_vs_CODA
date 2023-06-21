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

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
# library("ALDEx2")
# if (!requireNamespace("optparse", quietly = TRUE)){ install.packages("optparse") }
# library("optparse")
# if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
# library("rgr")
# if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
# library("phyloseq")
# if (!requireNamespace("vegan", quietly = TRUE)) BiocManager::install("vegan")
# library("vegan")
# if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
# library("DESeq2")
# if (!requireNamespace("philr", quietly = TRUE)) BiocManager::install("philr")
# library("philr")
# if (!requireNamespace("ape", quietly = TRUE)) BiocManager::install("ape")
# library("ape")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("ggplot2")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Zeller"
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

my_ds_names <- c( "raw counts table", "clr", "alr", "lognorm", "PhILR Silva Tree", "DESeq2", "ALDEx2.clr")
min_seq_depths <- c(0, 500, 1000, 5000, 10000, 20000, 30000, 50000, 70000, 90000, 10000)
mds_depth <- 5

result_df <- read.table(file = file.path(output_dir, "tables", paste0(project, "_PCA_seqdep_filt_results.csv")),
                        sep = ",")

pdf(file = file.path(output_dir, "graphics", "seq_depth_artifact_PCA12345_bar.pdf"))
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only,
                       aes(x=as.factor(seq_depth), y=spear_cor^2, fill=ds_nam)) +
    ggplot2::geom_bar(width = 0.8, position=position_dodge(width = 1), stat="identity",) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total reads per sample")) +
    ggplot2::xlab("Min read depth per sample") +
    ggplot2::ylab("Spearman Rsq")
  print(g)
}
dev.off()
print("made bar charts")

pdf(file = file.path(output_dir, "graphics", "seq_depth_artifact_PCA12345_line.pdf"))
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[result_df$mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=seq_depth, y=spear_cor^2, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::annotate("text", x = head(pca_only$seq_depth, n = length(my_ds_names)), 
                      y = head(c(pca_only$spear_cor^2), n = length(my_ds_names)), 
                      label = head(pca_only$ds_nam, n = length(my_ds_names)),
                      hjust = -0.1) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total reads per sample")) +
    ggplot2::xlab("Min read depth per sample (log scale)") +
    ggplot2::ylab(paste0("R squared: PC", i, " vs read depth")) + 
    ggplot2::labs(fill = "Transformations") +
    ggplot2::scale_x_continuous(trans='log10',
                                breaks=trans_breaks('log10', function(x) 10^x),
                                labels=trans_format('log10', math_format(10^.x))) +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text=element_text(size=20), #change font size of all text
              axis.text=element_text(size=20), #change font size of axis text
              axis.title=element_text(size=20), #change font size of axis titles
              plot.title=element_text(size=20), #change font size of plot title
              legend.text=element_text(size=20), #change font size of legend text
              legend.title=element_text(size=20)) #change font size of legend title
  
  print(g)
  
  #plot remaining taxa
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=seq_depth, y=taxa_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total reads per sample")) +
    ggplot2::xlab("Min read depth per sample") +
    ggplot2::ylab("Taxa") + 
    ggplot2::labs(fill = "Transformations") +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
  
  #plot remaining samples
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=seq_depth, y=samples_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total reads per sample")) +
    ggplot2::xlab("Min read depth per sample") +
    ggplot2::ylab("Samples") + 
    ggplot2::labs(fill = "Transformations") +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
  
  #plot percentage of zeros
  g <- ggplot2::ggplot(pca_only,
                       aes(x=seq_depth, y=zero_count/samples_left*taxa_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total reads per sample")) +
    ggplot2::xlab("Min read depth per sample") +
    ggplot2::ylab("Percentage of Zeros") +
    ggplot2::labs(fill = "Transformations") +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
}
dev.off()
print("made line chart")