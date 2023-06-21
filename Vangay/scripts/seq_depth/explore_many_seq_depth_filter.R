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
if (!requireNamespace("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
library("compositions")
library("phyloseq")
library("vegan")
library("DESeq2")
library("philr")
library("ape")
library("ALDEx2")
library("ggplot2")
library("edgeR")

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

#Plot shannon diversity against log10(total_seqs)
#Plot other normalization methods: otu log 100 and alr and clr 
my_ds_names <- c( "raw seqs", "clr(raw seqs)", "lognorm raw seqs", "philr ref", "DESeq2", "ALDEx2.clr")
min_seq_depths <- c(0, 500, 1000, 5000, 10000, 20000, 40000)
mds_depth <- 5

total_seqs <- rowSums(asv_table)
total_seqs <- data.frame(total_seqs, row.names = row.names(asv_table))

kend <- vector(mode = "numeric", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
perma_r2 <- vector(mode = "numeric", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
ds_num <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
ds_nam <- vector(mode = "character", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
mds_lev <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
seq_depth <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
var_exp <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
spear_cor <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
samples_left <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
taxa_left <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)
zero_count <- vector(mode = "integer", length = length(my_ds_names) * length(min_seq_depths) * mds_depth)


pdf(file = file.path(output_dir, "graphics", "seq_depth_artifact_PCA12345_scatter.pdf"))
counter <- 1
for(s in 1:length(min_seq_depths)){
  seq_d <- min_seq_depths[s]#new sequencing depth
  sd_filt_asv <- asv_table[total_seqs$total_seqs >= seq_d,]#dataset 1
  print(paste("sd_filtered dim:", paste(dim(sd_filt_asv))))
  safe_rns <- intersect(row.names(ref_ps@otu_table), row.names(sd_filt_asv)) #rows for this iterate
  ts <- rowSums(sd_filt_asv[safe_rns,])
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
  
  ref_philr <- philr::philr(new_ref_ps@otu_table, new_ref_ps@phy_tree,
                            part.weights='enorm.x.gm.counts',
                            ilr.weights='blw.sqrt')#dataset 4
  print(paste("made new philr", dim(as.data.frame(ref_philr))))

  ln_asv <- lognorm(sd_filt_asv)#dataset 6
  
  ald <- ALDEx2::aldex.clr(sd_filt_asv, mc.samples=12, denom="all", verbose=F)
  ald <- data.frame(ald@analysisData)
  print(paste("size of ald:", object.size(ald)))
  print(paste("ald dim:", paste(dim(ald))))
  
  my_datasets <- list(sd_filt_asv, my_clr,  ln_asv, ref_philr, new_DESeq2, ald)
  
  print(paste("finished seq depth filter:", s))
  
  for( ds in 1:length(my_datasets)){
    print(my_ds_names[ds])
    my_table <- as.data.frame(my_datasets[ds])
    zeros <- sum(my_table == 0)
    print(dim(my_table))
    ##-Create a PCA-----------------------------------------------------##
    my_prcmp <- prcomp(my_table, 
                       center = TRUE,
                       rank = mds_depth)#,
    # scale = TRUE)
    ##-Extract PCA matrix and convert to dataframe----------------------##
    myPCA <- data.frame(my_prcmp$x)
    my_var_exp <- my_prcmp$sdev^2/sum(my_prcmp$sdev^2)
    print(paste("finished ds:", ds))
    for (md in 1:mds_depth){
      # kend[counter] <- cor.test(log10(total_seqs[total_seqs > seq_d]), myPCA[,md], method = "kendall")$estimate
      # perma_r2[counter] <- adonis2(log10(total_seqs[total_seqs > seq_d]) ~ myPCA[,md])$R2[1]
      if (md == 1){
        my_spear <- cor.test(ts, myPCA[,md],
                             method = "spearman")
        plot(log10(ts), myPCA[,2],
             main = paste0( project, my_datasets[ds], ", PCA", md,"\nr_sq: ", my_spear$estimate^2),
             sub = paste0("spear: ", my_spear$estimate),
             xlab = paste0("Sequencing depth"),
             ylab = paste0("PCA", md)
            )
      }
      ds_num[counter] <- ds
      ds_nam[counter] <- my_ds_names[ds]
      mds_lev[counter] <- md
      seq_depth[counter] <- seq_d
      var_exp[counter] <- my_var_exp[md]
      spear_cor[counter] <- cor(total_seqs[total_seqs >= seq_d], myPCA[,md], method = "spearman")
      samples_left[counter] <- nrow(my_table)
      taxa_left[counter] <- ncol(my_table)
      zero_count[counter] <- zeros
      counter <- counter + 1
    }
  }
}
dev.off()
result_df <- data.frame(ds_num, ds_nam, perma_r2, mds_lev, seq_depth, var_exp, spear_cor, samples_left, zero_count, taxa_left)
print("created resulting DF")


write.table(result_df, 
            file = file.path(output_dir, "tables", paste0(project, "_PCA_seqdep_filt_results.csv")),
            sep = ",")

pdf(file = file.path(output_dir, "graphics", "seq_depth_artifact_PCA12345_bar.pdf"))
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only,
                       aes(x=as.factor(seq_depth), y=spear_cor^2, fill=ds_nam)) +
    geom_bar(width = 0.8, position=position_dodge(width = 1), stat="identity",) +
    ggtitle(paste0(project, ": PCA",  i, " vs total sequences per sample")) +
    xlab("Min sequence depth per sample") +
    ylab("Spearman Rsq")
  print(g)
}
dev.off()
print("made bar charts")

pdf(file = file.path(output_dir, "graphics", "seq_depth_artifact_PCA12345_line.pdf"))
for (i in 1:max(result_df$mds_lev)){
  pca_only <- result_df[mds_lev == i, ]
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=seq_depth, y=spear_cor^2, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::annotate("text", x = head(pca_only$seq_depth, n = length(my_ds_names)), 
                      y = head(c(pca_only$spear_cor^2), n = length(my_ds_names)), 
                      label = head(pca_only$ds_nam, n = length(my_ds_names)),
                      hjust = -0.1) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total sequences per sample")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("paste(italic(R) ^ 2, \" = .75\")") + 
    ggplot2::labs(fill = "Transformations") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  print(g)
  
  #plot remaining taxa
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=seq_depth, y=taxa_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total sequences per sample")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("Taxa") + 
    ggplot2::labs(fill = "Transformations") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
  
  #plot remaining samples
  g <- ggplot2::ggplot(pca_only, 
                       aes(x=seq_depth, y=samples_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total sequences per sample")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("Samples") + 
    ggplot2::labs(fill = "Transformations") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
  
  #plot percentage of zeros
  g <- ggplot2::ggplot(pca_only,
                       aes(x=seq_depth, y=zero_count/samples_left*taxa_left, group = ds_nam)) +
    ggplot2::geom_point(aes(color = factor(ds_nam))) +
    ggplot2::geom_line(aes(color = factor(ds_nam))) +
    ggplot2::ggtitle(paste0(project, ": PCA",  i, " vs total sequences per sample")) +
    ggplot2::xlab("Min sequence depth per sample") +
    ggplot2::ylab("Percentage of Zeros") +
    ggplot2::labs(fill = "Transformations") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::theme_minimal()
  g
  print(g)
}
dev.off()
print("made line chart")

