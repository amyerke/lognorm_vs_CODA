# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making some plots for a presentation about coda data transformations

rm(list = ls()) #clear workspace

####-Functions--------------------------------------------------------####
print("Loading custom functions.")

####-Load Depencencies------------------------------------------------####
options(repos = list(CRAN="http://cran.rstudio.com/"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("compositions", quietly = TRUE)) BiocManager::install("compositions")
library("compositions")
print("finished loading libraries")

####-Establish directory layout---------------------------------------####
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Zeller"
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))
source(file.path(home_dir, "lib", "table_manipulations.R"))

initial_table <- data.frame(readRDS(file.path(output_dir, "r_objects", "ForwardReads_DADA2.rds")))

filt_it <- initial_table[,sapply(initial_table, function(x) mean(x == 0) < 0.9)]

sample_2D <- filt_it[,1:2]

plot(sample_2D[,1], sample_2D[,2],
     main = "Raw counts",
     xlab = paste0(substring(colnames(sample_2D)[1], 1,10),"..."),
     ylab = paste0(substring(colnames(sample_2D)[2], 1,10),"..."))

sample_2D_prop <- proportions_transform(sample_2D)

plot(sample_2D_prop[,1], sample_2D_prop[,2],
     main = "Relative abundances",
     xlab = paste0(substring(colnames(sample_2D)[1], 1,10),"..."),
     ylab = paste0(substring(colnames(sample_2D)[2], 1,10),"..."))

clr_sample_2D <- as.data.frame(compositions::clr(as.matrix(sample_2D + 1)))

plot(clr_sample_2D[,1], clr_sample_2D[,2],
     main = "CLR",
     xlab = paste0(substring(colnames(sample_2D)[1], 1,10),"..."),
     ylab = paste0(substring(colnames(sample_2D)[2], 1,10),"..."))