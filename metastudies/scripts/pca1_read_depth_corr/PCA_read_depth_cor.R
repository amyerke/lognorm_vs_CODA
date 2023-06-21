# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for demonstrating that read depth drives PCA1

rm(list = ls()) #clear workspace

# --------------------------------------------------------------------------
print("Defining functions")
# --------------------------------------------------------------------------

# Loading dependencies------------------------------------------------------
print("Loading dependencies")
# --------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("ggplot2")
if (!requireNamespace("rgr", quietly = TRUE)) install.packages("rgr")
library("rgr")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")

# Establishing directory layout and other constants-------------------------
print("Establishing directory layout and other constants.")
# --------------------------------------------------------------------------
home_dir <- file.path('~','git',"lognorm_vs_CODA")
projects <- c("Jones", "Vangay", "Noguera-Julian", "Zeller")
output_dir <- file.path(home_dir, "metastudies", 'output')
mds_depth <- 5

# Creating empty vectors to fill in main loop-------------------------------
print("Creating empty vectors to fill in main loop.")
# --------------------------------------------------------------------------

kendal <- vector(mode = "numeric", length = length(projects) * mds_depth)
prj_nam <- vector(mode = "character", length = length(projects) * mds_depth)
mds_lev <- vector(mode = "integer", length = length(projects) * mds_depth)

wide_header <- unlist(lapply(1:mds_depth, function(x) paste0("PCA",x)))
wide_dframe <- data.frame(matrix(nrow = length(projects), ncol = length(wide_header)),
                          row.names = projects)
colnames(wide_dframe) <- wide_header

# Running main loop---------------------------------------------------------
print("Running main loop.")
# --------------------------------------------------------------------------
counter <- 1
for (pr in 1:length(projects)){
  prj <- projects[pr]
  print(paste("Starting", prj))
  proj_out <- file.path(home_dir, prj, "output")
  my_table <- readRDS(file.path(proj_out, "r_objects", "ForwardReads_DADA2.rds"))
  read_depth <- base::rowSums(my_table)
  ##-Create a PCA-----------------------------------------------------##
  my_prcmp <- prcomp(my_table, 
                     center = TRUE,
                     rank = mds_depth)
  ##-Extract PCA matrix and convert to dataframe----------------------##
  myPCA <- data.frame(my_prcmp$x)
  my_var_exp <- my_prcmp$sdev^2/sum(my_prcmp$sdev^2)
  wide_row <- c()
  for (md in 1:mds_depth){
    my_corr <- stats::cor(read_depth, myPCA[,md], method = "kendal")
    wide_row <- c(wide_row, my_corr)
    kendal[counter] <- my_corr
    prj_nam[counter] <- prj
    mds_lev[counter] <- md
    counter <- counter + 1
  }#for loop mds
  wide_dframe[pr,] <- wide_row
}

long_dframe <- data.frame(kendal, prj_nam, mds_lev)

# Saving tables-------------------------------------------------------------
print("Saving tables.")
# --------------------------------------------------------------------------

write.table(long_dframe,
            sep = ",",
            row.names = FALSE,
            file = file.path(output_dir, "long_pca_readdepth_corr.csv"))

write.table(wide_dframe,
            sep = ",",
            file = file.path(output_dir, "wide_pca_readdepth_corr.csv"))

write.table(format(wide_dframe,  digits=3),
            sep = ",",
            file = file.path(output_dir, "abrev_wide_pca_readdepth_corr.csv"))



