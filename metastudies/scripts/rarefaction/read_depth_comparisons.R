# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making looking at read depth levels for rarefaction

rm(list = ls()) #clear workspace

####-Load dependencies------------------------------------------------####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
if (!requireNamespace("vegan", quietly = TRUE)) BiocManager::install("vegan")
library("vegan")
print("finished loading libraries")

####-Establish directory layout---------------------------------------####
home_dir <- file.path('~','git',"lognorm_vs_CODA")
projects <- c("Jones", "Zeller", "Vangay", "Noguera-Julian")
output_dir <- file.path(home_dir, "metastudies", "output")
#asv tables are in each project in "project"/output/tables/ForwardReads_DADA2.txt
asv_tables <- file.path("output", "tables", "ForwardReads_DADA2.txt")
print("Finished establishing directory layout")


####-Information determined on first run of this script---------------####
lowest_ranks <- c(11,2,12,10)#lowest ranks to use downstream, in the same order as "projects"

####-Read in sequence tables------------------------------------------####
for (proj in 1:length(projects)){
  project <- projects[proj]
  asv_path <- file.path(home_dir,project,asv_tables)
  asv_table <- data.frame(data.table::fread(file = asv_path, check.names=FALSE,
                                            header=TRUE, data.table=FALSE),
                          row.names = 1, check.names=FALSE)
  my_sums <- rowSums(asv_table)
  min_reads <- min(my_sums)
  max_reads <- max(my_sums)
  rank_reads <- rank(my_sums)
  # print(paste0("Read depth ", project,  "\nmin: ", min_reads, ", max:", max_reads))
  plot(rank_reads, my_sums, 
       main = paste0("Read depth ", project,  "\nmin: ", min_reads, ", max: ",max_reads),
       log = "y", xlab = "Sample rank by read depth", ylab = "Read depth"
       )
  top_50_rank <- sort(rank_reads[rank_reads < 51])
  top_50_reads <- my_sums[order(factor(names(my_sums), levels=names(top_50_rank)))]
  top_50_reads <- top_50_reads[1:50]
  dif_between_read_depth <- diff(top_50_reads)
  plot(top_50_rank, top_50_reads, 
       main = paste0(project, " lowest 50 samples'", "\nRD at elbow:",
                     top_50_reads[lowest_ranks[proj]], " Samples lost: ",
                     top_50_rank[lowest_ranks[proj]]),
       log = "y", xlab = "Sample rank by read depth", ylab = "Read depth", xaxt = "n"
  )
  axis(1, at = seq(1,50,2), las=3, #customize tick marks, las turns the x ticks sideways
       tck = 1, lty = 2, col = "gray")#adds grid lines for visualization
  print(paste(project, "desired read depth:", top_50_reads[lowest_ranks[proj]]))
}
####-Identify results by hand-----------------------------------------####
# Based on the results in the graphs generated above, by hand the following ranks will be 
# used to set read depth threshold.
projects <- c("Jones", "Zeller", "Vangay", "Noguera-Julian")
lowest_ranks <- c(11,2,12,10)#lowest ranks to use downstream


# vegan::rarecurve(asv_table, step = 100)
