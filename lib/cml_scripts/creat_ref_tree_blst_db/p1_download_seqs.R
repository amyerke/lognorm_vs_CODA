# Author: Aaron Yerke
print("Script for downloading scripts from silva tree and 
building fasta for blast data base. Expects tree to be in 
home_dir/lib/ref_tree_objs/silva.")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("ape")
}
library("ape")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
library("optparse")
# --------------------------------------------------------------------------
print("Reading cml arguments")
# --------------------------------------------------------------------------
option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character",
                        default=file.path("~","git","lognorm_vs_CODA"),
                        help="dataset dir path"),
  optparse::make_option(c("-t", "--tree_fname"), type="character",
	default="viFy10M5J2nvIBpCLM-QMQ_newick.txt", 
                        help="project folder")
);

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

home_dir <- opt$homedir
output_dir <- file.path(home_dir, "lib", "ref_tree_objs")

setwd(file.path(output_dir))

fastaFile <- "treeFasta.fasta"

# --------------------------------------------------------------------------
print(paste("Build tree fasta at", fastaFile))
# --------------------------------------------------------------------------

tree <- read.tree(file.path(output_dir, "silva",opt$tree_fname))
for (i in 1:seq_along()(tree$tip.label)){
  lab <- tree$tip.label[i]
  id <- strsplit(lab, "_")[[1]][1]
  myGen <- ape:::read.GenBank(id, as.character=T, species.names = FALSE)
  myDna <- lapply(myGen, function(x) paste0(x, collapse = ''))
  cat(paste0(">", id, "\n", paste0(myDna), "\n"), file = fastaFile, append=TRUE)
  Sys.sleep(.45)
}
