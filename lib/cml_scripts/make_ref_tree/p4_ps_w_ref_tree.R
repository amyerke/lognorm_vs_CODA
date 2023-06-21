#!/usr/bin/env Rscript
# Author: Aaron Yerke (aaronyerke@gmail.com)
# Convert output from p1_dada2_rd1.R to silva based taxonomy. Loosely following :
#   http://benjjneb.github.io/dada2/training.html

##-cml argument processing------------------------------------------##
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git',"lognorm_vs_CODA"), 
                        help="dataset dir path", metavar="character"),
  optparse::make_option(c("-p", "--project"), type="character", default=NULL, 
                        help="project folder", metavar="character"),
  optparse::make_option(c("-m", "--metadata"), type="character", default=NULL,
                        help="metadata file path with filename", metavar="character"),
  optparse::make_option(c("-l", "--metadata_delim"), type="character", default=NULL,
                        help="metadata file deliminator", metavar="character"),
  optparse::make_option(c("-r", "--metadata_rowname"), type="character", default=NULL,
                        help="metadata file row to use for row names", metavar="character")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

##-load other dependencies------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("dada2",type = "source", checkBuilt = TRUE)
  BiocManager::install("DECIPHER")
}

library("ape")
library("DECIPHER")
library("dada2")
library("phyloseq")
print("external libraries loaded")

print("Establishing directory layout and other constants.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
silva_outf <- file.path(output_dir, "graphics", "silva_ref_tree.pdf")
silva_tree_plot_path <- file.path(output_dir, "graphics", paste0(project, "_ref_tree.pdf"))
silva_ps_robj <- file.path(output_dir, "r_objects","ref_tree_phyloseq_obj.rds")

##-import tables----------------------------------------------------##
myMeta = read.table(opt$metadata,
                    sep=opt$metadata_delim,
                    header=TRUE,
                    row.names = opt$metadata_rowname,
                    check.names = FALSE,
                    stringsAsFactors=FALSE)

print("Imported tables")


##-Import R objects and data preprocessing--------------------------##
con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2.rds"))
seqtab <- readRDS(con)
close(con)

con <- gzfile(file.path( output_dir, "r_objects", "ForwardReads_DADA2_taxonomy.rds"))
taxaTab <- readRDS(con)
close(con)

print("Imported R objects")


##-Build tree------------------------------------------------------##
tree <- phyloseq::read_tree(file.path(home_dir, "taxonomy" , "silva","viFy10M5J2nvIBpCLM-QMQ_newick.txt"))

print(paste("outputting file to", silva_outf))
pdf(file = silva_outf)
phyloseq::plot_tree(tree, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left")
dev.off()

tree_key <- read.table(file.path(output_dir, "tree_process_blast", "parsed_output.csv"),
                      sep = ",",
                      row.names = 1,
                      header = T)

match_df <- data.frame(
                      seqs=character(),
                      num_seqs=integer(),
                      stringsAsFactors=FALSE)

unmatched_ids <- c()

new_labels <- c()

for (i in 1:length(tree$tip.label)){
  old_lab = tree$tip.label[i]
  id = unlist(strsplit(old_lab, "_"))[1]
  if (id %in% tree_key$sseqid){
    
    index = which(tree_key$sseqid == id)[1]
    new_lab = row.names(tree_key)[index]
    tree$tip.label[i] = new_lab
    
    if (id %in% row.names(match_df)){
      print("in if")
      seqs = match_df[id,"seqs"]
      count = match_df[id,"num_seqs"]
      new_seqs = paste(seqs, new_lab)
      match_df[id,"seqs"] = new_seqs
      match_df[id,"num_seqs"] = count + 1
      
    }else{
    new_row = data.frame(
      seqs=new_lab,
      num_seqs=1,
      stringsAsFactors=FALSE)
    row.names(new_row) = c(id)
    match_df = rbind(match_df, new_row)
    }

  }else{
    unmatched_ids = c(unmatched_ids, id)
  }
  # index = which(tree_key$sseqid == id)
  # print("not in")
}

print(paste("unmatched ids:", length(unmatched_ids)))

print(paste("num tree tips pre pruning:", length(tree$tip.label)))
  
tree <- phyloseq::prune_taxa(row.names(tree_key), tree)

print(paste("num tree tip.label post pruning:", length(tree$tip.label)))

pdf(file = silva_tree_plot_path)
phyloseq::plot_tree(tree, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left")
dev.off()

print(paste("tree plotted to", silva_tree_plot_path))

print(paste("num duplicated tips:", sum(duplicated(tree$tip.label))))

print(paste("duplicated tips:", tree$tip.label[duplicated(tree$tip.label)]))

if( sum(duplicated(tree$tip.label)) > 0){
  print("removing duplicated tip")
  tree <- ape::drop.tip(tree, tree$tip.label[duplicated(tree$tip.label)], trim.internal = TRUE, subtree = FALSE,
                root.edge = 0, rooted = is.rooted(tree), collapse.singles = TRUE,
                interactive = FALSE)
}

print(paste("num duplicated rows taxatab: ", sum(duplicated(row.names(taxaTab)))))
print(paste("num duplicated cols taxatab: ", sum(duplicated(colnames(taxaTab)))))

ps <- phyloseq::phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(myMeta),
               tax_table(taxaTab),
               phy_tree(tree)
)

print(paste("outputing Silva ref phyloseq object to:", silva_ps_robj))
saveRDS(ps, file = silva_ps_robj)

print("Rscript complete.")
