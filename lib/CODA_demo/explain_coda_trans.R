# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script to demonstrate compsitional data transformations
# Fasta file was run through https://www.ebi.ac.uk/Tools/msa/clustalo/ using the phylip option
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-functions--------------------------------------------------------##
make_1d_graph <- function(df_1d, f_path, m_title){
  dat <- data.frame("data" = as.numeric(df_1d[1,]), 
                    "colrs" = colnames(df_1d),
                    "x_ax" = rep(1,ncol(df_1d)))
  g <- ggplot2::ggplot(data = dat,
                       aes(x="x_ax", y=data, col = as.factor(colrs),
                           label=as.factor(colrs))) +
    ggplot2::geom_point(show.legend = FALSE, size = 5) +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::ggtitle(m_title) +
    # ggplot2::geom_text(aes(color = "black")) +
    ggplot2::theme(axis.ticks.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank()) +
    ggplot2::theme_bw()
  ggplot2::ggsave(filename = f_path,
                  height = 4, width = 1, plot = g)
  print(g)
}
g.rowMeans <- function(y, p=rep(1, nrow(y))){
  sp <- sum(p) # as given in text on page 4
  exp(1/sp*rowSums(log(y)%*%diag(p)))
}

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
# library("ALDEx2")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library("RColorBrewer")
library("compositions")
library("phyloseq")
# library("vegan")
# library("DESeq2")
library("philr")
library("ape")
library("ggplot2")
library("ape")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git','balance_tree_exploration')
project <- "demo"
#home_dir <- file.path('cloud','project')
output_dir <- file.path(home_dir, project, 'output')
setwd(file.path(home_dir))

##-Functions--------------------------------------------------------##
source(file.path(home_dir, "lib", "statistical_functions.R"))
source(file.path(home_dir, "lib", "table_manipulations.R"))

##-Import tables and data preprocessing-----------------------------##

# The palette with black:
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


##-Create demo data-------------------------------------------------##
set.seed(77)
samp_size <- 8
alr_denom <- samp_size
raw_reads <- base::sample(0:100, size=samp_size, replace = T)
names(raw_reads) <- cbbPalette
alr_reads <- compositions::alr(raw_reads, ivar=alr_denom)
##------------------------------------------------------------------##
print(paste0("The alr denominator is index ", alr_denom, " and the raw read valule is ", raw_reads[alr_denom], "."))
##------------------------------------------------------------------##
print(paste0("Manual check of first value of alr transform: ", log((raw_reads[1]/raw_reads[alr_denom])) == alr_reads[1],"."))
##------------------------------------------------------------------##

clr_reads <- compositions::clr(raw_reads)
##-Test the testing data----------------------------------------------##

n_alr_reads <- raw_reads[! raw_reads %in% raw_reads[alr_denom]]

for(i in 1:length(raw_reads)-1){
  n_alr_reads[i] <- log(n_alr_reads[i]/raw_reads[alr_denom])
}

##------------------------------------------------------------------##
print(paste0("Another check of alr transform: ",base::identical(alr_reads, n_alr_reads),"."))
##------------------------------------------------------------------##

##------------------------------------------------------------------##
print(paste0("Check of clr transform (these should all be the same):"))
##------------------------------------------------------------------##
for(i in 1:length(raw_reads)){
  print(paste( log((raw_reads[i]/compositions::geometricmean(raw_reads))), clr_reads[i] ))
}

##-Make some graphs-------------------------------------------------##
palette( RColorBrewer::brewer.pal(8,"Dark2") )

##------------------------------------------------------------------##
print(paste0("Finding min and max for graphs to keep them on the same scale."))
##------------------------------------------------------------------##
min_y <- min(c(raw_reads,alr_reads,clr_reads))
max_y <- max(c(raw_reads,alr_reads,clr_reads))
#I ultimately decided not to use min/max y as I didn't like the aesthetic. Might go back later.
make_dot_plots <- function(named_num_vec, 
                           title_text, 
                           pdf_name,
                           min_y = min_y,
                           max_y = max_y) {
  dat <- data.frame("data" = as.numeric(named_num_vec), 
                    "colrs" = names(named_num_vec),
                    "x_ax" = rep(1,length(named_num_vec)))
  g <- ggplot2::ggplot(data = dat,
                       aes(x=x_ax, y=data, col=colrs)) +
    ggplot2::geom_point(show.legend = FALSE, size = 5) +
    # ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::ggtitle(title_text) +
    ggplot2:: scale_colour_manual("named_num_vec", values=as.character(names(named_num_vec))) +
    # ggplot2::geom_text() +
    ggplot2::scale_x_continuous(breaks = NULL) +
    scale_fill_manual(values=as.character(names(named_num_vec))) + 
    ggplot2::theme(axis.ticks.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank()) +
    ggplot2::theme_bw()
  ggplot2::ggsave(filename = file.path(output_dir, pdf_name),
                  height = 4, width = 1, plot = g)
  # print(g)
}

##------------------------------------------------------------------##
print(paste0("Finding min and max for graphs to keep them on the same scale."))
##------------------------------------------------------------------##
make_dot_plots(named_num_vec = raw_reads,
               title_text = "raw",
               pdf_name = paste0(project, "_new1_raw_demo.pdf"))
make_dot_plots(named_num_vec = alr_reads,
               title_text = "alr",
               pdf_name = paste0(project, "_new1_alr_demo.pdf"))
make_dot_plots(named_num_vec = clr_reads,
               title_text = "clr",
               pdf_name = paste0(project, "_new1_clr_demo.pdf"))


##-Import tables and data preprocessing-----------------------------##

demo_tree_fasta <- "
>a
TGCTAAGGCCGATGGCGACCGGCGCA
>b
CGCTAAGTTTGATGGCGACCGGCGCA
>c
TGCTAAGTTTGATGGCGACCGGCGCA
>d
TGCTAAGGCTGATGGCGACCGGCGCA
>e
TGCTTTCTCTTGCTGGCGACCGGCGCA
>f
TTGCAAACCAAAGCTGGCGACCAGCG
>g
TCCTTCGGGACTGATTATTTTGTGAC
>h
TGCTAAGGCTGATGGCGACCGGCGCT
"
#deom_tree was made by passing the demo_fasta through https://www.ebi.ac.uk/Tools/msa/clustalo/ using the phylip option

demo_tree <- ape::read.tree(file.path(home_dir, "demo", "clustalo-p2m-phylip.dnd"))
demo_tree <- ape::makeNodeLabel(demo_tree, method="number", prefix='n')

demo_taxa <- read.table(file.path(home_dir, "demo", "demo_tax.csv"), sep=",",
                        row.names = 1, header = TRUE)

demo_phylo <- phyloseq::phyloseq(otu_table(as.matrix(raw_reads), taxa_are_rows = T),
                                 tax_table(as.matrix(demo_taxa)),
                                 phy_tree(demo_tree))

demo_phylo@phy_tree <- makeNodeLabel(demo_phylo@phy_tree, method="number", prefix='n', )

pdf(file = file.path(output_dir, "demo_tree.pdf"), width = 6, height = 6)
ape::plot.phylo(demo_phylo@phy_tree, 
                show.node.label = TRUE, 
                show.tip.label = TRUE,
                use.edge.length = TRUE,
                root.edge = TRUE,
                # no.margin = TRUE,
                font = 100,
                direction = "upwards"
                )
# edgelabels(round(demo_phylo@phy_tree$edge.length, 2), 
#            font = -1, 
#            adj = c(-0.3,0.5),
#            col = "blue",
#            frame = "none")
dev.off()

name.balance(demo_phylo@phy_tree, demo_phylo@tax_table, 'n1')

##-philr transform--------------------------------------------------##
demo_philr <- philr::philr(t(demo_phylo@otu_table), demo_phylo@phy_tree,
                  part.weights='uniform',
                  ilr.weights='uniform',
                  return.all = FALSE)

##-plot the philr transform-----------------------------------------##
dat <- data.frame("data" = as.numeric(demo_philr[1,]),
                  "colrs" = colnames(demo_philr),
                  "x_ax" = rep(1,7))
g <- ggplot2::ggplot(data = dat,
                     aes(x="x_ax", y=data,
                         label=as.factor(colrs))) +
  ggplot2::geom_point(show.legend = FALSE, size = 5, shape = 0) +
  ggplot2::scale_color_brewer(palette="Dark2") +
  ggplot2::ggtitle() +
  ggplot2::geom_text(aes(color = "red"),
                     show.legend = FALSE) +
  ggplot2::theme(axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank()) +
  ggplot2::theme_bw()
ggplot2::ggsave(filename = file.path(output_dir, 
                                     paste0(project, "")),
                height = 4, width = 1, plot = g)
print(g)

##-recreating steps of philr----------------------------------------##
#change everything to a proportion that adds to 1
demo_minicl <- philr::miniclo(t(demo_phylo@otu_table))
make_1d_graph( df_1d = demo_minicl@.Data, 
               m_title = "minclos", 
               f_path = file.path(output_dir, paste0(project, "_minclos_demo.pdf")))

#shift everything by the part weight (in our case, just 1)
demo_mincl_shift <- philr::shiftp(demo_minicl, rep(1, ncol(t(demo_phylo@otu_table))))
make_1d_graph( df_1d = demo_mincl_shift, 
               m_title = "minclos+1", 
               f_path = file.path(output_dir, paste0(project, "_minclosureUnif_demo.pdf")))



s <- 0.1534247 
r <- 0.04657534
k <- 0.1726027

sqrt((r * s)/(r+s)) * log((geometricmean(s)) / (geometricmean(r)))

#n7
sqrt((1 * 1)/(1+1)) * log((geometricmean(17)) / (geometricmean(56)))
#n6 h,a,d
sqrt((1 * 2)/(1+2)) * log((geometricmean(63)) / (geometricmean(c(17,56))))
#n4 l = h,a,d r = c,b | l = 63, 17, 56 r = 68, 44
sqrt((3 * 2)/(3+2)) * log((geometricmean(c(68,44))) / (geometricmean(c(63,17,56))))

sbp <- phylo2sbp(demo_phylo@phy_tree)

ilrBase <- buildilrBasep(sbp, rep(1,8))
write.table(round(ilrBase, 2), file = file.path(output_dir, "ilrBase.csv"), sep=",")

ilrp(demo_mincl_shift, rep(1,8), ilrBase)

#philr transform (trying to break it down into simple steps)
pt1 <- log(demo_mincl_shift/g.rowMeans(demo_mincl_shift,rep(1,8)))%*%diag(rep(1,8))%*%ilrBase

pt2 <- log(demo_mincl_shift/g.rowMeans(demo_mincl_shift,rep(1,8))) %*% ilrBase

pt3 <- log(demo_minicl/exp(1/sum(rep(1,8)) * rowSums(log(demo_minicl)) )) %*% ilrBase


just_clr <- log(demo_minicl/exp(1/sum(rep(1,8)) * rowSums(log(demo_minicl)) ))


x <- seq(0,1,0.05)
y <- seq(1,0,-0.05)

plot(x,y)

# https://www.rdocumentation.org/packages/phytools/versions/0.7-70/topics/rotateNodes
# http://blog.phytools.org/2013/07/rotating-all-nodes-in-tree-or-set-of.html

