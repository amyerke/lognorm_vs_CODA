#!/usr/bin/env Rscript
# script for parsing the blast results to build a tree

# Columns follow this pattern:
# qseqid sseqid pident length evalue bitscore score ppos
# 
# 1.	qseqid	 query (e.g., unknown gene) sequence id
# 2.	sseqid	 subject (e.g., reference genome) sequence id
# 3.	pident	 percentage of identical matches
# 4.	length	 alignment length (sequence overlap)
# 5.	evalue	 expect value
# 6.	bitscore	 bit score
# 7.  score     Raw score
# 8.  ppos      Percentage of positive-scoring matches
# Comparison stategy: compare ppos (if tie, go with biggest alignment length)

print("Begining cml argument processing.")
####-read in external libraries---------------------------------------####
if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")
}
library("optparse")

####-cml argument processing------------------------------------------####
option_list <- list(
  make_option(c("-d", "--homedir"), type="character", 
              default=file.path('~','git',"lognorm_vs_CODA"), 
              help="dataset dir path"),
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="project folder name in homedir"),
  make_option(c("-i", "--input_file"), type="character", default="output.txt", 
              help="blast output from p2"),
  make_option(c("-o", "--output"), type="character", 
              default="parsed_output.csv", 
              help="name of output file, will be found in tree_process_blast")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

print(opt)

####-Establish directory constants------------------------------------------####
print("Establishing directory layout.")
home_dir <- opt$homedir
project <- opt$project
output_dir <- file.path(home_dir, project, 'output')
text_output_fn <- paste0("blast_stats_", project, "_", basename(opt$input_file), ".log")

setwd(file.path(output_dir, "tree_process_blast"))
print(paste0("Creating ", text_output_fn, "in", getwd()))
file.create(text_output_fn)

input_file <- opt$input_file

df <- data.frame(#qseqid=character(),#empty df to fill with parsed results
                 sseqid=character(),
                 bitscore=integer(),
                 count_matched=integer(),
                 stringsAsFactors=FALSE)

print("Parsing blast output.")
all_qseq <- c()
con <- file(input_file, "r")
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  result = unlist(strsplit(line, '\t'))
  qseq = result[1]
  sseq = unlist(strsplit(result[2],"\\|"))[2] #dbj|AB064923|
  evalue = as.numeric(result[5])
  bitsc = as.integer(result[6])
  # print(paste(sseq, bitsc))
  all_qseq <- unique(c(all_qseq, qseq))
  if (evalue < 10^-10){
    if ( qseq %in% row.names(df)){ #check if query already in df
      if (bitsc > df[qseq, "bitscore"]){#higher bitscore is better
        count = df[qseq,"count_matched"]
        df[qseq,] = list(sseq, bitsc, count + 1)#increase bitsc & count of existing row
      }
    }else{#if qseq not already in df
      newRow <- data.frame(sseqid=sseq,
                           bitscore=bitsc,
                           count_matched=1) 
      row.names(newRow) = qseq
      df = rbind(df, newRow)
    }
  }#close if eval
}# close while
close(con)

####-Record results---------------------------------------------------------####
print(paste("original nrow:", nrow(df)))
write(paste("original nrow:", nrow(df)),
      file=text_output_fn,
      append=TRUE)

print(paste("number of unique hits:", length(unique(rownames(df)))))
write(paste("number of unique hits:", length(unique(rownames(df)))),
      file=text_output_fn,
      append=TRUE)

print(paste("number of unique queries:", length(unique(all_qseq))))
write(paste("number of unique queries:", length(unique(all_qseq))),
      file=text_output_fn,
      append=TRUE)

mapped_fraction <- length(unique(rownames(df)))/length(unique(all_qseq))
print(paste("mapped_fraction:", mapped_fraction))
write(paste("mapped_fraction:", mapped_fraction),
      file=text_output_fn,
      append=TRUE)

myT <- table(df[,"sseqid"])

print(paste("ave seq/node:", mean(myT), "\nmax seq/node:", max(myT)))
write(paste("ave seq/node:", mean(myT), "\nmax seq/node:", max(myT)),
      file=text_output_fn,
      append=TRUE)

pdf(file = file.path(output_dir, "graphics", "p3_parse_blast.pdf"))
hist(myT, breaks = 150, xlab = "Sequences per node tip", main = "Histogram of seqs per node tip")
barplot(myT, las = 2, xlab = "Sequences per node tip", main = "Histogram of seqs per node tip")
dev.off()

write.csv(df, file = opt$output)

print("Script complete!")
