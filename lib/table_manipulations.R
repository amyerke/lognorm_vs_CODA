# Author: Aaron Yerke (aaronyerke@gmail.com)
# Library of table manipulating functions

##-Functions--------------------------------------------------------##
philr_tutorial_normalization <- function(df) {
  # The philr tutorial https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R uses this block for normalizations: 
  # ps <-  filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
  # ps <-  filter_taxa(ps, function(x) sd(x)/mean(x) > 3.0, TRUE)
  # ps <- transform_sample_counts(ps, function(x) x+1)
  df <- df[apply(df, 2, function(x) sum(x > 3) > (0.2*length(x)))]
  df <-  df[apply(df, 2, function(x) sd(x)/mean(x) > 3.0)]
  df <- df + 1
  return(df)
}

lognorm <- function(table){
  # This is Fodor's custom lognormal method.
  # table<-table[rowSums(table)>1000,]
  average<-sum(rowSums(table))/nrow(table)
  table<-sweep(table,1,rowSums(table),"/")
  table<-log10(table*average + 1)
  return(table)
}

proportions_transform <- function(table){
	# RC/n where RC = raw count and n = number of
	# sequences in a sample
  table<-sweep(table,1,rowSums(table),"/")
  return(table)
}

heilinger_transform <- function(table){
	# √(RC/n) where RC = raw count and n = number of
	# sequences in a sample
  table<-sweep(table,1,rowSums(table),"/")
  table<-sqrt(table)
  return(table)
}

filt_seq_dpth <- function(min_read_depth, df) {
	#function to remove samples below min read depth
	#because vegan's rrarefy doesn't seem to be doing it
  df <- df[rowSums(df) > min_read_depth, ]
	return(df)
}
