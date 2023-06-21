#!/usr/bin/env Rscript
# Author: Aaron Yerke
# Script for creating new subproject directory layout
rm(list = ls()) #clear workspace
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  # args[2] = "out.txt"
}

home_dir <- file.path('~','git',"lognorm_vs_CODA")
#project <- "Noguera-Julian"
project <- args[1]

new_dir <- file.path(home_dir, project)
main_dirs <- c("downloaded_seqs", "output", "scripts")
output_subdirs <- c("graphics", "r_objects", "tables", "taxonomy", "tree_process_blast")
script_subdirs <- c("data_preprocessing", "philr_statistics", "random_forest", "make_ref_tree","transformations", "download_scripts")

{
  if (!dir.exists(home_dir)) {stop(paste("This script expects this project to be located at", home_dir, "."))}
  
  print("home_dir found!")
}

if (!dir.exists(new_dir)){
  dir.create(new_dir)
} else {
  print(paste(new_dir, "already exists!"))
}

for(d in main_dirs){
  my_dir <- file.path(new_dir, d)
  if (!dir.exists(my_dir)){
    dir.create(my_dir)
  } else {
    print(print(paste(my_dir, "already exists!")))
  }
}
print("added main dirs")

for(d in output_subdirs){
  my_dir <- file.path(new_dir, "output", d)
  if (!dir.exists(my_dir)){
    dir.create(my_dir)
  } else {
    print(print(paste(my_dir, "already exists!")))
  }
}
print("added output subdirs")

for(d in script_subdirs){
  my_dir <- file.path(new_dir, "scripts", d)
  if (!dir.exists(my_dir)){
    dir.create(my_dir)
  } else {
    print(print(paste(my_dir, "already exists!")))
  }
}
print("added scripts subdirs")

