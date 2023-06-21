#!/usr/bin/env Rscript
# Running the p1_dada2_rd1.R with command line args

rm(list = ls()) #clear workspace

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Jones"

r_script <- file.path(home_dir, "r_libraries", "cml_scripts", "data_preprocessing","p1_dada2_rd1.R")

my_args <- paste(
  "-d", home_dir,
  "-p", project,
  "-f", "_1.fastq.gz",
  "-r", "_2.fastq.gz",
  "-t", "240"
)

sys_command <- paste(r_script, my_args)

system(sys_command,
       intern = FALSE,
       ignore.stdout = FALSE, ignore.stderr = FALSE,
       wait = TRUE, input = NULL,
       minimized = FALSE, invisible = TRUE, timeout = 0)
