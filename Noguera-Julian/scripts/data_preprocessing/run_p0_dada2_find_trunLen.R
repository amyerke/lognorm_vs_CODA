#!/usr/bin/env Rscript
# Running the p0_dada2_find_trunLen.R with command line args
# Note: download result with scp amy@hpc.uncc.edu:~/git/lognorm_vs_CODA/McDonald/output/graphics/plotQualF.png .


rm(list = ls()) #clear workspace

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Noguera-Julian"

r_script <- file.path(home_dir, "lib", "cml_scripts", "data_preprocessing","p0_dada2_find_trunLen.R")

my_args <- paste(
  "-d", home_dir,
  "-p", project,
  "-f", "_1.fastq.gz",
  "-r", "_2.fastq.gz"
)

sys_command <- paste(r_script, my_args)

tryCatch(
  { 
    system(sys_command,
           intern = FALSE,
           ignore.stdout = FALSE, ignore.stderr = FALSE,
           wait = TRUE, input = NULL,
           minimized = FALSE, invisible = TRUE, timeout = 0)
  },
  error=function(cond) {
    print('Opps, an error is thrown')
    message(cond)
  },
  warning=function(cond) {
    print('Opps! warning is thrown')
    message(cond)
    # Choose a return value in case of warning
    #return(NULL)
  }
)
print("Subscript complete!")

