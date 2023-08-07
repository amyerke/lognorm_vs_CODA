#!/usr/bin/env Rscript
# Running the p3_parse_blast.R with command line args

rm(list = ls()) #clear workspace

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Jones"
cml_scripts <- file.path(home_dir, "r_libraries", "cml_scripts")
r_script <- file.path(cml_scripts, "make_ref_tree", "p3_parse_blast.R")

##-Make args for cml script-----------------------------------------##
my_args <- paste(
  "-d", home_dir,
  "-p", project,
  "-i", "output_90prcnt_filt.txt",
  "-o", "parsed_output_90prcnt_filt.csv"
)

##-Make and run command---------------------------------------------##
sys_command <- paste(r_script, my_args)
tryCatch({
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

