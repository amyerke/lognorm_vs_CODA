#!/usr/bin/env Rscript
# Running the p2_denovo_tree_UPGMA.R with command line args

rm(list = ls()) #clear workspace

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Jones"
cml_scripts <- file.path(home_dir, "r_libraries", "cml_scripts")
r_script <- file.path(cml_scripts, "transformations", "p2_denovo_tree_UPGMA.R")

metad <- file.path(home_dir,project, "SraRunTable.txt")
medtad_delim <- ","

##-Make args for cml script-----------------------------------------##
my_args <- paste(
  "-d", home_dir,
  "-p", project,
  "-m", metad,
  "-l", ",",
  "-r", "Run",
  "-s", "_1000_sd_filtered"
)

# Options:
#   -d CHARACTER, --homedir=CHARACTER
# dataset dir path
# 
# -p CHARACTER, --project=CHARACTER
# project folder
# 
# -m CHARACTER, --metadata=CHARACTER
# metadata file path with filename
# 
# -l CHARACTER, --metadata_delim=CHARACTER
# metadata file deliminator
# 
# -r CHARACTER, --metadata_rowname=CHARACTER
# metadata file row to use for row names
# 
# -s CHARACTER, --outputfilesuffix=CHARACTER
# output_file_suffix
# 
# -f NUMERIC, --filter_level=NUMERIC
# taxonimic level for making otu table 1-6
# 
# -h, --help
# Show this help message and exit

##-Make and run command---------------------------------------------##
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
    print('Opp sa warning is thrown')
    message(cond)
    # Choose a return value in case of warning
    #return(NULL)
  }
)
print("Subscript complete!")

