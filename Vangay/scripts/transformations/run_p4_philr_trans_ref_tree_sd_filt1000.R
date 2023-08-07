#!/usr/bin/env Rscript
# Author: Aaron Yerke (aaronyerke@gmail.com)
# This is a pipeline that was created to run the philr transform p4_philr_transform.R in the cml scripts.
# Many resources were used to create this pipeline:
#   Processing the sequences through DADA2, making trees with phangorn, and putting them into phyloseq objects:
#     This is by far the biggest source:
#     https://github.com/spholmes/F1000_workflow/blob/master/MicrobiomeWorkflow/MicrobiomeWorkflowII.Rmd
#     or (in article format)
#     https://f1000research.com/articles/5-1492/v2
#     
#   DADA2:
#     https://benjjneb.github.io/dada2/tutorial.html
#   
#   Using philr:
#     https://bioconductor.org/packages/release/bioc/html/philr.html
#     Code example:
#       https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.R
rm(list = ls()) #clear workspace

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Vangay"
cml_scripts <- file.path(home_dir, "lib", "cml_scripts")
r_script <- file.path(cml_scripts, "transformations", "p4_philr_transform.R")

metad <- file.path(home_dir,project, "SraRunTable.txt")
medtad_delim <- ","

##-Make args for cml script-----------------------------------------##
my_args <- paste(
  "-d", home_dir,
  "-p", project,
  "-l", "1000",
  "-f", "ref_tree_phyloseq_obj.rds",
  "-s", "_ref_tree_1000_sd_flt"
)

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
    print('Opps, a warning is thrown')
    message(cond)
    # Choose a return value in case of warning
    #return(NULL)
  }
)
print("Subscript complete!")

