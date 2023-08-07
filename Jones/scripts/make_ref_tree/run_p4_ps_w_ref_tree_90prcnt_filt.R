#!/usr/bin/env Rscript
# Running the p4_ps_w_ref_tree.R with command line args

rm(list = ls()) #clear workspace

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Jones"
cml_scripts <- file.path(home_dir, "lib", "cml_scripts")
r_script <- file.path(cml_scripts, "make_ref_tree", "p4_ps_w_ref_tree.R")

metad <- file.path(home_dir,project, "SraRunTable.txt")

##-Make args for cml script-----------------------------------------##
my_args <- paste(
  "--homedir", home_dir,
  "--project", project,
  "--metadata", metad,
  "--metadata_delim", ",",
  "--metadata_rowname", "Run",
  "--modified_LTP_fn", "ltp_90prcnt_filt_tree.pdf",
  "--ps_LTP_rds", "ltp_90prcnt_filt.rds",
  "--seq_count_table", "filtered_90prcnt_dada2.rds",
  "--tree_key", "parsed_output_90prcnt_filt.csv"
)

# Options:
# 	-d HOMEDIR, --homedir=HOMEDIR
# 		dataset dir path

# 	-p PROJECT, --project=PROJECT
# 		project folder

# 	-m METADATA, --metadata=METADATA
# 		metadata file path with filename

# 	-l METADATA_DELIM, --metadata_delim=METADATA_DELIM
# 		metadata file deliminator

# 	-r METADATA_ROWNAME, --metadata_rowname=METADATA_ROWNAME
# 		metadata file row to use for row names

# 	-g REF_TREE_PDF, --ref_tree_pdf=REF_TREE_PDF
# 		output file name for the unmodified LTP tree pdf
#                         - will be in graphics folder

# 	-u MODIFIED_LTP_FN, --modified_LTP_fn=MODIFIED_LTP_FN
# 		output file name for the modified LTP tree pdf
#                         - will be in graphics folder

# 	-o PS_LTP_RDS, --ps_LTP_rds=PS_LTP_RDS
# 		output file name for for final phyloseq object
#                         -will go in r_objects folder

# 	-s SEQ_COUNT_TABLE, --seq_count_table=SEQ_COUNT_TABLE
# 		Counts table r object in r_objects folder

# 	-t TAXONOMY, --taxonomy=TAXONOMY
# 		Taxonomy table

# 	-k TREE_KEY, --tree_key=TREE_KEY
# 		metadata file deliminator

# 	-h, --help
# 		Show this help message and exit

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
    print('Opps! warning is thrown')
    message(cond)
    # Choose a return value in case of warning
    #return(NULL)
  }
)
print("Subscript complete!")

