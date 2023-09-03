# Author: Aaron Yerke
# Script to download from the SRA

rm(list = ls()) #clear workspace

####-Load dependencies------------------------------------------------####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
print("finished loading libraries")

####-Establish directory layout---------------------------------------####
home_dir <- file.path("~", "git", "lognorm_vs_CODA")
project <- "Vangay"
download_dir <- file.path(home_dir, project, "downloaded_seqs")
print(paste("Download destination:", download_dir))
if (!dir.exists(download_dir)) dir.create(download_dir) #create dir if needed
sra_path <- file.path(home_dir, project, "fullMetadata.tsv")
sra_run_table <- read.table(sra_path,
                            sep = "\t", check.names = FALSE,
                            header = TRUE)

print("Established constants.")

####-Processing data and downloading----------------------------------####

print(paste("There are", nrow(sra_run_table),
"rows and", ncol(sra_run_table), "columns in", sra_path))

print("Creating SRR list.")
my_accessions <- sra_run_table$run_accession

print(paste("number of rows:", length(my_accessions)))

downloaded_files <- list.files(path = download_dir)
print(paste("number of files already downloaded:", length(downloaded_files)))

print(paste("need to download this many more:", 
length(my_accessions) - length(downloaded_files)))

##-Download SRA files ----------------------------------------------##
setwd(download_dir)
for (run in my_accessions) {
  my_file <- paste0(run, "_1.fastq.gz")
  if (!my_file %in% downloaded_files){
    print(paste("attempting download of:", run))
    my_command <- paste("module load sra-tools ;",
                        "fasterq-dump -S", run)
    print(paste("my command:", my_command))
    system(command = my_command, wait = TRUE)
  }else {
    print(paste(run, "was already there!"))
  }
  Sys.sleep(1)
}

downloaded_files <- list.files(path = download_dir)
print(paste("number of files after running:", length(downloaded_files)))

print("Script completed.")
