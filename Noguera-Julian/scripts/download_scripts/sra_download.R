# Author: Aaron Yerke
# Script to download from the SRA
# This was helepful: https://github.com/jsilve24/philr/blob/master/vignettes/philr-intro.Rmd#L142

rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##

print("Establishing constants")
home_dir <- file.path('~','git',"lognorm_vs_CODA")
project <- "Noguera-Julian"
download_dir <- file.path(home_dir, project, "downloaded_seqs")
if (!dir.exists(download_dir)) dir.create(download_dir) #create dir if needed
print(paste("Download destination:", download_dir))

sra_run_table <- read.table(file.path(home_dir, project, "SraRunTable.txt"),
                            sep = ",",
                            header = TRUE)

print("Creating SRR list.")
my_rows <- sra_run_table$Assay.Type == "AMPLICON"
my_accessions <- sra_run_table$Run[my_rows]

print(paste("number of rows:", length(my_accessions)))

downloaded_files <- list.files(path = download_dir)
print(paste("number of files already downloaded:", length(downloaded_files)))

print(paste("need to download this many more:", length(my_accessions) - length(downloaded_files)))

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
  }else{
    print(paste(run, "was already there!"))
  }
  Sys.sleep(1)  
}  

downloaded_files <- list.files(path = download_dir)
print(paste("number of files after running:", length(downloaded_files)))

print("Script completed.")
# bash download commands:
# module load sra-tools
# nohup fasterq-dump $SRR -O ../../downloaded_seqs
# 
# Output format:
#   SRR5799852_1.fastq.gz



