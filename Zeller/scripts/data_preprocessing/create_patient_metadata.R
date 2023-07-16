# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for munging metdata files. Requires SraRunTable.csv to create
# patient_metadata.csv

rm(list = ls()) #clear workspace

####-Functions--------------------------------------------------------####
print("Loading custom functions.")


####-Load Depencencies------------------------------------------------####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!requireNamespace("optparse", quietly = TRUE)){
  install.packages("optparse")}
library("optparse")
print("finished loading libraries")

#### Read command line objects ####
option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git',"lognorm_vs_CODA"), 
                        help="dataset dir path", metavar="character")
); 

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

#### Set up constants ####
print("Establishing directory layout and other constants.")
home_dir <- opt$homedir
project <- "Zeller"
output_dir <- file.path(home_dir, project, 'output')

#### Importing and preprocessing tables ####
print("Importing and preprocessing tables")

df <- read.table(file.path(home_dir, project, "SraRunTable.csv"), sep = ",",
                            header = T, check.names = F)
df <- df[df$`Assay Type` == "AMPLICON",]



selected_columns <- c("Run", "Age",	"host_subject_id",
                      "geographic_location_(country_and/or_sea region)",
                      "Collection_date", "AJCC_Stage", "localization",
                      "tissue_type")

df <- df[,selected_columns]

my_missing <- sort(colSums(df == ""), decreasing = TRUE)/nrow(df)*100
print(my_missing)

write.csv(df, file = file.path(home_dir, project, "patient_metadata.csv"), 
          row.names = F)

