# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making table of metadata characterics

rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
print("finished loading libraries")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
projects <- c("Jones", "Zeller","Vangay", "Noguera-Julian")
proj_metaname <- c("patient_metadata.tsv", "patient_metadata.csv", "patient_metadata.tsv", "patient_metadata.tsv")
proj_meta_delim <- c("\t",",","\t","\t")
proj_meta_colname <- c("Run", "Run", "run_accession", "Run")
final_output_dir <- file.path(home_dir, "metastudies", "output")
##-Functions--------------------------------------------------------##

##-Create tree attribute vectors------------------------------------##
binary_meta <- c()
multi_meta <- c()
continous_meta <- c()

for (proj in 1:length(projects)) {
	project <- projects[proj]
	print(project)
	proj_out <- file.path(home_dir, prj, "output")
	metadata <- data.frame(data.table::fread(file = file.path(proj_out, proj_metaname[proj]),
																					header=TRUE, data.table=FALSE), row.names = proj_meta_colname[proj])
	my_table <- readRDS(file.path(proj_out, "r_objects", "ForwardReads_DADA2.rds"))
	read_depth <- base::rowSums(my_table)
	metadata <- metadata[row.names(my_table),]
	
	#Correlations for continuous data and pvalues for factors
	corr_meta <- c()
	corr_value <- c()
	pval_meta <- c()
	pval_val <- c()
	pval_adj_pval <- c()
	
	meta_class <- sapply(metadata,class)
	for (colmn in 1:ncol(metadata)){
		if (meta_class[colmn] == "integer" || meta_class[colmn] == "numeric"){
		  corr_meta <- meta_class
			continous_meta <- c(continous_meta, colnames(metadata)[colmn])
		}else{
		  # print(unique(metadata)[,colmn])
			

			}
		}
	}
}


