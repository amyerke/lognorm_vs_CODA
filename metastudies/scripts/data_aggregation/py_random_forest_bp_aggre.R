# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for making table of metadata characterics

rm(list = ls()) #clear workspace

##-Load Depencencies------------------------------------------------##
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("qpdf", quietly = TRUE)) install.packages("qpdf")
print("finished loading libraries")

##-Establish directory layout---------------------------------------##
home_dir <- file.path('~','git',"lognorm_vs_CODA")
projects <- c("Jones", "Zeller", "Vangay", "Noguera-Julian")
rf_bp_filename <- "bp_sklearn_random_forest_manual_0.75train.pdf"
output_path <- file.path(home_dir, "dissertation_supplement", "rf_bp_py","sklearn_rf_bp.pdf")
##-Functions--------------------------------------------------------##

##-Create tree attribute vectors------------------------------------##
pdf_paths <- c()
for (proj in 1:length(projects)) {
  project <- projects[proj]
  my_pdf <- file.path(home_dir, project, "output", "graphics", rf_bp_filename)
	pdf_paths <- c(pdf_paths, my_pdf)
}

qpdf::pdf_combine(input = pdf_paths, output_path)



