# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing lognorm impmentation
# Test 1 uses the "FormwardReads" DADA2 output from Zeller
      
###Defining functions -----------------------------------------------------####
print("Defining functions")

#Loading dependencies ------------------------------------------------------
print("Loading dependencies")
# --------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("testthat")) install.packages("testthat")
library("testthat")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
library("optparse")
####-Reading cml arguments -------------------------------------------------####
print("Reading cml arguments")
# --------------------------------------------------------------------------
option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git','balance_tree_exploration'), 
                        help="dataset dir path"),
  optparse::make_option(c("-p", "--project"), type="character", default="Zeller", 
                        help="project folder")
);

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

#Establishing directory layout and other constants--------------------------
print("Establishing directory layout and other constants.")
# --------------------------------------------------------------------------
home_dir <- opt$homedir
project <- opt$project
working_dir <- file.path(home_dir, "lib", "unit_tests")
base::setwd(working_dir)

#Functions-----------------------------------------------------------------
source(file.path(home_dir, "lib", "table_manipulations.R"))

#Set up constants----------------------------------------------------------
print("Setting up other constants")
# --------------------------------------------------------------------------

if (!file.exists("ForwardReads_DADA2.txt")) unzip("lognorm.zip")

untransformed <- data.frame(data.table::fread(file = file.path("ForwardReads_DADA2.txt"),
                                           header=TRUE, data.table=FALSE), row.names = 1)
# untransformed <- read.table(file = file.path("logNorm","ForwardReads_DADA2.txt"),
#                                               header=TRUE, data.table=FALSE), row.names = 1)
# confirmed_transformed <- read.csv(file = file.path("logNorm","ForwardReads_DADA2LogNorm.txt.zip"),
#                                               header=TRUE, row.names = 1, sep = "\t")
confirmed_transformed <- data.frame(data.table::fread(file = file.path("ForwardReads_DADA2LogNorm.txt"),
                                                      header=TRUE, data.table=FALSE), row.names = 1)
#Testing lognorm-----------------------------------------------------------
test_that("Testing lognorm",
          { expect_equal(lognorm(untransformed), confirmed_transformed)
            expect_equivalent(lognorm(untransformed), confirmed_transformed)
            # expect_identical(lognorm(untransformed), confirmed_transformed)
          })

print("Completed script!")

