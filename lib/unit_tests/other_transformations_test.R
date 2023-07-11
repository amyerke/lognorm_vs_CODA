# Author: Aaron Yerke (aaronyerke@gmail.com)
# Script for testing lognorm impmentation
# Test 1 uses the "FormwardReads" DADA2 output from Zeller

####Loading dependencies ------------------------------------------------------####
print("Loading dependencies")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("testthat")) install.packages("testthat")
library("testthat")
if (!requireNamespace("data.table", quietly = TRUE)) BiocManager::install("data.table")
library("data.table")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
library("optparse")
####-Reading cml arguments -------------------------------------------------####
print("Reading cml arguments")

option_list <- list(
  optparse::make_option(c("-d", "--homedir"), type="character", 
                        default=file.path('~','git','lognorm_vs_CODA'), 
                        help="dataset dir path"),
  optparse::make_option(c("-p", "--project"), type="character", default="Zeller", 
                        help="project folder")
);

####Defining functions -----------------------------------------------------####
print("Defining functions")
source(file.path(home_dir, "lib", "table_manipulations.R"))

opt_parser <- optparse::OptionParser(option_list=option_list);

opt <- optparse::parse_args(opt_parser);

print(opt)

#Establishing directory layout and other constants--------------------------
print("Establishing directory layout and other constants.")

home_dir <- opt$homedir
project <- opt$project
working_dir <- file.path(home_dir, "lib", "unit_tests")
base::setwd(working_dir)

####Set up constants for lognorm test-----------------------------------------####
print("Setting up constants for testing lognorm")

if (!file.exists("ForwardReads_DADA2.txt")) unzip("lognorm.zip")

untransformed <- data.frame(data.table::fread(file = file.path("ForwardReads_DADA2.txt"),
                                           header=TRUE, data.table=FALSE), row.names = 1)
# untransformed <- read.table(file = file.path("logNorm","ForwardReads_DADA2.txt"),
#                                               header=TRUE, data.table=FALSE), row.names = 1)
# confirmed_transformed <- read.csv(file = file.path("logNorm","ForwardReads_DADA2LogNorm.txt.zip"),
#                                               header=TRUE, row.names = 1, sep = "\t")
# confirmed_transformed was transformed by Dr. Fodor's java program
confirmed_transformed <- data.frame(data.table::fread(file = file.path("ForwardReads_DADA2LogNorm.txt"),
                                                      header=TRUE, data.table=FALSE), row.names = 1)
####Testing lognorm-----------------------------------------------------------####
test_that("Testing lognorm",
          { expect_equal(lognorm(untransformed), confirmed_transformed)
            expect_equivalent(lognorm(untransformed), confirmed_transformed)
            # expect_identical(lognorm(untransformed), confirmed_transformed)
          })

#### Testing proportions and heilinger transformation ####
small_untransformed <- data.frame(matrix(data = c(0,1,2,3,4,4,3,2,1,0), nrow = 2, ncol = 5))

small_confirmed_proportions <- data.frame(
  matrix(data = c(0/10,1/10,2/10,3/10,4/10,4/10,3/10,2/10,1/10,0/10), 
         nrow = 2, ncol = 5))
small_confirmed_heilinger <- data.frame(
  matrix(data = c(sqrt(0/10),sqrt(1/10),sqrt(2/10),sqrt(3/10),sqrt(4/10),sqrt(4/10),sqrt(3/10),sqrt(2/10),sqrt(1/10),sqrt(0/10)), 
         nrow = 2, ncol = 5))

test_that("Testing proportions",
          { expect_equal(proportions_transform(small_untransformed), small_confirmed_proportions)
            expect_equivalent(proportions_transform(small_confirmed_proportions), small_confirmed_proportions)
            # expect_identical(lognorm(untransformed), confirmed_transformed)
          })
test_that("Testing Heilinger",
          { expect_equal(heilinger_transform(small_untransformed), small_confirmed_heilinger)
            expect_equivalent(heilinger_transform(small_untransformed), small_confirmed_heilinger)
            # expect_identical(lognorm(untransformed), confirmed_transformed)
          })
print("Completed script!")

