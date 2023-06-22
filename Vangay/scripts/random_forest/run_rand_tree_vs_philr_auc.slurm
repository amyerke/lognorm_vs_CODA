#!/bin/bash

#SBATCH --partition=Orion                                    		# Partition/queue requested on server
#SBATCH --job-name=ilr_vs_philr_slurm                                   		# Job name
#SBATCH --time=50:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=50gb                                          		# Memory limit
#SBATCH --mail-type=ALL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/Vangay_%x.%j.log     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/Vangay_%x.%j.log      # Set directory for error log
##SBATCH --uid=amyerke
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load R

Rscript ~/git/lognorm_vs_CODA/lib/cml_scripts/random_forest/random_forest_philr_comparison.R \
  -d ~/git/lognorm_vs_CODA \
  -p Vangay \
  -m ~/git/lognorm_vs_CODA/Vangay/patient_metadata.tsv \
  -r run_accession \
  -n 20

  # optparse::make_option(c("-d", "--homedir"), type="character", 
  #                       default=file.path('~','git',"lognorm_vs_CODA"), 
  #                       help="dataset dir path", metavar="character"),
  # optparse::make_option(c("-p", "--project"), type="character", default=NULL, 
  #                       help="project folder", metavar="character"),
  # optparse::make_option(c("-m", "--metadata"), type="character", default=NULL,
  #                       help="metadata file path with filename", metavar="character"),
  # optparse::make_option(c("-l", "--metadata_delim"), type="character", default="\t",
  #                       help="metadata file deliminator", metavar="character"),
  # optparse::make_option(c("-r", "--metadata_rowname"), type="character", default=NULL,
  #                       help="metadata file row to use for row names", metavar="character"),
  # optparse::make_option(c("-n", "--num_cycles"), type="numeric", default=20,
  #                       help="Number of times to shuffle data and run loop again", 
  #                       metavar="character")

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"


