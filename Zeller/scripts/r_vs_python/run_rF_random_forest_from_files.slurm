#!/bin/bash

#SBATCH --partition=Pisces                                    		# Partition/queue requested on server
#SBATCH --job-name=Zeller_rF_random_forest_from_file                                   		# Job name
#SBATCH --time=200:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=150gb                                          		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/%x.%j.log     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/%x.%j.log      # Set directory for error log
##SBATCH --uid=amyerke
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load R

Rscript ~/git/lognorm_vs_CODA/lib/cml_scripts/r_vs_python/rF_random_forest_from_files.R \
  -d ~/git/lognorm_vs_CODA \
  -p Zeller 

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

