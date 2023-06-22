#!/bin/bash

#SBATCH --partition=Orion                                    		# Partition/queue requested on server
#SBATCH --job-name=testing_pyCaret                                   		# Job name
#SBATCH --time=25:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=50gb                                          		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/vand_%x.%j.out     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/vand_%x.%j.err      # Set directory for error log
##SBATCH --uid=amyerke
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

#module load anaconda3/2020.11
#conda activate testing_pyCaret

srun python ~/git/lognorm_vs_CODA/Jones/scripts/mlm_weighting_py/ref_tree_caret.py

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

