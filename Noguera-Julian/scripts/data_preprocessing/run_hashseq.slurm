#!/bin/bash

#SBATCH --partition=Orion                                    		# Partition/queue requested on server
#SBATCH --job-name=hashseq                                   		# Job name
#SBATCH --time=5:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=500gb                                          		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/Noguera-Julian_%x.%j.out     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/Noguera-Julian_%x.%j.out      # Set directory for error log
##SBATCH --uid=amyerke
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load R

srun Rscript ~/git/lognorm_vs_CODA/lib/cml_scripts/data_preprocessing/hashseq.R \
	-i ~/git/lognorm_vs_CODA/Noguera-Julian/downloaded_seqs \
	-o ~/git/lognorm_vs_CODA/Noguera-Julian/output/hashseq

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
