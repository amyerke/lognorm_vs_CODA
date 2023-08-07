#!/bin/bash

#SBATCH --partition=Orion                                    		# Partition/queue requested on server
#SBATCH --job-name=run_p3_parse_blast                                   		# Job name
#SBATCH --time=24:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=50gb                                          		# Memory limit
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

Rscript ~/git/lognorm_vs_CODA/lib/cml_scripts/make_ref_tree/p3_parse_blast.R \
  -d ~/git/lognorm_vs_CODA \
  -p Noguera-Julian \
  -i output_90prcnt_filt.txt \
  -o parsed_output_90prcnt_filt.csv

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
