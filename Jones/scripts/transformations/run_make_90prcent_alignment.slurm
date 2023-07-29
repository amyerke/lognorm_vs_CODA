#!/bin/bash

#SBATCH --partition=Orion                              		# Partition/queue requested on server
#SBATCH --job-name=run_make_90_filtered_alignment_Jones		# Job name
#SBATCH --time=24:00:00                                		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                    			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=50gb                                     		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                       	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                     	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/%x.%j.out       # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/%x.%j.out        # Set directory for error log
##SBATCH --uid=amyerke
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load R

Rscript ~/git/lognorm_vs_CODA/lib/cml_scripts/transformations/make_alignment.R
--homedir ~/git/lognorm_vs_CODA \
--project Jones \
--input_table filtered_90prcnt_dada2.rds \
--output_file filtered_90prcnt_aln

# Options:
#         -d CHARACTER, --homedir=CHARACTER
#                 dataset dir path
# 
#         -p CHARACTER, --project=CHARACTER
#                 project folder
# 
#         -i INITIAL_TABLE, --initial_table=INITIAL_TABLE
#                 initial table, relative to project path, with filename, 
#                         should be an r object
# 
#         -o OUTPUT_FILE, --output_file=OUTPUT_FILE
#                 output file name - will be in output/tables and output/r_objects, don't add extension
# 
#         -h, --help
#                 Show this help message and exit

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
