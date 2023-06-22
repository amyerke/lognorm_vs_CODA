#!/bin/bash

#SBATCH --partition=Pisces                                    		# Partition/queue requested on server
#SBATCH --job-name=inter_dataset_rand                                   		# Job name
#SBATCH --time=20:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=4                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=200gb                                          		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/%x.%j.log     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/%x.%j.log      # Set directory for error log
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load anaconda3/2020.11
conda activate testing_pyCaret
python ~/git/lognorm_vs_CODA/metastudies/scripts/compare_rand_forest/python_accuracy_v_accuracy.py\
  -d ~/git/lognorm_vs_CODA

conda deactivate
module unload anaconda3/2020.11

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

