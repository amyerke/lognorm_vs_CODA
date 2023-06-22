#!/bin/bash

#SBATCH --partition=Pisces                                    		# Partition/queue requested on server
#SBATCH --job-name=rf_all                                   		# Job name
#SBATCH --time=30:00:00                                      		# Time limit (hrs:min:sec)
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

export PYTHONUNBUFFERED=TRUE #lets print statements show up without flush=True
module load anaconda3/2020.11
conda activate testing_pyCaret

echo launching Jones

srun --nodes=1 python ~/git/lognorm_vs_CODA/lib/cml_scripts/random_forest/random_forest_manual_train_ntree_mtry.py\
  -p Jones \
  -a True \
  -f patient_metadata.tsv \
  -i Run \
  -t 0.75 \
  -n 500

echo launching Vangay

srun --nodes=1 python ~/git/lognorm_vs_CODA/lib/cml_scripts/random_forest/random_forest_manual_train_ntree_mtry.py\
  -p Vangay \
  -a True \
  -f patient_metadata.tsv \
  -i run_accession \
  -t 0.75 \
  -n 500

echo launching Zeller

srun --nodes=1 python ~/git/lognorm_vs_CODA/lib/cml_scripts/random_forest/random_forest_manual_train_ntree_mtry.py\
  -p Zeller \
  -a True \
  -f patient_metadata.csv \
  -i Run \
	-l , \
  -t 0.75 \
  -n 500

echo launching Noguera-Julian

srun --nodes=1 python ~/git/lognorm_vs_CODA/lib/cml_scripts/random_forest/random_forest_manual_train_ntree_mtry.py\
  -p Noguera-Julian \
  -a True \
  -f patient_metadata.tsv \
  -i Run \
  -t 0.75 \
  -n 500

conda deactivate
module unload anaconda3/2020.11

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

