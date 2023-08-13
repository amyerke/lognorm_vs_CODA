#!/bin/bash

#SBATCH --partition=Pisces                                    		# Partition/queue requested on server
#SBATCH --job-name=Jones_rev_python_random_forest_manual   		# Job name
#SBATCH --time=48:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=100gb                                          		# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                              	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                            	# Specified email address
#SBATCH --output=/users/amyerke/slurmLogs/%x.%j.log     # Set directory for standard output
#SBATCH --error=/users/amyerke/slurmLogs/%x.%j.log      # Set directory for error log
##SBATCH --uid=amyerke
#SBATCH --get-user-env

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

export PYTHONUNBUFFERED=TRUE #lets print statements show up without flush=True

module load anaconda3/2020.11
conda activate testing_transforms_conda
python ~/git/lognorm_vs_CODA/lib/cml_scripts/random_forest/revision_random_forest_manual_train.py\
  --project Jones \
  --use_all_meta True \
  --metada_fn patient_metadata.tsv \
  --meta_index_col Run \
  --training 0.75
conda deactivate
module unload anaconda3/2020.11
# usage: random_forest_manual_train.py [-h] [-m] [-d homedir] [-p project]
#                                      [-a use_all_meta] [-f meta_fn] [-l delim]
#                                      [-i meta_index_col] [-t training]

# This script runs a random forest test on various datasets.

# optional arguments:
#   -h, --help            show this help message and exit
#   -m, --metadata_cols   Metadata columns to analyse
#   -d homedir, --homedir homedir
#                         path to git balance treee exploration git repository
#   -p project, --project project
#                         project folder
#   -a use_all_meta, --use_all_meta use_all_meta
#                         use all metadata
#   -f meta_fn, --metada_fn meta_fn
#                         Name of file at the top of the project folder to use
#                         as metadata.
#   -l delim, --delimiter delim
#                         File delimiting symbol for metadata. Default is tab.
#   -i meta_index_col, --meta_index_col meta_index_col
#                         Name of column to use as row name for metadata
#   -t training, --training training
#                         Percentating of table to use for training. The rest
#                         will be used for testing.
echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"



