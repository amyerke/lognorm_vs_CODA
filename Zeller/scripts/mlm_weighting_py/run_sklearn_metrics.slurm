#!/bin/bash

#SBATCH --partition=Pisces                           # Partition/queue requested on server
#SBATCH --job-name=model_metrics_sklearn_Zeller      # Job name
#SBATCH --time=2:40:00                               # Time limit (hrs:min:sec)
#SBATCH --nodes=1                                    # Number of nodes requested
#SBATCH --ntasks-per-node=1                          # Number of CPUs (processor cores/tasks)
#SBATCH --mem=70gb                                      # Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL                      # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=amyerke@uncc.edu                    # Specified email address
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
srun python ~/git/lognorm_vs_CODA/lib/cml_scripts/mlm_metrics/sklearn_model_metrics.py\
  -p Zeller \
  -a True \
  -f patient_metadata.csv \
  -l , \
  -i Run

conda deactivate
module unload anaconda3/2020.11

# usage: sklearn_model_metrics.py [-h] [-m] [-d homedir] [-p project]
#                                 [-a use_all_meta] [-f meta_fn] [-l delim]
#                                 [-i meta_index_col]

# Process some integers.

# optional arguments:
#   -h, --help            show this help message and exit
#   -m, --metadata_cols   Metadata columns to analyse
#   -d homedir, --homedir homedir
#                         dataset dir path
#   -p project, --project project
#                         project folder
#   -a use_all_meta, --use_all_meta use_all_meta
#                         use all metadata
#   -f meta_fn, --metada-fn meta_fn
#                         Name of file at the top of the project folder to use
#                         as metadata.
#   -l delim, --delimiter delim
#                         File delimiting symbol for metadata.
#   -i meta_index_col, --meta_index_col meta_index_col
#                         Name of column to use as row name for metadata

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

