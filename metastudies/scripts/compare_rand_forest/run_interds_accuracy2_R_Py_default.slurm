#!/bin/bash

#SBATCH --partition=Pisces                                    		# Partition/queue requested on server
#SBATCH --job-name=inter_dataset_rand                                   		# Job name
#SBATCH --time=2:00:00                                      		# Time limit (hrs:min:sec)
#SBATCH --nodes=1                                         			# Number of nodes requested
#SBATCH --ntasks-per-node=1                          			# Number of CPUs (processor cores/tasks)
#SBATCH --mem=50gb                                          		# Memory limit
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

module load anaconda3/2020.11
conda activate testing_pyCaret
python ~/git/balance_tree_exploration/lib/cml_scripts/metastudy_comparisons/inter_ds_accuracy_v_accuracy.py \
  --homedir ~/git/balance_tree_exploration \
  --x_axis_fname wide_random_forest_score_R_20.csv \
  --y_axis_fname sklearn_random_forest_manual_0.75train.csv \
  --x_axis_label R_randomForest \
  --y_axis_label Sklearn_RF_default \
  --x_axis_col_name trans_group \
  --y_axis_col_name dataset \
  --project_list Jones Zeller Noguera-Julian Vangay

conda deactivate
module unload anaconda3/2020.11

# parser.add_argument("-d", "--homedir",
#                   default=os.path.expanduser(os.path.join("~", "git", "lognorm_vs_CODA")),
#                   help="path to git balance treee exploration git repository", dest="homedir", metavar="homedir")
# parser.add_argument("-x", "--x_axis_fname", default=None,
#             			help="File name of the csv in the \{project\}/output/tables section in each project for the x axis")
# parser.add_argument("-y", "--y_axis_fname", default=None,
#             			help="File name of the csv in the \{project\}/output/tables section in each project for the y axis")
# parser.add_argument("-a", "--x_axis_col_name", help="column name from the x axis csv", default="trans_group")
# parser.add_argument("-b", "--y_axis_col_name", help="column name from the y axis csv", default="dataset")
# parser.add_argument("-t", "--x_axis_metadata_col_name", help="metadata column name from the x axis csv", default="metadata_col")
# parser.add_argument("-u", "--y_axis_metadata_col_name", help="metadata column name from the y axis csv", default="metadata")
# parser.add_argument('-p', '--project_list', nargs='+', default=[],
# 										help="Project folders that have the x and y files.")
# options, unknown = parser.parse_known_args()
# parser.add_argument("-v", "--x_axis_label", help="label for the x axis", default=os.path.basename(options.x_axis_fname))
# parser.add_argument("-w", "--y_axis_label", help="label for the x axis", default=os.path.basename(options.y_axis_fname))
# options, unknown = parser.parse_known_args()

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

