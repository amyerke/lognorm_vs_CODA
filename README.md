# Read depth and balance trees

This repository holds the code for the manuscript "Log-normalizing to read depth outperforms compositional data transformations in machine learning applications", in which in which we compared and contrasted different data transformations across different datasets. It was designed to be extensible and ensure that each transformation was treated fairly.

## Repository structure
### Common command line scripts
In order to ensure that each dataset is handled in the same way, common scripts that were run at the commandline were created. These are stored at [lib/cml_scripts](https://github.com/amyerke/lognorm_vs_CODA/tree/main/lib/cml_scripts). The scripts at this level are organized into the following folders based on similar function:
- [creat_ref_tree_blst_db](https://github.com/amyerke/lognorm_vs_CODA/tree/main/lib/cml_scripts/creat_ref_tree_blst_db) holds scripts for building the blast database that is used to build a phylogenetic tree for the Silva/LTP PhILR transformation.
- [data_preprocessing](https://github.com/amyerke/lognorm_vs_CODA/tree/main/lib/cml_scripts/data_preprocessing) holds scripts for running the datasets through DADA2.
- [make_ref_tree](https://github.com/amyerke/lognorm_vs_CODA/tree/main/lib/cml_scripts/make_ref_tree) uses the database created by the scripts in creat_ref_tree_blst_db to build trees for PhILR input.
- [metastudy_comparisons](https://github.com/amyerke/lognorm_vs_CODA/tree/main/lib/cml_scripts/metastudy_comparisons) houses scripts that use the data from the random forests to compare the output of all datasets.
- [mlm_metrics](https://github.com/amyerke/lognorm_vs_CODA/tree/main/lib/cml_scripts/mlm_metrics) holds scripts for comparing machine learning algorthims (MLAs).
- [random_forest](https://github.com/amyerke/lognorm_vs_CODA/tree/main/lib/cml_scripts/random_forest) holds scripts for running random forest.
- [transformations](https://github.com/amyerke/lognorm_vs_CODA/tree/main/lib/cml_scripts/transformations) holds scripts for building transformations other than the reference tree.

### Lib directory layout:
	lib
	├── CODA_demo
	├── cml_scripts
	│   ├── creat_ref_tree_blst_db
	│   ├── data_preprocessing
	│   ├── make_ref_tree
	│   ├── metastudy_comparisons
	│   ├── mlm_metrics
	│   ├── random_forest
	│   └── transformations
	├── one_off_scripts
	├── readme_images
	├── ref_tree_objs
	│   └── silva
	├── taxonomy
	└── unit_tests

### Project directory layout:
	[project]
	├── output
	│   ├── graphics
	│   ├── r_objects
	│   ├── tables
	│   ├── taxonomy
	│   ├── tree_process_blast
	│   └── trees
	└── scripts
			├── data_preprocessing
			├── download_scripts
			├── make_ref_tree
			├── mlm_weighting_py
			├── random_forest
			├── read_depth
			└── transformations


Many of the scripts used in this project assume that this repository is located at ~/git/lognorm_vs_CODA, where ~ represents the computer's home directory.

### Dataset/project scripts
Each project has it's own "scripts" subdirectory at the top level that closely matches the common command line scripts. This is where slurm and R scripts are located that launch jobs using common command line scripts. The ultimate goal of these scripts is to build the transformations and find the accuracy of each one with the random forest.
Most of the project scripts are wrapped in either slurm batch script or an R script that makes calls to the common scripts.

### Datasets used in this project:
All the datasets in the project are publicly available 16S datasets. Each dataset is given a folder at the root level that holds the data, shell scripts, and the output that is unique to that dataset - output that uses all the datasets collectively is house elsewhere (see below):
- [Noguera-Julian](https://github.com/amyerke/lognorm_vs_CODA/tree/main/Noguera-Julian) 
  * SRA Accession: [PRJNA307231](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA307231)
- [Vangay](https://github.com/amyerke/lognorm_vs_CODA/tree/main/Vangay)
  * SRA Accession: [PRJEB28687](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB28687)
- [Jones](https://github.com/amyerke/lognorm_vs_CODA/tree/main/Jones)
  * SRA Accession: [PRJNA397450](https://www.ebi.ac.uk/ena/browser/view/PRJNA397450)
- [Zeller](https://github.com/amyerke/lognorm_vs_CODA/tree/main/Zeller)
  * SRA Accession: [PRJNA397450](https://www.ebi.ac.uk/ena/browser/view/PRJEB6070)
- [metastudies]

### Output of the scripts
The various scripts used in this repository output various types of files. Each dataset/project in this repository has its own output directory at the top level.
- **graphics**: Contains PDFs, png, and other graphic output
- **r_objects**: Contains r objects, which often load faster than tables into R.
- **tables**: comma or tab seperated data tables
- **taxonomy**: files for building trees for PhILR input
- **tree_process_blast**: files for building Silva/LTP trees for PhILR input

The exception is that for output the metastudies is all in the metastudies/output folder and not seperated.

## Processing the sequencing data
In order to create downstream figures, the sequcencing data from each project must be processed and run through the random forest so that the accuracies can be recorded.
### Downloading sequencing data
 Each project contains a comma seperated file that holds the SRA accessions called "SraRunTable.txt" for that projects sequences. In each project's "scripts/download_scripts" folder there is an R script for downloading the fastQ files. This can be run with:

`Rscript sra_download.R`

- **Output**: fastq files in [project]/downloaded_seqs/

For this repository, only the fastQ files of the forward reads were used. 

### Running DADA2
#### Step 0: Determining truncation length
In order to determine where to truncate the sequences: 
- **common script**: [lib/cml_scripts/data_preprocessing/p0_dada2_find_trunLen.R](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/data_preprocessing/p0_dada2_find_trunLen.R)
- **project specific script**: `[project]/scripts/data_preprocessing/run_p0_dada2_find_trunLen.slurm`
- **Output**: [project]/output/graphics/plotQualF.png

#### Step 1: Running DADA2
Truncation as determined in previous script
- **common script**: [lib/cml_scripts/data_preprocessing/p1_dada2_rd1.R](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/data_preprocessing/p1_dada2_rd1.R)
- **project specific script**: `[project]/scripts/data_preprocessing/run_p0_dada2_find_trunLen.slurm`
- **Outputs**:
	- [project]/output/r_objects/ForwardReads_DADA2.rds
	- [project]/output/tables/ForwardReads_DADA2.txt
	- [project]/output/r_objects/ForwardReads_DADA2_taxonomy.rds
	- [project]/output/tables/ForwardReads_DADA2_taxonomy.aln
	- [project]/output/r_objects/ForwardReads_DADA2_alignment.rds
	- filtered fastQ files in [project]/downloaded_seqs/filtered

## Transformations used in this project
### Non-Compositional transformations
- Proportions: $$\frac{x_i}{\sum{x_{ij}}}$$
- Heilinger: $$\sqrt{\frac{x_i}{\sum{x_{ij}}}}$$
- lognorm: $$log(\frac{x_i}{\sum{x_{ij}}}*\frac{n}{N} + PC)$$
Where RC = raw counts in a cell, Σx = number of sequences in a sample, n = total number of counts in the table, N = total number of samples, PC = pseudo-count, for this project was equal to 1.
- Rarefaction to read depth of 1000 
<!-- log_10⁡((RC/n)×(∑(x/N)+PC))t -->
### Compositional transformations
- Additive log ratio (alr)
- Centered log ratio (clr)
- Phylogentic isometric log ratio
	- Made from:
		- UPGMA tree,
		- IQTREE2
		- Silva/Living Tree Project

### The following transformations are made using the same script:
- Additive log ratio (alr)
- Centered log ratio (clr)
- Lognorm
- Proportions
- Heilinger
- Rarefaction
- Prevalence filtered counts table

Each of these were created from with two counts tables that had each undergone different filtering methods:
1. No filter (the raw DADA2 output)
2. A prevalence filter where sequence varients in less than 10% were removed (prev_filt).

- **common script**: [lib/cml_scripts/transformations/make_nonilr_datasets.R](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/transformations/make_nonilr_datasets.R)
- **project specific script**: `[project]/scripts/transformations/run_make_nonilr_dataset.slurm`
- **Output**: 
	- [project]/output/r_objects/lognorm_asv.rds
	- [project]/output/tables/lognorm_dada2.csv
	- [project]/output/r_objects/alr_asv.rds
	- [project]/output/tables/alr_asv.csv
	- [project]/output/r_objects/clr_asv.rds
	- [project]/output/tables/clr_asv.csv
	- [project]/output/r_objects/propotions_asv.rds
	- [project]/output/tables/propotions_asv.csv
	- [project]/output/r_objects/heilinger_asv.rds
	- [project]/output/tables/heilinger_asv.csv
	- [project]/output/r_objects/1000_rrarefy_asv.rds
	- [project]/output/tables/1000_rrarefy_asv.csv
	- [project]/output/r_objects/filtered_90prcnt_dada2.rds
	- [project]/output/tables/filtered_90prcnt_dada2.csv
	- [project]/output/r_objects/filtered_90prcnt_lognorm_dada2.rds
	- [project]/output/tables/filtered_90prcnt_lognorm_dada2.csv
	- [project]/output/r_objects/filtered_90prcnt_alr_asv.rds
	- [project]/output/tables/filtered_90prcnt_alr_asv.csv
	- [project]/output/r_objects/filtered_90prcnt_clr_asv.rds
	- [project]/output/tables/filtered_90prcnt_clr_asv.csv
	- [project]/output/r_objects/filtered_90prcnt_propotions_asv.rds
	- [project]/output/tables/filtered_90prcnt_propotions_asv.csv
	- [project]/output/r_objects/filtered_90prcnt_heilinger_asv.rds
	- [project]/output/tables/filtered_90prcnt_heilinger_asv.csv
	- [project]/output/r_objects/filtered_90prcnt_1000_rrarefy_asv.rds
	- [project]/output/tables/filtered_90prcnt_1000_rrarefy_asv.csv

## Building Silva/Living Tree Project trees
The 16S reference tree that was used in this project can be found in [lib/ref_tree_objs/silva/](https://github.com/amyerke/lognorm_vs_CODA/blob/main/lib/ref_tree_objs/silva).

### Building the BLAST database
**The BLAST only has to be created once for the entire repository.** All the projects can use the same reference tree and therefore the same BLAST database.

#### Step 1: Downloading sequences from the ref tree
The reference tree has GenBank accession numbers for tip labels, but our counts tables use sequence variants, so the first step is to download the sequences. 

- **common script**: [lib/cml_scripts/creat_ref_tree_blst_db/p1_download_seqs.R](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/creat_ref_tree_blst_db/p1_download_seqs.R)
- **project specific script**: **N/A**
- **Output**: 
	- lib/ref_tree_objs/treeFasta.fasta

#### Step 2: Downloading sequences from the ref tree
This step requires BLAST.

- **common script**: [lib/cml_scripts/creat_ref_tree_blst_db/p2_makeDB.sh](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/creat_ref_tree_blst_db/p2_makeDB.sh)
- **project specific script**: **N/A**
- **Output**: 
	- blast database in lib/ref_tree_objs/db with the following files:
		- tree.nhr
		- tree.nin
		- tree.nog
		- tree.nsd 
		- tree.nsi
		- tree.nsq

### Create tree that is the union between project trees and the ref tree using BLAST

This needs to be done in each project for each counts table from which a Silva/LTP reference tree is made.

#### Reference tree summary figure
![Silva steps](lib/readme_images/ref_tree_steps.jpg)

#### Step 1: Convert project unique sequences from sequence variant counts table to fasta
In order find the union, the sequences will be BLASTed against the BLAST database that was made from the reference tree. 
The first step is to make proper input for the BLAST database in the form of a fasta file.

- **common script**: [lib/cml_scripts/make_ref_tree/p1_make_asv_fasta.R](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/make_ref_tree/p1_make_asv_fasta.R)
- **project specific scripts**:
	- `[project]/scripts/make_ref_tree/run_p1_make_asv_fasta_90prcent_filtered.slurm`
	- `[project]/scripts/make_ref_tree/run_p1_make_asv_fasta.slurm`
- **Output**: 
	fasta file with name designated in input - will be found in [project]output/tree_process_blast, ex:
	[project]output/tree_process_blast/dada2seqs.fasta

#### Step 2: Blast project fasta against ref tree database

- **common script**: [lib/cml_scripts/make_ref_tree/p2_blast.sh
](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/make_ref_tree/p2_blast.sh)
- **project specific scripts**:
	- `[project]/scripts/make_ref_tree/run_p2_blast_90prcnt_filt.slurm`
	- `[project]/scripts/make_ref_tree/run_p2_blast.slurm`
- **Output**: BLAST outout
	Found in [project]output/tree_process_blast/

#### Step 3: Parse BLAST output and build tree
- **common script**: [lib/cml_scripts/make_ref_tree/p3_parse_blast.R.sh
](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/make_ref_tree/p3_parse_blast.R)
- **project specific scripts**:
	- `[project]/scripts/make_ref_tree/run_p3_parse_blast_90prcnt_filt.slurm
	- [project]/scripts/make_ref_tree/run_p3_parse_blast.R`
- **Output**: `parsed BLAST outout
	Found in [project]output/tree_process_blast/
	ex: parsed_output.csv
`
#### Step 4: Create tree and phyloseq object from parsed output

- **common script**: [lib/cml_scripts/make_ref_tree/p4_ps_w_ref_tree.R
](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/make_ref_tree/p4_ps_w_ref_tree.R)
- **project specific scripts**:
	- `[project]/scripts/make_ref_tree/run_p4_ps_w_ref_tree_90prcnt_filt.slurm`
	- `[project]/scripts/make_ref_tree/run_p4_ps_w_ref_tree.slurm`
- **Output**: phyloseq object with reference tree

## For building the IQTREE trees:
  IQTREE requires aligned sequences. We experimented with both ClustalW and the DECIPHER R library, though ultimately decided to use DECIPHER because it can just be added as a line to the DADA2 script. However, the follow P0 will create an alignment for anyone who used an older version of p1_dada2_rd1.R.
### ~~P0 ClustalW~~
  ~~Using input from the build fasta step from~~ `p1_make_asv_fasta.R`
`/path/to/clustalo -in input.fasta -out output.fasta -fasta`
### Step 0 Make DECIEPHER Alignment*
- **common script**: [lib/cml_scripts/transformations/make_alignment.R](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/transformations/make_alignment.R)
- **project specific scripts**:
	- `[project]/scripts/make_ref_tree/Noguera-Julian/scripts/transformations/run_make_90prcent_alignment.slurm`
	- `[project]/scripts/make_ref_tree/Noguera-Julian/scripts/transformations/run_make_alignment.slurm`
- **Output**: Trees found in output/r_objects and output/tables. IE. ForwardReads_DADA2_alignment.rds and ForwardReads_DADA2_taxonomy.aln

### Step 1 use IQTREE2 to create tree
- **common script**: [lib/cml_scripts/transformations/make_alignment.R](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/transformations/make_alignment.R)
- **project specific scripts**:
	- `[project]/scripts/transformations/run_p5_make_iqtree.slurm`
	- `[project]/scripts/run_p5_make_90prcnt_filtered_iqtree.slurm`
- **Output**: Trees found in output/r_objects and output/tables. IE. ForwardReads_DADA2_alignment.rds and ForwardReads_DADA2_taxonomy.aln

#### Note on command for making IQTREE
`iqtree -s example.phy -T AUTO`
***For the Jones dataset, the GTR+F+R5 was selected because it had the lowest BIC after the modelfinder timing out at 48 hours. This was used for every other dataset afterwards.
```
No. 	Model         	-LnL         	df  AIC          	AICc       	 BIC
8	GTR+F+R5      	1688897.219	56065 3489924.438	6290170504.438	3806380.010
```
## Build UPGMA tree

### Step 1: Build UPGMA tree in R

- **common script**: [lib/cml_scripts/transformations/p2_denovo_tree_UPGMA.R](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/transformations/p2_denovo_tree_UPGMA.R)
- **project specific scripts**:
	- `[project]/scripts/transformations/run_p2_denovo_tree_UPGMA.slurm`
	- `[project]/scripts/transformations/run_p2_denovo_tree_UPGMA_90prevFilt.slurm`
- **Output**: 
 - phyloseq object in [project]/output/r_objects

## Create PhILR transformations with each tree

This script creates PhILR transformations from each tree, as well as shuffles the nodes of each tree and repeats the PhILR transformation.

- **common script**: [lib/cml_scripts/transformations/PhILR_random_trees_and_counts_tables.R](https://github.com/amyerke/lognorm_vs_CODA/lib/cml_scripts/transformations/PhILR_random_trees_and_counts_tables.R)
- **project specific script**:
	- `[project]/scripts/transformations/run_PhILR_random_trees_and_counts_tables.R.slurm`
- **Outputs**: 
	In [project]/output/tables a folder was made for each tree: 
	- IQTREE2 From Prevalence Filtered Counts table,
	- IQTREE2 From Unfiltered Counts Table,
	- UPGMA From Prevalence Filtered Counts table,
	- UPGMA From Unfiltered Counts Table,
	- Silva/LTP From Prevalence Filtered Counts table,
	- Silva/LTP From Unfiltered Counts table,

	In each folder was created:
	- The untransformed counts table from which the PhILR data was transformed,
	-	PhILR transformed data from the "true" tree,
	- 3 (or any amount selected) PhILR transforms with shuffled data,

	For Each PhILR transform, each combination of branch weighting schemes (uniform, blw, blw.sqrt, mean.descendants), AND part weighting scheme (uniform, gm.counts, anorm, anorm.x.gm.counts, enorm, and enorm.x.gm.counts). **Except the Prevalence filtered trees, which were only filtered with enorm and blw.sqrt because this was already determined to be the weighting scheme that we would use for downstream analysis and for some trees, PhILR was unable to perform all the transforms.**

## All common scripts for transformations
![Transformations](lib/readme_images/transfomrations_map.jpg)
Blue boxes are different transformations that are tested.

## Notes about transformations:
- For the Jones, Noguera, and Vangay datasets, the prevenalance filter, Silva filter, and the PhILR tutorial filtering removed all the sequences when used together.
- For Jones, Vangay, and Noguera-Julian, not all PhILR weighting schemes worked with prevalence filtering alone as the tree was too big for ape's node distance function, which is internally used by PhILR.
- For Jones dataset UPGMA, prevelance filtering and PhILR tutorial filtering removed all the sequences.

## Steps for Accuracy Metastudies

### 1. Random forest accuracy for each dataset
#### For each dataset:
```
lib/cml_scripts/random_forest/random_forest_manual_train.py
```
   This requires iqtree, upgma, silva, alr, clr, and lognorm transformations in the project out/tables/ directory. 
#### Example command: 
```
python ~/git/lognorm_vs_CODA/lib/cml_scripts/random_forest/random_forest_manual_train.py \
  --project Jones \
  --use_all_meta True \
  --metada_fn patient_metadata.tsv \
  --meta_index_col Run \
  --training 0.75
```
#### Output
* filename: sklearn_random_forest_manual_0.75train.csv
* found in: [project]/output/tables
* Provides: accuracy scores

### 2. Metastudy accuracy vs accuracy
#### Combines data from all studies
* requires: sklearn_random_forest_manual_0.75train.csv in output/tables
* Compares each transformation to each other one.
```
lib/cml_scripts/random_forest/random_forest_manual_train.py
```



