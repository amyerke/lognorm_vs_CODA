# Read depth and balance trees

> This repository holds the code for my PhD dissertation, which is a comparison of balance tree methods across different datasets and a study on how read depth affects the downstream results.

## Datasets used in this project:
> All the datasets in the project are publicly available 16S datasets
* Noguera-Julian
  - Accession: [PRJNA307231](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA307231)
* Vangay
  - Accession: [PRJEB28687](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB28687)
* Jones
  - Accession: [PRJNA397450](https://www.ebi.ac.uk/ena/browser/view/PRJNA397450)
* Zeller
  - Accession: [PRJNA397450](https://www.ebi.ac.uk/ena/browser/view/PRJEB6070)

## Transformations used in this project

## Building Tree of Life Project trees

## For building the IQTREE trees:
  IQTREE requires aligned sequences, therefore we first must align the sequences with ClustalW.
### ClustalW
  Using input from the build fasta step from `p1_make_asv_fasta.R`
`/path/to/clustalo -in input.fasta -out output.fasta -fasta`

### Command for making IQTREE
`iqtree -s example.phy -T AUTO`
***For the Jones dataset, the GTR+F+R5 was selected because it had the lowest BIC after the modelfinder timing out at 48 hours.
```
No. 	Model         	-LnL         	df  AIC          	AICc       	 BIC
8	GTR+F+R5      	1688897.219	56065 3489924.438	6290170504.438	3806380.010
```
## Transformations
![Transformations](lib/readme_images/transfomrations_map.jpg)
Blue boxes are different transformations that are tested.

### 1. Silva LTP Reference Tree
![Silva steps](lib/readme_images/ref_tree_steps.jpg)
### 2. IQTREE2
### 3. UPGMA
### 3. Lognorm, ALR, and CLR

## Steps for Accuracy Metastudies

### 1. Random forest accuracy for each dataset
#### For each dataset:
```
lib/cml_scripts/random_forest/random_forest_manual_train.py
```
   This requires iqtree, upgma, silva, alr, clr, and lognorm transformations in the project out/tables/ directory. 
#### Example command: 
```
python ~/git/balance_tree_exploration/lib/cml_scripts/random_forest/random_forest_manual_train.py \
  --project Jones \
  --use_all_meta True \
  --metada_fn patient_metadata.tsv \
  --meta_index_col Run \
  --training 0.75
```
#### Output
* filename: sklearn_random_forest_manual_0.75train.csv
* found in: project/output/tables
* Provides: accuracy scores

### 2. Metastudy accuracy vs accuracy
#### Combines data from all studies
* requires: sklearn_random_forest_manual_0.75train.csv in output/tables
* Compares each transformation to each other one.
```
lib/cml_scripts/random_forest/random_forest_manual_train.py
```
