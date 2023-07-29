#Script from creating IQTree on cluster - requires CLUSTALW and IQTREE modules

#Getting fancy with color in print statements
RED='\033[0;31m' #for print red
NC='\033[0m' # No Color

home_dir=$1 #first comandline argument is the project name
project=$2
aligned=$3

echo "Found arguments ${home_dir}, ${project}, and ${aligned}."

#my_fasta is made from
if [ "$aligned" = ""]; then
    echo "aligned will default to"
    aligned=$home_dir/$project/output/trees/ForwardReads_DADA2_taxonomy.aln
    echo "$aligned"
fi

#check for fastq file that is our starting point
if [ -f "$aligned" ]; then
    echo "$aligned exists."
else 
    printf "${RED}$aligned does not exist,\nrun your sequences through Dada2 and create the fastq with p1 and p2${NC}"
fi

cd $home_dir/$project/output/trees/

module load iqtree

iqtree2 -s $aligned -T AUTO -m GTR+F+R5 

printf "${RED}Reached end of script.${NC}"

