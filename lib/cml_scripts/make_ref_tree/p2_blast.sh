#Author: Aaron Yerke (aaronyerke@gmail.com)
echo "Loading blast module"
module load blast/2.9.0+

home_dir=$1 #first comandline argument is the project name
project=$2
fasta=$3
output=$4
echo "Found arguments ${home_dir} and ${project}."

db_path=${home_dir}/lib/ref_tree_objs/db/tree

echo "Cd to ~/git/lognorm_vs_CODA/${project}/output/tree_process_blast."
cd ~/git/lognorm_vs_CODA/${project}/output/tree_process_blast
# cd ~/git/lognorm_vs_CODA/Fodor/output/tree_process_blast

echo "Attempting to blastn."
blastn -query ${fasta} \
  -db ${db_path} \
  -out ${output} \
  -outfmt "6 qseqid sseqid pident length evalue bitscore score ppos"

echo "Reached end of script"
 # 1.	 qseqid	 query (e.g., unknown gene) sequence id
 # 2.	 sseqid	 subject (e.g., reference genome) sequence id
 # 3.	 pident	 percentage of identical matches
 # 4.	 length	 alignment length (sequence overlap)
 # 5.	 evalue	 expect value
 # 6.	 bitscore	 bit score
 # 7.  score     Raw score
 # 8.  ppos      Percentage of positive-scoring matches