#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# ustacks of Stacks needs to be included in $PATH (v2.53; http://catchenlab.life.illinois.edu/stacks/)

## Command-line args:
read_dir=$1
in_file=$2
indv=$(sed -n "$SLURM_ARRAY_TASK_ID"p $in_file)
out_dir=$3
m=$4
M=$5
nt=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### ustacks_denovo.sh: Starting script."
echo -e "#### ustacks_denovo.sh: Read directory: $read_dir"
echo -e "#### ustacks_denovo.sh: Input file: $in_file"
echo -e "#### ustacks_denovo.sh: Individual: $indv"
echo -e "#### ustacks_denovo.sh: Output directory: $out_dir"
echo -e "#### ustacks_denovo.sh: Parameter value of m: $m"
echo -e "#### ustacks_denovo.sh: Parameter value of M: $M"
echo -e "#### ustacks_denovo.sh: Number of threads: $nt \n\n"

################################################################################
#### BUILD LOCI FOR SAMPLE ####
################################################################################
echo -e "#### ustacks_denovo.sh: Building loci for sample $indv with ustacks ...\n"
ustacks -f $read_dir/$indv.1.fq.gz -i $SLURM_ARRAY_TASK_ID -o $out_dir -m $m -M $M -p $NC -t gzfastq --name $indv

## Report:
echo -e "\n#### ustacks_denovo.sh: Done with script."
date







