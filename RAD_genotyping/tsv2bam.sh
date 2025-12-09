#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# tsv2bam of Stacks needs to be included in $PATH (v2.53; http://catchenlab.life.illinois.edu/stacks/)

## Command-line args:
catalog_dir=$1
read_dir=$2
popmap=$2
nt=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### tsv2bam_denovo.sh: Starting script."
echo -e "#### tsv2bam_denovo.sh: Catalog directory: $catalog_dir"
echo -e "#### tsv2bam_denovo.sh: Read directory: $read_dir"
echo -e "#### tsv2bam_denovo.sh: Population map: $popmap"
echo -e "#### tsv2bam_denovo.sh: Number of threads: $nt \n\n"

################################################################################
#### TRANSPOSE DATA ####
################################################################################
echo -e "#### tsv2bam_denovo.sh: Transposing data with tsv2bam ...\n"
tsv2bam -P $catalog_dir -R $read_dir -M $popmap -t $nt

## Report:
echo -e "\n#### tsv2bam_denovo.sh: Done with script."
date







