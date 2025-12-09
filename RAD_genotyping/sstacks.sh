#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# sstacks of Stacks needs to be included in $PATH (v2.53; http://catchenlab.life.illinois.edu/stacks/)

## Command-line args:
catalog_dir=$1
popmap=$2
nt=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### sstacks_denovo.sh: Starting script."
echo -e "#### sstacks_denovo.sh: Catalog directory: $catalog_dir"
echo -e "#### sstacks_denovo.sh: Population map: $popmap"
echo -e "#### sstacks_denovo.sh: Number of threads: $nt \n\n"

################################################################################
#### MATCH SAMPLES AGAINST CATALOG ####
################################################################################
echo -e "#### sstacks_denovo.sh: Matching samples against catalog with sstacks ...\n"
sstacks -P $catalog_dir -M $popmap -p $nt

## Report:
echo -e "\n#### sstacks_denovo.sh: Done with script."
date







