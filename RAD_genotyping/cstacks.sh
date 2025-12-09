#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# cstacks of Stacks needs to be included in $PATH (v2.53; http://catchenlab.life.illinois.edu/stacks/)

## Command-line args:
stacks_dir=$1
catalog_dir=$2
popmap=$3
m=$4
M=$5
n=$6
nt=$7

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### cstacks_denovo.sh: Starting script."
echo -e "#### cstacks_denovo.sh: Stacks directory: $stacks_dir"
echo -e "#### cstacks_denovo.sh: Catalog directory: $catalog_dir"
echo -e "#### cstacks_denovo.sh: Population map: $popmap"
echo -e "#### cstacks_denovo.sh: Parameter value of n: $m"
echo -e "#### cstacks_denovo.sh: Parameter value of n: $M"
echo -e "#### cstacks_denovo.sh: Parameter value of n: $n"
echo -e "#### cstacks_denovo.sh: Number of threads: $nt \n\n"

################################################################################
#### BUILD CATALOG OF LOCI ####
################################################################################
echo -e "#### cstacks_denovo.sh: Building catalog of loci with cstacks ...\n"
ln -s $stacks_dir/stacks.m$m.M$M/*gz $catalog_dir
cstacks -P $catalog_dir -M $popmap -n $n -p $nt --report_mmatches

## Report:
echo -e "\n#### cstacks_denovo.sh: Done with script."
date







