#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# gstacks of Stacks needs to be included in $PATH (v2.53; http://catchenlab.life.illinois.edu/stacks/)

## Command-line args:
main=$1
popmap=$2
nt=$3


## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### gstacks.sh: Starting script."
echo -e "#### gstacks.sh: Main parameters: $main"
echo -e "#### gstacks.sh: Population map: $popmap"
echo -e "#### gstacks.sh: Suffix for BAM files: $suffix \n\n"

################################################################################
#### CREATE STACKS WITH GSTACKS ####
################################################################################
echo -e "#### gstacks.sh: Creating stacks with gstacks for individuals in $popmap ...\n"
gstacks $main -M $popmap -t $nt


## Report:
echo -e "\n#### gstacks.sh: Done with script."
date
