#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# ADMIXTURE needs to be included in $PATH (v1.3.0; https://dalexander.github.io/admixture/)

## Command-line args:
in_file=$1
K=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### clustering.sh: Starting script."
echo -e "#### clustering.sh: Input file: $in_file"
echo -e "#### clustering.sh: Number of clusters: $K \n\n"

##############################################################################
#### CLUSTERING ANALYSIS WITH ADMIXTURE ####
##############################################################################
echo -e "#### clustering.sh: Running clustering analysis ...\n"
cd $(dirname $in_file)
admixture --cv $in_file $K

## Report:
echo -e "\n#### clustering.sh: Done with script."
date
