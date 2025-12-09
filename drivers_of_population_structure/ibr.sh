#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args
scripts_dir=$1
wd=$2
prefix=$3
model=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### ibr.sh: Starting script."
echo -e "#### ibr.sh: Scripts directory: $scripts_dir"
echo -e "#### ibr.sh: Working directory: $wd"
echo -e "#### ibr.sh: Taxon set: $prefix"
echo -e "#### ibr.sh: Model: $model \n\n"

################################################################################
#### EVALUATE ISOLATION-BY-RESISTANCE MODEL ####
################################################################################
echo -e "#### ibr.sh: Evaluating isolation-by-resistance model ...\n"
Rscript $scripts_dir/ibr_radish.R $wd $prefix $model

## Report:
echo -e "\n#### ibr.sh: Done with script."
date