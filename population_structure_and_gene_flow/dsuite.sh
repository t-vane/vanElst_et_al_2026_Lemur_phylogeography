#!/bin/bash
#SBATCH -p medium96s

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# Dsuite needs to be included in $PATH (v0.4 r38; https://github.com/millanek/Dsuite)

## Command-line args:
tree=$1
out=$2
vcf=$3
sets=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### dsuite.sh: Starting script."
echo -e "#### dsuite.sh: Input tree file: $tree"
echo -e "#### dsuite.sh: Output prefix: $out"
echo -e "#### dsuite.sh: VCF file: $vcf"
echo -e "#### dsuite.sh: Taxon sets: $sets \n\n"

#############################################################################
#### TEST FOR EXCESS ALLELE SHARING ####
#############################################################################
echo -e "#### dsuite.sh: Testing for excess allele sharing ...\n"
Dsuite Dtrios -t $tree -o $out $vcf $sets

## Report:
echo -e "\n#### dsuite.sh: Done with script."
date


