#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# IQ-TREE needs to be included in $PATH (v2.2.0; https://iqtree.github.io/)

## Command-line args:
in_file=$1
prefix=$2
ufboot=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### iqtree.sh: Starting script."
echo -e "#### iqtree.sh: Input file: $in_file"
echo -e "#### iqtree.sh: Prefix: $prefix"
echo -e "#### iqtree.sh: Number of ultrafast bootstrap replicates: $ufboot \n\n"

################################################################################################################
#### INFER ML PHYLOGENY WITH ASCERTAINMENT BIAS CORRECTION AND APPROXIMATE LIKELIHOOD RATIO TEST IN IQ-TREE ####
################################################################################################################
echo -e "#### iqtree.sh: Maximum likelihood phylogenetic inference ...\n"
iqtree2 -T AUTO -s $in_file --seqtype DNA -m GTR+G+ASC -nstop 200 -B $ufboot -wbt -bnni -alrt 1000 --prefix $prefix

## Report:
echo -e "\n#### iqtree.sh: Done with script."
date
