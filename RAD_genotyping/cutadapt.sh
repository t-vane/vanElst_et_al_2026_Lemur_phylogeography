#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# cutadapt needs to be included in $PATH (v4.2; https://cutadapt.readthedocs.io/en/stable/)

## Command-line args:
length=$1
ind=$2
read_dir=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### cutadapt.sh: Starting script."
echo -e "#### cutadapt.sh: Target read length: $length"
echo -e "#### cutadapt.sh: Individual: $ind"
echo -e "#### cutadapt.sh: Read dir: $read_dir \n\n"

################################################################################
#### CUT READS TO UNIFORM LENGTH ####
################################################################################
echo -e "#### cutadapt.sh: Cutting reads to uniform length of $length ...\n"
gzip -d $read_dir/$ind.trimmed.1.fq.gz $read_dir/$ind.trimmed.2.fq.gz
cutadapt -l $length -m $length -o $read_dir/$ind.1.fq -p $read_dir/$ind.2.fq $read_dir/$ind.trimmed.1.fq $read_dir/$ind.trimmed.2.fq
gzip $read_dir/$ind.1.fq $read_dir/$ind.2.fq
gzip $read_dir/$ind.trimmed.1.fq $read_dir/$ind.trimmed.2.fq

## Report:
echo -e "\n#### cutadapt.sh: Done with script."
date







