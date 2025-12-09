#!/bin/sh
#SBATCH -p standard96

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
ascbias=/home/nibtve93/software/raxml_ascbias/ascbias.py # (https://github.com/btmartin721/raxml_ascbias)
amas=/home/nibtve93/software/AMAS/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)
tabtofasta=/home/nibtve93/software/vcf_tab_to_fasta_alignment_TVE.pl
export PERL5LIB=/home/nibtve93/software/vcftools/src/perl/

## Command-line args:
in_file=$1
working_dir=$2
prefix=$(basename $in_file .vcf)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### prepareAlignment.sh: Starting script."
echo -e "#### prepareAlignment.sh: Input file: $in_file"
echo -e "#### prepareAlignment.sh: Working directory: $working_dir"
echo -e "#### prepareAlignment.sh: Prefix: $prefix \n\n"

################################################################################
#### CREATING ALIGNMENT FROM VCF FILE ####
################################################################################
echo "#### prepareAlignment.sh: Creating tab from VCF file"
vcf-to-tab < $in_file > $working_dir/$prefix.tab

echo "#### prepareAlignment.sh: Converting tab file to FASTA format"
perl $tabtofasta -i $working_dir/$prefix.tab > $working_dir/$prefix.fasta

echo "#### prepareAlignment.sh: Converting FASTA to PHYLIP format"
python $amas convert -i $working_dir/$prefix.fasta -f fasta -u phylip -d dna
mv $working_dir/$prefix.fasta-out.phy $working_dir/$prefix.phy

echo "#### prepareAlignment.sh: Removing invariant sites"
python $ascbias -p $working_dir/$prefix.phy -o $working_dir/$prefix.noinv.phy

echo "#### prepareAlignment.sh: Calculating alignment summary"
python $amas summary -f phylip -d dna -i $working_dir/$prefix.noinv.phy -o $working_dir/$prefix.noinv.phy.summary

## Report:
echo -e "#### prepareAlignment.sh: Done with script."
date
