################################
#### PHYLOGENETIC INFERENCE ####
################################
script_dir=/home/nibtve93/scripts/phylogeneticInference


################
## Microcebus ##
################
phyl_dir=/scratch/projects/nib00015/northeastProject/microcebus/all/phylogeneticInference
vcf_dir=/scratch/projects/nib00015/northeastProject/microcebus/all/stacks/vcfFiltering

mkdir -p $phyl_dir/logFiles $phyl_dir/iqtree

## Prepare for alignment
sbatch --output=$phyl_dir/logFiles/prepareForRaxml.oe $scripts_dir/prepare_alignment.sh $vcf_dir/populations.snps.filtered06.mac3.vcf $phyl_dir

## Run IQ-TREE
ufboot=1000
sbatch --output=$phyl_dir/logFiles/iqtree.oe --account=nib00015 $scripts_dir/iqtree_checkpoint.sh $phyl_dir/populations.snps.filtered06.mac3.noinv.phy $phyl_dir/iqtree/populations.snps.filtered06.mac3.noinv $ufboot


################
## Avahi ##
################
phyl_dir=/scratch/usr/nibtve93/northeastProject/avahi/all/phylogeneticInference
vcf_dir=/scratch/usr/nibtve93/northeastProject/avahi/all/stacks/final/vcfFiltering

mkdir -p $phyl_dir/logFiles $phyl_dir/iqtree

## Prepare for alignment
sbatch --output=$phyl_dir/logFiles/prepareForRaxml.oe $scripts_dir/prepare_alignment.sh $vcf_dir/populations.snps.filtered06.mac3.vcf $phyl_dir

## Run IQ-TREE
ufboot=1000
sbatch --output=$phyl_dir/logFiles/iqtree.oe --account=nib00015 $scripts_dir/iqtree_checkpoint.sh $phyl_dir/populations.snps.filtered06.mac3.noinv.phy $phyl_dir/iqtree/populations.snps.filtered06.mac3.noinv $ufboot
