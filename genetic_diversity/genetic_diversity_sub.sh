###################################
#### HETEROZYGOSITY ESTIMATION ####
###################################
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)


################
## Microcebus ##
################
vcf_file=/scratch/projects/nib00015/northeastProject/microcebus/all/stacks/vcfFiltering/populations.snps.filtered06.mac3.vcf
out_dir=/scratch/usr/nibtve93/northeastProject/avahi/all/heterozygosity
mkdir -p $out_dir
vcftools --vcf $vcf_file --het --out $out_dir/microcebus_GC_based


###########
## Avahi ##
###########
vcf_file=/scratch/usr/nibtve93/northeastProject/avahi/all/stacks/final/vcfFiltering/populations.snps.filtered06.mac3.vcf
out_dir=/scratch/usr/nibtve93/northeastProject/avahi/all/heterozygosity
mkdir -p $out_dir
vcftools --vcf $vcf_file --het --out $out_dir/avahi_GC_based

