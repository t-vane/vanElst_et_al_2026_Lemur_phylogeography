#############################
#### CLUSTERING ANALYSIS ####
#############################
## PLINK needs to be included in $PATH (v1.9; https://vcftools.github.io/index.html)

scripts_dir=/home/nibtve93/scripts/popGen
populations=14

################
## Microcebus ##
################
## Set up an run clustering analysis
for prefix in jonahi lehilahytsara simmonsi
do
    # Set directories and files
    vcf_dir=/scratch/projects/nib00015/northeastProject/microcebus/$prefix/stacks/vcfFiltering
    vcf_file=$vcf_dir/populations.snps.filtered06.mac3.vcf
    clust_dir=/scratch/projects/nib00015/northeastProject/microcebus/$prefix/admixture
    mkdir -p $clust_dir/logFiles
    # Rename chromosomes
    grep "^#" $vcf_file > $vcf_dir/header.txt 
    grep -v "^#" $vcf_file > $vcf_dir/body.txt
    for i in $(seq 60 91)
    do
        echo "## Processing chromosome $i ..."
        j=$(( $i - 59 ))
        sed -i "s/NC_0336${i}.1/${j}/g" $vcf_dir/body.txt
    done
    cat $vcf_dir/header.txt $vcf_dir/body.txt > $clust_dir/$(basename $vcf_file .vcf).admix.vcf
    rm $vcf_dir/header.txt $vcf_dir/body.txt
    # Convert SNP to bed file
    plink --vcf $admix_dir/$(basename $vcf_file .vcf).admix.vcf --chr-set 32 --make-bed --double-id --out $clust_dir/$(basename $vcf_file .vcf).admix
    # Run admixture
    for K in $(seq 2 $populations)
    do
        sbatch --account=nib00015 --output=$clust_dir/logFiles/clustering.K$K.oe $scripts_dir/clustering.sh $clust_dir/$(basename $vcf_file .vcf).admix.bed $K
    done
done

## Produce cross validation file to determine best K
for prefix in jonahi lehilahytsara simmonsi
do
    # Generate file
    clust_dir=/scratch/projects/nib00015/northeastProject/microcebus/$prefix/admixture
    rm $clust_dir/crossval.txt; touch $clust_dir/crossval.txt
    for K in $(seq 2 $populations)
    do
        cross=$(grep "CV error" $clust_dir/logFiles/admixture.K$K.oe | awk '{print $4'})
        echo -e "$cross\t$K" >> $clust_dir/crossval.txt
    done
    # Output best K
    best=$(sort $clust_dir/crossval.txt | head -1 | awk '{print $2}')
    echo "## Best K for $prefix: $best"
done

###########
## Avahi ##
###########
## Set up an run clustering analysis
# Set directories and files
vcf_dir=/scratch/projects/nib00015/northeastProject/microcebus/lanigerMooreorum/stacks/final/vcfFiltering
vcf_file=$vcf_dir/populations.snps.filtered06.mac3.vcf
clust_dir=/scratch/projects/nib00015/northeastProject/microcebus/lanigerMooreorum/admixture
mkdir -p $clust_dir/logFiles
# Rename chromosomes
grep "^#" $vcf_file > $vcf_dir/header.txt 
grep -v "^#" $vcf_file > $vcf_dir/body.txt
sed -i 's/[^\t]*/1/1' $vcf_dir/body.txt
cat $vcf_dir/header.txt $vcf_dir/body.txt > $clust_dir/$(basename $vcf_file .vcf).admix.vcf
rm $vcf_dir/header.txt $vcf_dir/body.txt
# Convert SNP to bed file
plink --vcf $admix_dir/$(basename $vcf_file .vcf).admix.vcf --make-bed --double-id --out $clust_dir/$(basename $vcf_file .vcf).admix
# Run admixture
for K in $(seq 2 $populations)
do
    sbatch --account=nib00015 --output=$clust_dir/logFiles/clustering.K$K.oe $scripts_dir/clustering.sh $clust_dir/$(basename $vcf_file .vcf).admix.bed $K
done

## Produce cross validation file to determine best K
# Generate file
clust_dir=/scratch/projects/nib00015/northeastProject/microcebus/lanigerMooreorum/admixture
rm $clust_dir/crossval.txt; touch $clust_dir/crossval.txt
for K in $(seq 2 $populations)
do
    cross=$(grep "CV error" $clust_dir/logFiles/admixture.K$K.oe | awk '{print $4'})
    echo -e "$cross\t$K" >> $clust_dir/crossval.txt
done
# Output best K
best=$(sort $clust_dir/crossval.txt | head -1 | awk '{print $2}')
echo "## Best K for $prefix: $best"


######################################
#### PRINCIPAL COMPONENT ANALYSIS ####
######################################
## Principal component analyses were conducted with the R script pca.R


################################
#### INFERENCE OF GENE FLOW ####
################################
## Inference of gene flow rates was conducted with the R script gene_flow_inference.R


###############################
#### EXCESS ALLELE SHARING ####
###############################

################################
## M. jonahi - M. macarthurii ##
################################
vcf_file=/scratch-emmy/projects/nib00015/northeastProject/microcebus/all/stacks/vcfFiltering/populations.snps.filtered06.mac3.vcf
out_dir=/scratch-emmy/projects/nib00015/northeastProject/microcebus/all/dstat
mkdir -p $out_dir
sbatch --account=nib00015 --output=$out_dir/dsuite.oe $SCRIPTS_DIR/dsuite.sh $out_dir/tree.nwk $out_dir/out $vcf_file $out_dir/SETS.txt


###############################
## A. laniger - A. mooreorum ##
###############################
vcf_file=/scratch/usr/nibtve93/northeastProject/avahi/all/stacks/final/vcfFiltering/populations.snps.filtered06.mac3.final.vcf
out_dir=/scratch/usr/nibtve93/northeastProject/avahi/all/dstat
mkdir -p $out_dir
sbatch --output=$out_dir/dsuite.oe $scripts_dir/dsuite.sh $out_dir/tree.nwk $out_dir/out $vcf_file $out_dir/SETS.txt

