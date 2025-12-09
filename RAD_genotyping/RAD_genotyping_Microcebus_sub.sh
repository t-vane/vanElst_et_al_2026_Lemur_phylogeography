#########################################
#### REFERENCE MAPPING AND FILTERING ####
#########################################
scripts_dir=/home/nibtve93/scripts/referenceMapping

set_id=northeastProject_microcebus
reference_dir=$PWORK/references/mmur3
reference=$reference_dir/GCF_000165445.2_Mmur_3.0_genomic.fna # Reference genome in fasta format

## Index reference genome if not done yet
mkdir -p $reference_dir/logFiles

bwa=true # Boolean specifying whether to create index with BWA
bwa_index=mmur3 # Output name of bwa index
samtools=true # Boolean specifying whether to create index with SAMtools
gatk=true # Boolean specifying whether to create index with GATK (Picard)

sbatch --account=nib00015 --output=$reference_dir/logFiles/indexing.$set_id.oe $scripts_dir/indexing.sh $reference $bwa $bwa_index $samtools $gatk

## Align trimmed reads to reference genome and filter
## The script assumes that trimmed read files are named *.trimmed.1.fq.gz and *.trimmed.2.fq.gz
nt=6
in_dir=$PWORK/trimmedReads/$set_id 
out_dir=$PWORK/bamFiles/$set_id
mkdir -p $out_dir/logFiles

for i in pe se
do
	in_file=$out_dir/map_$i.txt # List of samples (without file extensions) for which reference alignment shall be conducted
	no_inds=$(cat $in_file | wc -l)
	minmapq=20
	
	# Reference mapping, filtering and extraction of genomic regions
	sbatch --account=nib00015 --job-name=map_filter_pip_$i --array=1-$no_inds --output=$out_dir/logFiles/reference_mapping_$i.%A_%a.$set_id.oe $scripts_dir/reference_mapping.sh $i $nt $reference_dir/$bwa_index $in_dir $out_dir $in_file
	
	# Sort and quality filter
	sbatch --account=nib00015 --job-name=map_filter_pip_$i --dependency=singleton --array=1-$no_inds --output=$out_dir/logFiles/quality_filter_$i.%A_%a.$set_id.oe $scripts_dir/quality_filter.sh $i $nt $out_dir $in_file $minmapq
	
	# Deduplicate 
	[[ $i == pe ]] && sbatch --account=nib00015 --job-name=map_filter_pip_$i --dependency=singleton --array=1-$no_inds --output=$out_dir/logFiles/deduplicate_$i.%A_%a.$set_id.oe $scripts_dir/deduplicate.sh $out_dir $in_file $minmapq
	
	# Extract genomic regions
	bed=$out_dir/regionFile_autosomes.bed # BED file with genomic regions that shall be extracted
	exclude="NW_|NC_028718.1|NC_033692.1" # String of chromosomes that are no longer represented (separator: "|")
	suffix=auto # Suffix for naming of final BAM files
	sbatch --account=nib00015  --job-name=map_filter_pip_$i --dependency=singleton --array=1-$no_inds --output=$out_dir/logFiles/extract_regions_$i.%A_%a.$set_id.oe $scripts_dir/extract_regions.sh $i $out_dir $in_file $minmapq $bed "$exclude" $suffix
done


##########################
#### GENOTYPE CALLING ####
##########################
scripts_dir=/home/nibtve93/scripts/stacks
bam_dir=$PWORK/bamFiles/$set_id
suffix=auto # Suffix for final BAM files (see scripts for reference mapping)
nc=16

## Run gstacks
for prefix in all jonahiMacarthurii lehilahytsara simmonsi
do
    stacks_dir=/scratch/projects/nib00015/northeastProject/microcebus/$prefix/stacks
    popmap=$stacks_dir/popmap.txt
    mkdir -p $stacks_dir/logFiles

    main="-I $bam_dir -O $stacks_dir -S .$suffix.bam"
    sbatch --account=nib00015 -c $nc --output=$stacks_dir/logFiles/gstacks.oe $scripts_dir/gstacks.sh "$main" $popmap $nc
done

## Run populations
for prefix in all jonahiMacarthurii lehilahytsara simmonsi
do
    stacks_dir=/scratch/projects/nib00015/northeastProject/microcebus/$prefix/stacks
    pop_dir=$stacks_dir/populations; mkdir -p $pop_dir

    filters="-p 1 -r 0.7"
    output="--vcf"

    sbatch --account=nib00015 -c $nc --output=$stacks_dir/logFiles/populations.oe $scripts_dir/populations.sh $nc $stacks_dir $pop_dir $popmap "$filters" "$output"
done


#######################
#### VCF FILTERING ####
#######################
## Pipeline adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)
scripts_dir=/home/nibtve93/scripts/vcfFiltering
bam_dir=$PWORK/bamFiles/$set_id
reference=$PWORK/references/mmur3/GCF_000165445.2_Mmur_3.0_genomic_reduced.fna
suffix=auto

## Main filtering pipeline
for prefix in all jonahiMacarthurii lehilahytsara simmonsi
do
    vcf_raw=/scratch/projects/nib00015/northeastProject/microcebus/$prefix/stacks/populations/populations.snps.vcf
    vcf_dir=/scratch/projects/nib00015/northeastProject/microcebus/$prefix/stacks/vcfFiltering
    mkdir -p $vcf_dir/logFiles
    ln -s $vcf_raw $vcf_dir

    # Filter for minimum depth
    min_dp=5
    mean_dp=15
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --output=$vcf_dir/logFiles/01_filter_min-dp.oe $scripts_dir/01_filter_min-dp.sh $vcf_raw $min_dp $mean_dp $vcf_dir/$set_id.populations.snps.filtered01.vcf
    # Apply three rounds of filtering for missing data across individuals and genotypes
    maxmiss_geno1=0.5 # One minus maxmimum missingness across genotypes (filtering round 1), i.e., maximum missingness = 1-$maxmiss_geno1
    maxmiss_geno2=0.6 # One minus maxmimum missingness across genotypes (filtering round 2), i.e., maximum missingness = 1-$maxmiss_geno2
    maxmiss_geno3=0.7 # One minus maxmimum missingness across genotypes (filtering round 3), i.e., maximum missingness = 1-$maxmiss_geno3
    filter_inds=FALSE # Boolean specifying whether to filter for missingness across individuals
    maxmiss_ind1=0.9 # Maxmimum missingness across individuals (filtering round 1)
    maxmiss_ind2=0.7 # Maxmimum missingness across individuals (filtering round 2)
    maxmiss_ind3=0.5 # Maxmimum missingness across individuals (filtering round 3)
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --dependency=singleton --output=$vcf_dir/logFiles/02_filter_missing-1.oe $scripts_dir/02_filter_missing-1.sh $vcf_dir/populations.snps.filtered01.vcf $vcf_dir/populations.snps.filtered02.vcf $maxmiss_geno1 $maxmiss_geno2 $maxmiss_geno3 $filter_inds $maxmiss_ind1 $maxmiss_ind2 $maxmiss_ind3
    # Annotate with INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --dependency=singleton --output=$vcf_dir/logFiles/03_annot_gatk.oe $scripts_dir/03_annot_gatk.sh $vcf_dir/populations.snps.filtered02.vcf $vcf_dir/populations.snps.filtered03.vcf $bam_dir $reference $suffix
    # Filter for INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --dependency=singleton --output=$vcf_dir/logFiles/04_filter_gatk.oe $scripts_dir/04_filter_gatk.sh $vcf_dir/populations.snps.filtered03.vcf $vcf_dir/populations.snps.filtered04-soft.vcf $vcf_dir/populations.snps.filtered04-hard.vcf $reference
    # Filter for maximum depth as (mean depth + 2 * standard deviation) / number of individuals
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --dependency=singleton --output=$vcf_dir/logFiles/05_filter_max-dp.oe $scripts_dir/05_filter_max-dp.sh $vcf_dir/populations.snps.filtered04-hard.vcf $vcf_dir/populations.snps.filtered05.vcf 
    # Apply final round of filtering for missing data across individuals and genotypes
    maxmiss_geno=0.9 # One minus maxmimum missingness across genotypes, i.e., maximum missingness = 1-$maxmiss_geno
    filter_inds=FALSE # Boolean specifying whether to filter for missingness across individuals
    maxmiss_ind=0.25 # Maxmimum missingness across individuals
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --dependency=singleton --output=$vcf_dir/logFiles/06_filter_missing-2.oe $scripts_dir/06_filter_missing-2.sh $vcf_dir/populations.snps.filtered05.vcf $vcf_dir/populations.snps.filtered06.vcf $maxmiss_geno $filter_inds $maxmiss_ind
    # Apply minor allele count filter
    mac=3
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --dependency=singleton --output=$vcf_dir/logFiles/07_filter_mac.$set_id.oe $scripts_dir/07_filter_mac.sh $vcf_dir/populations.snps.filtered06.vcf $vcf_dir/populations.snps.filtered06.mac$mac.vcf $mac
    # Create VCF file for M. jonahi (excluding M. macarthurii)
    if [[ "$base" == "jonahiMacarthurii" ]]; then
    rem_string="--remove-indv Mmac_01y06_hely_S12 --remove-indv Mmac_01y07_hely_S12 --remove-indv Mmac_04y06_hely_S12 --remove-indv Mmac_04y13_hely_S12 --remove-indv Mmac_05y08_hely_S12 --remove-indv Mmac_06y08_hely_S12 --remove-indv Mmac_07y08_hely_S12 --remove-indv Mmac_08y08_hely_S12 --remove-indv Mmac_f32y21_vovo --remove-indv Mmac_f35y21_vovo --remove-indv Mmac_m01y06_hely --remove-indv Mmac_m31y21_vovo_S12 --remove-indv Mmac_m34y21_vovo --remove-indv Mmac_m38y21_vovo --remove-indv Mmac_m40y21_vovo"
    vcftools $rem_string --vcf $vcf_dir/populations.snps.filtered06.mac$mac.vcf --recode --recode-INFO-all --stdout > /scratch/projects/nib00015/northeastProject/microcebus/jonahi/stacks/vcfFiltering/populations.snps.filtered06.mac$MAC.vcf
    fi
done
