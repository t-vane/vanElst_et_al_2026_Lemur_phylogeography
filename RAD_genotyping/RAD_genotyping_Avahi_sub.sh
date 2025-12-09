scripts_dir=/home/nibtve93/scripts/stacks
bam_dir=$PWORK/trimmedReads/northeastProject_avahi
nc=32

########################################################################################################
#### STACKS PARAMETER TUNING FOLLOWING PARIS ET AL. (2017); https://doi.org/10.1111/2041-210X.12775 ####
########################################################################################################
## Script adapted and modified from Gabriele M. Sgarlata

for prefix in all lanigerMooreorum
do
    stacks_dir=/scratch/projects/nib00015/northeastProject/avahi/$prefix/stacks/parameterTuning
    in_file=$stacks_dir/../inds.txt
    no_inds=$(cat $in_file | wc -l)
    popmap=$stacks_dir/../popmap.txt
    mkdir $stacks_dir/logFiles

    # Step 0: Create reads with uniform length distribution
    length=142
    jid_cutadapt=""
    for ind in $(cat $in_file)
    do
        jid=$(sbatch --account=nib00015 --parsable --output=$stacks_dir/logFiles/cutadapt.$ind.oe $scripts_dir/cutadapt.sh $length $ind $bam_dir)
        jid_cutadapt="$jid_cutadapt:$jid"
    done 
    jid_cutadapt=${jid_cutadapt#:} 

    # Step 1: ustacks (build loci for each sample)
    jid_ustacks=""
    for M in 2 3 4 5 6 7 8
    do
        m=3
        out_dir=$stacks_dir/stacks.m$m.M$M
        mkdir -p $out_dir
        jid=$(sbatch --array=1-${no_inds} --account=nib00015 --parsable --dependency=afterok:$jid_cutadapt --output=$stacks_dir/logFiles/ustacks.m$m.M$M.%A_%a.oe $scripts_dir/ustacks.sh $bam_dir $in_file $out_dir $m $M $nc)
        jid_ustacks="$jid_ustacks:$jid"
    done

    for m in 2 4 5 6 7 8
    do
        M=2
        out_dir=$stacks_dir/stacks.m$m.M$M
        mkdir -p $out_dir
        jid=$(sbatch --array=1-${no_inds} --account=nib00015 --parsable --dependency=afterok:$jid_cutadapt --output=$stacks_dir/logFiles/ustacks.m$m.M$M.%A_%a.oe $scripts_dir/ustacks.sh $bam_dir $in_file $out_dir $m $M $nc)
        jid_ustacks="$jid_ustacks:$jid"
    done
    jid_ustacks=${jid_ustacks#:}

    # Step 2: cstacks (build catalog of loci)
    jid_cstacks=""
    for M in 2 3 4 5 6 7 8
    do
        m=3
        n=1
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n; mkdir -p $catalog_dir
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_ustacks --output=$stacks_dir/logFiles/cstacks.m$m.M$M.n$n.oe $scripts_dir/cstacks.sh $stacks_dir $catalog_dir $popmap $popmap $m $M $n $nc)
        jid_cstacks="$jid_cstacks:$jid"
    done

    for m in 2 4 5 6 7 8
    do
        M=2
        n=1
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n; mkdir -p $catalog_dir
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_ustacks --output=$stacks_dir/logFiles/cstacks.m$m.M$M.n$n.oe $scripts_dir/cstacks.sh $stacks_dir $catalog_dir $popmap $popmap $m $M $n $nc)
        jid_cstacks="$jid_cstacks:$jid"
    done

    for n in 2 3 4 5 6 7 8
    do
        M=2
        m=3
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n; mkdir -p $catalog_dir
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_ustacks --output=$stacks_dir/logFiles/cstacks.m$m.M$M.n$n.oe $scripts_dir/cstacks.sh $stacks_dir $catalog_dir $popmap $popmap $m $M $n $nc)
        jid_cstacks="$jid_cstacks:$jid"
    done
    jid_cstacks=${jid_cstacks#:}

    # Step 3: sstacks (match samples against catalog)
    jid_sstacks=""
    for M in 2 3 4 5 6 7 8
    do
        m=3
        n=1
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_cstacks --output=$stacks_dir/logFiles/sstacks.m$m.M$M.n$n.oe $scripts_dir/sstacks.sh $catalog_dir $popmap $nc)
        jid_sstacks="$jid_sstacks:$jid"
    done

    for m in 2 4 5 6 7 8
    do
        M=2
        n=1
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_cstacks --output=$stacks_dir/logFiles/sstacks.m$m.M$M.n$n.oe $scripts_dir/sstacks.sh $catalog_dir $popmap $nc)
        jid_sstacks="$jid_sstacks:$jid"
    done

    for n in 2 3 4 5 6 7 8
    do
        M=2
        m=3
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_cstacks --output=$stacks_dir/logFiles/sstacks.m$m.M$M.n$n.oe $scripts_dir/sstacks.sh $catalog_dir $popmap $nc)
        jid_sstacks="$jid_sstacks:$jid"
    done
    jid_sstacks=${jid_sstacks#:}

    # Step 4: tsv2bam (transpose data to store by locus and incorporate PE reads)
    jid_tsv2bam=""
    for M in 2 3 4 5 6 7 8
    do
        m=3
        n=1
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_sstacks --output=$stacks_dir/logFiles/tsv2bam.m$m.M$M.n$n.oe $scripts_dir/tsv2bam.sh $catalog_dir $bam_dir $popmap $nc)
        jid_tsv2bam="$jid_tsv2bam:$jid"
    done

    for m in 2 4 5 6 7 8
    do
        M=2
        n=1
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_sstacks --output=$stacks_dir/logFiles/tsv2bam.m$m.M$M.n$n.oe $scripts_dir/tsv2bam.sh $catalog_dir $bam_dir $popmap $nc)
        jid_tsv2bam="$jid_tsv2bam:$jid"
    done

    for n in 2 3 4 5 6 7 8
    do
        M=2
        m=3
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_sstacks --output=$stacks_dir/logFiles/tsv2bam.m$m.M$M.n$n.oe $scripts_dir/tsv2bam.sh $catalog_dir $bam_dir $popmap $nc)
        jid_tsv2bam="$jid_tsv2bam:$jid"
    done
    jid_tsv2bam=${jid_tsv2bam#:}

    # Step 5: gstacks (variant calling and genotyping)
    jid_gstacks=""
    for M in 2 3 4 5 6 7 8
    do
        m=3
        n=1
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        main="-P $catalog_dir --rm-pcr-duplicates --rm-unpaired-reads"
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_tsv2bam --output=$stacks_dir/logFiles/gstacks.m$m.M$M.n$n.oe $scripts_dir/gstacks.sh "$main" $popmap $nc)
        jid_gstacks="$jid_gstacks:$jid"
    done

    for m in 2 4 5 6 7 8
    do
        M=2
        n=1
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        main="-P $catalog_dir --rm-pcr-duplicates --rm-unpaired-reads"
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_tsv2bam --output=$stacks_dir/logFiles/gstacks.m$m.M$M.n$n.oe $scripts_dir/gstacks.sh "$main" $popmap $nc)
        jid_gstacks="$jid_gstacks:$jid"
    done

    for n in 2 3 4 5 6 7 8
    do
        M=2
        m=3
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        main="-P $catalog_dir --rm-pcr-duplicates --rm-unpaired-reads"
        jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_tsv2bam --output=$stacks_dir/logFiles/gstacks.m$m.M$M.n$n.oe $scripts_dir/gstacks.sh "$main" $popmap $nc)
        jid_gstacks="$jid_gstacks:$jid"
    done
    jid_gstacks=${jid_gstacks#:}

    # Step 6: populations (statistics and exporting)
    jid_populations=""
    for M in 2 3 4 5 6 7 8
    do
        n=1
        m=3
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        for r in 0.4 0.6 0.8
        do
            jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_gstacks --output=$stacks_dir/logFiles/populations.m$m.M$M.n$n.r$r.oe $scripts_dir/populations.sh $nc $catalog_dir $catalog_dir/populations.r$r $popmap)
            jid_populations="$jid_populations:$jid"
        done
    done

    for m in 2 4 5 6 7 8
    do
        M=2
        n=1
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        for r in 0.4 0.6 0.8
        do
            jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_gstacks --output=$stacks_dir/logFiles/populations.m$m.M$M.n$n.r$r.oe $scripts_dir/populations.sh $nc $catalog_dir $catalog_dir/populations.r$r $popmap)
            jid_populations="$jid_populations:$jid"
        done
    done

    for n in 2 3 4 5 6 7 8
    do
        M=2
        m=3
        catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n
        for r in 0.4 0.6 0.8
        do
            jid=$(sbatch --account=nib00015 --parsable --dependency=afterok:$jid_gstacks --output=$stacks_dir/logFiles/populations.m$m.M$M.n$n.r$r.oe $scripts_dir/populations.sh $nc $catalog_dir $catalog_dir/populations.r$r $popmap)
            jid_populations="$jid_populations:$jid"
        done
    done
    jid_populations=${jid_populations#:}

    # Step 7: Estimate and plot number of SNPs and loci as well as variant and fixed sites for each parameter combination
    sbatch --account=nib00015 --dependency=afterok:$jid_populations --output=$stacks_dir/logFiles/tuning_summary.oe $scripts_dir/tuning_summary.sh $scripts_dir $stacks_dir $stacks_dir/n_snps_per_locus.tsv $stacks_dir/count_fixed_variant.tsv $stacks_dir/plots
done


################################
#### FINAL GENOTYPE CALLING ####
################################
## Best parameter combination was estimated to be M=2, m=2, n=3

for prefix in all lanigerMooreorum
do
    stacks_dir=/scratch/projects/nib00015/northeastProject/avahi/$prefix/stacks/final
    in_file=$stacks_dir/../inds.txt
    no_inds=$(cat $in_file | wc -l)
    popmap=$stacks_dir/../popmap.txt
    M=2
    m=2
    n=3
    catalog_dir=$stacks_dir/stacks.m$m.M$M/catalog.m$m.M$M.n$n; mkdir -p $catalog_dir

    mkdir $stacks_dir/logFiles

    ## ustacks (build loci for each sample)
    # Can be taken directly from parameter tuning directory
    cp -r $stacks_dir/../parameterTuning/stacks.m$m.M$M $stacks_dir
    rm -r $stacks_dir/stacks.m2.M2/catalog.m$m.M$M.n1/
    ln -s $stacks_dir/stacks.m$m.M$M/*gz $catalog_dir

    ## cstacks (build catalog of loci)
    jid=$(sbatch --parsable --account=nib00015 --output=$stacks_dir/logFiles/cstacks.m$m.M$M.n$n.oe $scripts_dir/cstacks_denovo.sh $catalog_dir $popmap $n $nc)

    ## sstacks (match samples against catalog)
    jid=$(sbatch --parsable --account=nib00015 --dependency=afterok:$jid --output=$stacks_dir/logFiles/sstacks.m$m.M$M.n$n.oe $scripts_dir/sstacks_denovo.sh $catalog_dir $popmap $nc)

    ## tsv2bam (transpose data to store by locus and incorporate PE reads)
    jid=$(sbatch --parsable --account=nib00015 --dependency=afterok:$jid --output=$stacks_dir/logFiles/tsv2bam.m$m.M$M.n$n.oe $scripts_dir/tsv2bam_denovo.sh $catalog_dir $bam_dir $popmap $nc)

    ## gstacks (variant calling and genotyping)
    main="-P $catalog_dir --rm-pcr-duplicates --rm-unpaired-reads"
    jid=$(sbatch --parsable --account=nib00015 --dependency=afterok:$jid --output=$stacks_dir/logFiles/gstacks.m$m.M$M.n$n.oe $scripts_dir/gstacks_denovo.sh "$main" $popmap $nc)

    ## populations (statistics and exporting)
    r=0.7
    filters="-p 1 -r $r"
    output="--vcf --fasta-samples --fasta-samples-raw --plink"
    sbatch --account=nib00015 --dependency=afterok:$jid --output=$stacks_dir/logFiles/populations.m$m.M$M.n$n.r$r.oe $scripts_dir/populations.sh $nc $catalog_dir $catalog_dir/populations.r$r $popmap "$filters" "$output"


#######################
#### VCF FILTERING ####
#######################
## Pipeline adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)
scripts_dir=/home/nibtve93/scripts/vcfFiltering

## Main filtering pipeline
for prefix in all lanigerMooreorum
do
    vcf_raw=/scratch/projects/nib00015/northeastProject/avahi/$prefix/stacks/final/stacks.m2.M2/catalog.m2.M2.n3/populations.r0.7/populations.snps.vcf
    vcf_dir=/scratch/projects/nib00015/northeastProject/avahi/$prefix/stacks/final/vcfFiltering
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
    # Filter for maximum depth as (mean depth + 2 * standard deviation) / number of individuals
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --dependency=singleton --output=$vcf_dir/logFiles/05_filter_max-dp.oe $scripts_dir/05_filter_max-dp.sh $vcf_dir/populations.snps.filtered02.vcf $vcf_dir/populations.snps.filtered05.vcf 
    # Apply final round of filtering for missing data across individuals and genotypes
    maxmiss_geno=0.9 # One minus maxmimum missingness across genotypes, i.e., maximum missingness = 1-$maxmiss_geno
    filter_inds=FALSE # Boolean specifying whether to filter for missingness across individuals
    maxmiss_ind=0.25 # Maxmimum missingness across individuals
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --dependency=singleton --output=$vcf_dir/logFiles/06_filter_missing-2.oe $scripts_dir/06_filter_missing-2.sh $vcf_dir/populations.snps.filtered05.vcf $vcf_dir/populations.snps.filtered06.vcf $maxmiss_geno $filter_inds $maxmiss_ind
    # Apply minor allele count filter
    mac=3
    sbatch --account=nib00015 --job-name=${prefix}_vcffilt --dependency=singleton --output=$vcf_dir/logFiles/07_filter_mac.$set_id.oe $scripts_dir/07_filter_mac.sh $vcf_dir/populations.snps.filtered06.vcf $vcf_dir/populations.snps.filtered06.mac$mac.vcf $mac
    # Create VCF file for M. jonahi (excluding M. macarthurii)
    if [[ "$base" == "lanigerMooreorum" ]]; then
    rem_string="--remove-indv Amoo_AvahiFizono01_fizo_S12 --remove-indv Amoo_AvahiFizono02_fizo --remove-indv Amoo_AvahiFizono04_fizo --remove-indv Amoo_f02y21_habe_S12 --remove-indv Amoo_m01y21_habe_S12"
    vcftools $rem_string --vcf $vcf_dir/populations.snps.filtered06.mac$mac.vcf --recode --recode-INFO-all --stdout > /scratch/projects/nib00015/northeastProject/avahi/laniger/stacks/vcfFiltering/populations.snps.filtered06.mac$MAC.vcf
    fi
done

