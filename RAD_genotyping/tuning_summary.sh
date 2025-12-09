#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
scripts_dir=$1
stacks_dir=$2
snps_out=$3
fixed_out=$4
plot_dir=$5; mkdir -p $plot_dir

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### tuning_summary.sh: Starting script."
echo -e "#### tuning_summary.sh: Scripts directory: $scripts_dir"
echo -e "#### tuning_summary.sh: Stacks directory: $stacks_dir"
echo -e "#### tuning_summary.sh: Output for estimation of number of number of SNPs and loci: $snps_out"
echo -e "#### tuning_summary.sh: Output for estimation of number of variant and fixed sites: $fixed_out"
echo -e "#### tuning_summary.sh: Plotting directory: $plot_dir \n\n"

################################################################################
#### EXTRACT AND PLOT PARAMETER TUNING SUMMARY STATISTICS ####
################################################################################
echo -e "#### tuning_summary.sh: Extracting parameter tuning summary statistics ...\n"
echo -e 'par_set\tm\tM\tn\tr\tn_snps\tn_loci' > $snps_out
for stacks in $stacks_dir/stacks.m*
do
    echo "processing $stacks"
    for catalog in $stacks/catalog.m*
    do
        echo "...processing $catalog"
        for populations in $catalog/populations.r0.*
        do 
            echo "......processing $populations"
            sed -n '/BEGIN snps_per_loc_postfilters/,/END snps_per_loc_postfilters/ p' $populations/populations.log.distribs | grep -E '^[0-9]' > $populations/snps_per_loc.log 
            paramet_comb=$(echo $catalog | sed "s,${stacks}/catalog.,,g")
            vabe=$(echo $populations | sed "s,${catalog}/populations.r0.,,g")
            paramet_comb=$(echo $paramet_comb.r$vabe)
            paramet_comb_new=$(echo $paramet_comb | sed "s,\.,-,g")
            test=$(echo $paramet_comb_new | grep -o -E '[0-9]+')
            test=(${test[@]}) 
            line_prefix=$(echo $paramet_comb_new\\t${test[0]}\\t${test[1]}\\t${test[2]}\\t${test[3]}\\t)
            sed -r "s/^/$line_prefix/" $populations/snps_per_loc.log >> $snps_out
        done
    done
done

echo -e 'par_set\tm\tM\tn\tr\tTotalSites\tVariantSites\tFixedSites' > $fixed_out
for stacks in $stacks_dir/stacks.m*
do
    echo "processing $stacks"
    for catalog in $stacks/catalog.m*
    do
        echo "...processing $catalog"
        populations=$catalog/populations.r0.8
        echo "......processing $populations"
        paramet_comb=$(echo $catalog | sed "s,${stacks}/catalog.,,g")
        vabe=$(echo $populations | sed "s,${catalog}/populations.r0.,,g")
        paramet_comb=$(echo $paramet_comb.r$vabe)
        paramet_comb_new=$(echo $paramet_comb | sed "s,\.,-,g")
        test=$(echo $paramet_comb_new | grep -o -E '[0-9]+')
        test=(${test[@]})
        values=$(sed -e '1,5d' $populations/populations.sumstats_summary.tsv | awk '{print $3,$4}')
        values=(${values[@]})
        fixed="$((${values[0]}-${values[1]}))"
        echo -e $paramet_comb_new\\t${test[0]}\\t${test[1]}\\t${test[2]}\\t${test[3]}\\t${values[0]}\\t${values[1]}\\t${fixed} >> $fixed_out
    done
done

echo -e "#### tuning_summary.sh: Plotting polymorphic loci ...\n"
Rscript $scripts_dir/plotPolyLoci.R $stacks_dir $plot_dir

echo -e "#### tuning_summary.sh: Plotting SNPs per locus ...\n"
Rscript $scripts_dir/plotSNPsPerLocus.R $stacks_dir $plot_dir

echo -e "#### tuning_summary.sh: Plotting fixed variant sites ...\n"
Rscript $scripts_dir/plotCountFixVarSites.R $stacks_dir $plot_dir

## Report:
echo -e "\n#### tuning_summary.sh: Done with script."
date



