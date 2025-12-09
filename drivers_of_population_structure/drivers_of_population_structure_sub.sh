scripts_dir=/home/nibtve93/scripts/popGen

###############################
#### ISOLATION-BY-DISTANCE ####
###############################
## Tests for isolation-by-distance were conducted with the R script ibd.R


################################################
#### ESTIMATED EFFECTIVE MIGRATION SURFACES ####
################################################
## PLINK needs to be included in $PATH (v1.9; https://vcftools.github.io/index.html)
## bed2diffs needs to be included in $PATH (v1; https://github.com/dipetkov/eems/tree/master/bed2diffs)

chrom_file=/scratch/projects/nib00015/mmur3/renameChromosomes.txt
shape_file=$eems_out/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000/River_Mada_1 # Prefix of river shape file
seed=123

################
## Microcebus ##
################
## Before running the code, make sure to:
# Create coord file for individuals manually (in order as given by bcftools query) and save as $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.coord
# Create outer file with extent of migration surfaces manually and save as $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.outer
for prefix in jonahi lehilahytsara simmonsi
do
    vcf_dir=/scratch/projects/nib00015/northeastProject/microcebus/$prefix/stacks/vcfFiltering
    vcf_file=$vcf_dir/populations.snps.filtered06.mac3.vcf
    eems_dir=/scratch/projects/nib00015/northeastProject/microcebus/$prefix/eems/
    eems_out=$eems_dir/results; mkdir -p $eems_out

    # Rename chromosomes in SNP file and convert to bed
    bcftools annotate --rename-chrs $chrom_file $vcf_file > $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.vcf
    plink --vcf $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.vcf --make-bed --double-id --chr-set 32 --out $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom
    
    # Run bed2diffs
    bed2diffs_v1 --bfile $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom --nthreads 4

    # Create config file
    > $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
    echo "datapath = $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
    echo "mcmcpath = $eems_out/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
    echo "nIndiv = $(bcftools query -l $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.vcf | wc -l)" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
    echo "nSites = $(egrep -v "#" $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.vcf | wc -l)" >> $eems_dir/$(basename $VCF_FILE .vcf).eems.renameChrom.ndemes1000.ini
    echo "nDemes = 1000" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
    echo "diploid = true" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
    echo "numMCMCIter = 4000000" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
    echo "numBurnIter = 1000000" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
    echo "numThinIter = 9999" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini

    # Submit job for estimation of migration surfaces
    jid=$(sbatch --account=nib00015 --output=$eems_dir/eems.oe $scripts_dir/eems.sh $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini $seed)

    # Sumit job for plots
    pop_coords=$eems_out/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000/pop_coords.txt # File with coordinates to plot populations on EEMS
    sbatch --account=nib00015 --dependency=afterok:${jid##* } --output=$eems_dir/plot_eems.oe plot_eems.sh $scripts_dir $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000 $eems_out/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000 $pop_coords $shape_file
done

###########
## Avahi ##
###########
## Before running the code, make sure to:
# Create coord file for individuals manually (in order as given by bcftools query) and save as $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.coord
# Create outer file with extent of migration surfaces manually and save as $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.outer
vcf_dir=/scratch/projects/nib00015/northeastProject/avahi/$prefix/stacks/final/vcfFiltering
vcf_file=$vcf_dir/populations.snps.filtered06.mac3.vcf
eems_dir=/scratch/projects/nib00015/northeastProject/avahi/$prefix/eems/
eems_out=$eems_dir/results; mkdir -p $eems_out

# Rename chromosomes in SNP file and convert to bed
grep "^#" $vcf_file > $vcf_dir/header.txt 
grep -v "^#" $vcf_file > $vcf_dir/body.txt
sed -i -E 's/^([0-9]+)/1/g' $vcf_dir/body.txt
cat $vcf_dir/header.txt $vcf_dir/body.txt > $eems_dir/$(basename $vcf_file .vcf).eems.vcf
plink --vcf $EEMS_DIR/$(basename $VCF_FILE .vcf).eems.vcf --make-bed --double-id --out $EEMS_DIR/$(basename $VCF_FILE .vcf).eems

# Run bed2diffs
bed2diffs_v1 --bfile $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom --nthreads 4

# Create config file
> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
echo "datapath = $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
echo "mcmcpath = $eems_out/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
echo "nIndiv = $(bcftools query -l $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.vcf | wc -l)" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
echo "nSites = $(egrep -v "#" $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.vcf | wc -l)" >> $eems_dir/$(basename $VCF_FILE .vcf).eems.renameChrom.ndemes1000.ini
echo "nDemes = 1000" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
echo "diploid = true" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
echo "numMCMCIter = 4000000" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
echo "numBurnIter = 1000000" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini
echo "numThinIter = 9999" >> $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini

# Submit job for estimation of migration surfaces
jid=$(sbatch --parsable --account=nib00015 --output=$eems_dir/eems.oe $scripts_dir/eems.sh $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000.ini $seed)

# Sumit job for plots
pop_coords=$eems_out/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000/pop_coords.txt # File with coordinates to plot populations on EEMS
sbatch --account=nib00015 --dependency=afterok:$jid --output=$eems_dir/plot_eems.oe plot_eems.sh $scripts_dir $eems_dir/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000 $eems_out/$(basename $vcf_file .vcf).eems.renameChrom.ndemes1000 $pop_coords $shape_file
done


#################################
#### ISOLATION-BY-RESISTANCE ####
#################################
################
## Microcebus ##
################
for prefix in jonahi lehilahytsara simmonsi
do
    wd=/scratch-emmy/projects/nib00015/northeastProject/microcebus/$prefix/ibr_radish; mkdir -p $wd/logFiles
    for model in $(cat $wd/models.txt)
    do
        sbatch --account=nib00015 --output=$wd/logFiles/$prefix.$model.oe $scripts_dir/ibr.sh $scripts_dir $wd $prefix $model
    done	
done

###########
## Avahi ##
###########
wd=/scratch-emmy/projects/nib00015/northeastProject/microcebus/laniger/ibr_radish; mkdir -p $wd/logFiles
for model in $(cat $wd/models.txt)
do
    sbatch --account=nib00015 --output=$wd/logFiles/$prefix.$model.oe $scripts_dir/ibr.sh $scripts_dir $wd $prefix $model
done	
