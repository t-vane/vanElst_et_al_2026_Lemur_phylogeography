# Species-specific responses to paleoclimatic changes and landscape barriers drive contrasting phylogeography of co-distributed lemur species in northeastern Madagascar

This repository holds scripts for the following analyses conducted for the publication [van Elst et al. (2026), *Mol. Ecol.*](https://doi.org/10.1111/mec.70195):
- RAD genotyping
- Phylogenetic inference
- Population structure and gene flow
- Drivers of population structure
- Genetic diversity

Each directory contains one or more scripts from which the respective pipelines are executed, labeled by the suffix "_sub". Input and output files can be found at [Dryad](https://doi.org/10.5061/dryad.4xgxd25nt). 

### RAD genotyping
`./RAD_genotyping` contains scripts for:
- Alignment of trimmed *Microcebus* reads against the *Microcebus murinus* reference genome (Mmur 3.0; [Larsen et al. (2017), *BMC Biol.*](https://doi.org/10.1186/s12915-017-0439-6) with [BWA v0.7.17](https://github.com/lh3/bwa)) and subsequent filtering of BAM files with [SAMtools v1.11](http://www.htslib.org/)
- Genotype calling with [Stacks v2.53](https://catchenlab.life.illinois.edu/stacks/) using a reference-guided approach for *Microcebus* and a *de novo* approach for *Avahi* samples.
- VCF filtering with [VCFtools v0.1.17](https://vcftools.github.io/index.html) and [GATK v3.8.1/v4.1.9.0](https://gatk.broadinstitute.org/hc/en-us)

### Phylogenetic inference
`./phylogenetic_inference` contains scripts to infer maximum likelihood phylogenies with [IQ-TREE v2.2.0](https://iqtree.github.io/).

### Population structure and gene flow
`./population_structure_and_gene_flow` contains scripts for:
- Clustering analysis with [ADMIXTURE v1.3.0](https://dalexander.github.io/admixture/)
- Principal component analysis with the R package ['SNPRelate' v1.32.2](https://github.com/zhengxwen/SNPRelate)
- Estimation of gene flow rates under a coalescent model with the R package ['gene.flow.inference' v0.0.0.9](https://github.com/elundgre/gene-flow-inference)
- Testing for excess allele sharing with [Dsuite v0.4 r38](https://github.com/millanek/Dsuite)

### Drivers of population structure
`./drivers_of_population_structure` contains scripts for:
- Testing for isolation-by-distance with the R package ['vegan' v2.5-7](https://vegandevs.github.io/vegan/)
- Visualization of spatial deviations from a model of isolation-by-distance with [EEMS](https://github.com/dipetkov/eems)
- Isolation-by-resistance modeling with the R package ['radish' v0.0.2](https://github.com/nspope/radish) and associated inference of spatial rasters with the R packages ['terra' v1.6.47](https://rspatial.github.io/terra/), ['rasterdiv' v0.3.4](https://mattmar.github.io/rasterdiv/), ['whitebox' v2.2.0](https://whiteboxr.gishub.org/), ['ENMtools' v1.1.2](https://github.com/danlwarren/ENMTools) and ['ENMeval' v2.0.4](https://jamiemkass.github.io/ENMeval)

### Genetic diversity
`./genetic_diversity` contains scripts to infer individual observed heterozygosities with [VCFtoosl v0.1.17](https://vcftools.github.io/index.html).
