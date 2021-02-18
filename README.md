# Pipeline for analysis of Z-linked and autosomal genomic variation in barn swallows

Below are details about the data processing and analysis steps in our analysis of genetic diversity and differentiation across the genomes of barn swallow populations, focusing on comparisons between the Z chromosome and autosomes to understand the mechanisms governing sex-linked variation and the role of the Z chromosome in speciation. 

The steps described here rely on the following software:

* conda
* trimmomatic
* bwa
* GATK (v3.8-1-0 and v4.0.8.1)
* samtools/bcftools/htslib/bgzip/tabix
* vcftools
* bedtools
* [Pixy](https://pixy.readthedocs.io/en/latest/)
* ADMIXTURE
* [Genomics general scripts](https://github.com/simonhmartin/genomics_general)
* TWISST
* R

Lists and miscellaneous files are in the `processing_files` directory.
Shell and Python scripts are in the `scripts` directory.
Population genetic summary statistics output from `pixy` are in the `pixy_results` directory.
R scripts used in analyses are in the `R` directory.

Note that you may need to adjust the organization of file locations to suite your environment.

## Contents

* [Read processing](#read-processing)
* [Read filtering](#read-filtering)
* [Read mapping](#read-mapping)
* [Variant calling](#variant-calling)
* [Variant filtering](#variant-filtering)
* [Pixy analysis](#pixy-analysis)
* [Tajima's D analysis](#tajimas-d-analysis)
* [Fst and PBS analysis]($fst-pbs-analysis)
* [ADMIXTURE analysis](#admixture-analysis)
* [Topology weighting analysis](#topology-weighting-analysis)
* [ABBA-BABA analysis](#abba-baba-analysis)
* [Analysis in R](#analysis-in-r)
