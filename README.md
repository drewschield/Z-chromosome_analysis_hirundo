# Analysis of Z-linked genomic variation in barn swallows

Below are details about the data processing and analysis steps in our analysis of genetic diversity and differentiation across the genomes of barn swallow populations, focusing on comparisons between the Z chromosome and autosomes to understand the mechanisms governing sex-linked variation and the role of the Z chromosome in speciation. This workflow is a companion to the description in Schield et al. (in review).

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

Note that you will need to adjust the organization of file locations and paths to suit your environment.

## Contents

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

### Read filtering

Raw whole-genome resequencing read data are available at the NCBI [SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA323498/)

Illumina libraries were sequenced on two lanes, so we will concatenate raw read files and quality filter the concatenated input.

We will impose these filters to trim reads:

* Remove 5' end bases if quality is below 20
* Remove 3' end bases if quality is below 20
* Minimum read length = 32
* Remove reads if average quality is < 30

#### Set up environment

Get raw fastq data into `fastq` directory. <br /> Make a `fastq_filtered` directory for output.

```
mkdir fastq
mkdir fastq_filtered
```

#### Concatenate raw data and filter reads with `trimmomatic`

The script below will concatenate the data and run trimmomatic on the paired reads for samples in `processing_files/sample.orig.list`.

cat_trimmomatic.sh:

```
for line in `cat sample.orig.list`; do
	name=$line
	echo processing and filtering ${name}.
	cat /data1/hirundo_data_master/WGS/BarnSwallows/fastq/${name}_*_1.fq.gz > ./fastq/${name}_1.fq.gz
	cat /data1/hirundo_data_master/WGS/BarnSwallows/fastq/${name}_*_2.fq.gz > ./fastq/${name}_2.fq.gz
	trimmomatic PE -phred33 -threads 16 ./fastq/${name}_1.fq.gz ./fastq/${name}_2.fq.gz ./fastq_filtered/${name}_1_P.trim.fq.gz ./fastq_filtered/${name}_1_U.trim.fq.gz ./fastq_filtered/${name}_2_P.trim.fq.gz ./fastq_filtered/${name}_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30
	rm ./fastq/${name}_*.fq.gz
done
```

`sh cat_trimmomatic.sh`

### Read mapping

We will map filtered read data to the [Hirundo rustica (Chelidonia) reference genome](http://gigadb.org/dataset/view/id/100531) using `bwa`.