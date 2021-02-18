# Analysis of Z-linked genomic variation in barn swallows

Below are details about the data processing and analysis steps in our analysis of genetic diversity and differentiation across the genomes of barn swallow populations, focusing on comparisons between the Z chromosome and autosomes to understand the mechanisms governing sex-linked variation and the role of the Z chromosome in speciation. This workflow is a companion to the description in Schield et al. (in review).

The steps below depend on the following software and assume that dependencies are on the user path:

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
* [Quantifying mapping results](#quantifying-mapping-results)
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

__*Raw whole-genome resequencing read data are available at the NCBI [SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA323498/).*__

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

#### Run `trimmomatic` on outgroup *Hirundo smithii* sample

```
trimmomatic PE -phred33 -threads 16 ./fastq/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_1.fq.gz ./fastq/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_2.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_1_P.trim.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_1_U.trim.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_2_P.trim.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30
```

### Read mapping

We will map filtered read data to the [Hirundo rustica (Chelidonia) reference genome](http://gigadb.org/dataset/view/id/100531) using `bwa`.

#### Set up environment, map reads with `bwa`, sort with `samtools`

`mkdir bam`

The script below will map filtered reads for samples in `processing_files/sample.list` and also rename the output with helpful 'HR' prefixes.

bwa_mem.sh:

```
for line in `cat sample.list`; do
	name=$line
	out=`echo $name | sed 's/L_/HR/'`
	echo "Mapping filtered $out data to reference"
	bwa mem -t 16 -R "@RG\tID:$out\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:$out" Hirundo_rustica_Chelidonia.fasta ./fastq_filtered/${name}_1_P.trim.fq.gz ./fastq_filtered/${name}_2_P.trim.fq.gz | samtools sort -@ 16 -O bam -T temp -o ./bam/$out.bam -
done
```

`sh bwa_mem.sh ./processing_files/sample.list`

#### Run `bwa` on outgroup *Hirundo smithii* sample

```
bwa mem -t 16 -R "@RG\tID:RS_5\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:RS_5" Hirundo_rustica_Chelidonia.fasta ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_1_P.trim.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_2_P.trim.fq.gz | samtools sort -@ 16 -O bam -T temp -o ./bam/RS_5.bam -
```

### Quantify mapping results

#### Generate samtools index for reference genome

`samtools faidx Hirundo_rustica_Chelidonia.fasta`

#### Index mapping files

`for i in ./bam/*.bam; do samtools index -@ 8 $i; done`

#### Index outgroup *Hirundo smithii* sample

`samtools index -@ 8 ./bam/RS_5.bam`

#### Output mapping statistics using `samtools`

```
cd ./bam
mkdir stats
```

The script below will calculate mapping statistics per bam file.

samtools_stat.sh:

```
for i in *.bam; do
	echo calculating mapping statistics for $i
	samtools stat -@ 8 $i > ./stats/$i.stat.txt
done
```

`sh samtools_stat.sh`

#### Calculate rough coverage estimate (assuming a 1.21 Gb genome size)

```
for i in ./stats/*.stat.txt; do echo $i; grep 'bases mapped:' $i | awk '{print $4/121000000}'; done
```

### Variant calling

We will use `GATK` for variant discovery.

These steps use local installations of GATK 3.8-1-0 and 4.0.8.1

#### Set up environment

```
mkdir gvcf
mkdir vcf
```

#### Generate sequence dictionary 

`./gatk-4.0.8.1/gatk CreateSequenceDictionary -R Hirundo_rustica_Chelidonia.fasta`

#### Call individual variants / generate genomic VCF per sample using `HaplotypeCaller`

The script below will call `HaplotypeCaller` on bam files for samples in `processing_files/sample.list`.

GATK_HaplotypeCaller.sh:

```
list=$1
for i in `cat $list`; do
	./gatk-4.0.8.1/gatk HaplotypeCaller -R Hirundo_rustica_Chelidonia.fasta --ERC GVCF -I ./bam/$i.bam -O ./gvcf/$i.raw.snps.indels.g.vcf
	bgzip ./gvcf/$i.raw.snps.indels.g.vcf
done
```

`sh GATK_HaplotypeCaller.sh ./processing_files/sample.list`

*Note: GATK HaplotypeCaller will take a small eternity to run on all of these samples one-by-one. Consider breaking up the job into smaller lists of samples and running jobs in parallel.*

#### Run `HaplotypeCaller` on outgroup *Hirundo smithii* sample

```
./gatk-4.0.8.1/gatk HaplotypeCaller -R Hirundo_rustica_Chelidonia.fasta --ERC GVCF -I ./bam/RS_5.bam -O ./gvcf/RS_5.raw.snps.indels.g.vcf
bgzip ./gvcf/RS_5.raw.snps.indels.g.vcf
tabix -p vcf ./gvcf/RS_5.raw.snps.indels.g.vcf.gz
```

#### Call variants among cohort of samples using `GenotypeGVCFs`

This will call an 'all-sites' VCF among individuals in `processing_files/sample+smithii.gvcf.list`

```
cd vcf
../gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ../Hirundo_rustica_Chelidonia.fasta -V sample+smithii.gvcf.list -allSites -o hirundo_rustica+smithii.allsites.raw.vcf.gz
```

### Variant filtering

Variant filtration will proceed with a few general steps:
1. Impose hard quality filters
2. Remove indels and repeats
3. Remove sites with extreme read depths and on scaffolds not assigned to chromosomes
4. Remove female heterozygous sites on the Z chromosome
5. Set additional filters for specific analyses

#### Annotate variants not passing quality filters using `VariantFiltration`

We will impose these __hard filters__ to remove low-quality variants:

* QD < 2.0
* FS > 60.0
* MQ < 40.0
* MQRankSum < -12.5
* ReadPosRankSum < -8.0

```
../gatk-4.0.8.1/gatk VariantFiltration -V hirundo_rustica+smithii.allsites.raw.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" --mask ../genome_annotation/GCA_003692655.1_Chelidonia_genomic.repeat.sort.bed --mask-name REP -O hirundo_rustica+smithii.allsites.HardFilter.vcf.gz
```

#### Recode indels, repeats, and sites failing quality filters as missing genotypes using `bcftools`

```
bcftools filter --threads 20 -e 'TYPE="indel" || FILTER="REP" || FILTER="QD2" || FILTER="FS60" || FILTER="MQ40" || FILTER="MQRankSum-12.5" || FILTER="ReadPosRankSum-8"' --set-GTs . -O z -o hirundo_rustica+smithii.allsites.HardFilter.recode.vcf.gz hirundo_rustica+smithii.allsites.HardFilter.vcf.gz
tabix -p vcf hirundo_rustica+smithii.allsites.HardFilter.recode.vcf.gz
```

#### Filter sites with extreme read depths and on unassigned scaffolds

Use `vcftools` to calculate mean depth per site:

```
vcftools --gzvcf hirundo_rustica+smithii.allsites.HardFilter.recode.vcf.gz --site-mean-depth --out hirundo_rustica+smithii.allsites.HardFilter.recode
```

Note: there is likely a faster `bcftools` solution to this, using the `query` subcommand to output information from the FORMAT field.

Based on these data, the mean depth = 4, 2.5th depth quantile = 0, and 97.5th depth quantile = 9.54. We'll want to recode sites with super high read depths as missing genotypes to avoid effects of potential paralogous mappings.

We'll recode sites with mean depth above 9.5 as missing data and remove sites on unassigned scaffolds using `bcftools`.

This command will output sites on scaffolds listed in `processing_files/Hirundo_rustica_Barn2Flycatcher_ChromAssigned.bed` and recode sites with extremely high mean read depth:

```
bcftools view --threads 16 -R Hirundo_rustica_Barn2Flycatcher_ChromAssigned.bed hirundo_rustica+smithii.allsites.HardFilter.recode.vcf.gz | bcftools filter --threads 16 -e 'MEAN(FORMAT/DP)>9.5' --set-GTs . -O z -o hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.vcf.gz
tabix -p vcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.vcf.gz
```

#### Identify and remove __female heterozygous sites__ on the Z chromosome

Female barn swallows are hemizygous and cannot have heterozygous genotypes on the Z chromosome. Heterozygous variant calls on the Z chromosome in females are therefore spurious and should be removed prior to analysis. We'll conservatively recode any sites with heterozygous genotypes in females as missing data for all individuals.

Extract biallelic SNPs on the Z chromosome to query variants for female heterozygous sites:

```
bcftools view --threads 16 -m2 -M2 -U -v snps -r QRBI01000206.1,QRBI01000143.1,QRBI01000102.1,QRBI01000186.1,QRBI01000173.1,QRBI01000125.1,QRBI01000103.1,QRBI01000275.1,QRBI01000095.1,QRBI01000244.1,QRBI01000092.1,QRBI01000140.1,QRBI01000225.1,QRBI01000195.1,QRBI01000141.1,QRBI01000097.1,QRBI01000234.1,QRBI01000139.1,QRBI01000052.1 -O v -o hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.snps.chrZ.vcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.vcf.gz
```

The Python script `identify_female_Zhet_sites.py` will identify sites in the Z-linked VCF with any heterozygous calls in females, based on samples in `processing_files/sample.female.list`:

```
python identify_female_Zhet_sites.py ./processing_files/sample.female.list hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.snps.chrZ.vcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.snps.chrZ.female_Zhet_sites.txt
```

Format BED file for masking sites:

```
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.snps.chrZ.female_Zhet_sites.txt > annot.female_Zhet.bed
```

Index BED feature file using GATK `IndexFeatureFile`:

`../gatk-4.0.8.1/gatk IndexFeatureFile --feature-file annot.female_Zhet.bed`

Mask sites with `VariantFiltration`:

```
../gatk-4.0.8.1/gatk VariantFiltration -V hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.vcf.gz --mask annot.female_Zhet.bed --mask-name ZHET -O hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.female_Zhet.vcf.gz
```

Recode sites as missing genotypes using `bcftools`:

```
bcftools filter --threads 16 -e 'FILTER="ZHET"' --set-GTs . -O z -o hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.vcf.gz hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.female_Zhet.vcf.gz
tabix -p vcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.vcf.gz
```

















































