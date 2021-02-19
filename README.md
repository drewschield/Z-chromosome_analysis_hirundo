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
* PhyML
* TWISST
* R

Lists and miscellaneous files are in the `processing_files` directory.
Shell and Python scripts are in the `scripts` directory.
Population genetic summary statistics output from `pixy` are in the `pixy_results` directory.
Tajima's *D* statistics are in the `tajimas_d_results` directory.
Relative population differentiation statistics (*Fst* and *PBS*) are in the `fst_results` and `pbs_results` directories.
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

#### Identify and remove *female heterozygous sites* on the Z chromosome

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

#### Set additional filters for various analyses

Certain analyses we'll perform require further removal of sites not meeting specific criteria. Specifically, we'll produce VCFs where we retain:

* SNPs present in >= 60% of samples
	* Non-singleton SNPs
	* SNPs with minor-allele frequency (MAF) >= 0.05
	* SNPs with MAF >= 0.05 and no closer than 100 bp to the nearest SNP

First, extract biallelic SNPs using `bcftools`:

```
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.vcf.gz hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.vcf.gz
```

Keep SNPs meeting missing data threshold (i.e., --max-missing 0.4) using `vcftools`:

```
vcftools --gzvcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.vcf.gz --recode --stdout --max-missing 0.4 | bgzip -c > hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.vcf.gz
```

Remove singletons:

```
vcftools --gzvcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.vcf.gz --recode --stdout --mac 2 | bgzip -c > hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.mac2.vcf.gz
```

Remove SNPs with MAF < 0.05:

```
vcftools --gzvcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.vcf.gz --recode --stdout --mac 2 --maf 0.05  | bgzip -c > hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.vcf.gz
```

Thin SNPs with MAF >= 0.05 by 100 bp (and only keep ingroup samples based on `processing_files/sample.ingroup.list`):

```
vcftools --gzvcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.vcf.gz --recode --stdout --keep ./processing_files/sample.ingroup.list --min-alleles 2 --max-alleles 2 --maf 0.05 --thin 100  | bgzip -c > hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.vcf.gz
```

Bonus: extract autosomal and Z-linked SNPs using `bcftools` based on `processing_files/Hirundo_rustica_Barn2Flycatcher_AutoAssigned.bed` and `processing_files/Hirundo_rustica_Barn2Flycatcher_ZAssigned.bed`:

```
bcftools view -R Hirundo_rustica_Barn2Flycatcher_ZAssigned.bed -O v -o hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ.vcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.vcf.gz
bcftools view -R Hirundo_rustica_Barn2Flycatcher_AutoAssigned.bed -O v -o hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto.vcf hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.vcf.gz
```

### Pixy analysis

We'll use [pixy](https://pixy.readthedocs.io/en/latest/) to estimate nucleotide diversity (pi) for barn swallow populations.

I installed `pixy` in it's own `conda` environment with Python v3.6.

#### Set up environment

We'll make an `analysis` directory and run the various analysis steps below within associated subdirectories.

```
mkdir analysis
cd analysis
mkdir pixy
cd pixy
mkdir pixy_results
mkdir pixy_zarr
```

#### Subset all-sites VCF by scaffold

Pixy will run on each scaffold independently.

`mkdir vcf_all-sites_chrom`

The script below will parse an input VCf into component scaffold-specific VCFs in `processing_files/Hirundo_rustica_Barn2Flycatcher_ChromAssigned.list`.

parse_chrom_all-sites_VCF.sh:

```
for scaff in `cat Hirundo_rustica_Barn2Flycatcher_ChromAssigned.list`; do
	echo parsing $scaff VCF
	bcftools view --threads 8 -r $scaff -O z -o ./vcf_all-sites_chrom/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.$scaff.vcf.gz hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.vcf.gz
	tabix -p vcf ./vcf_all-sites_chrom/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.$scaff.vcf.gz
	echo indexed $scaff VCF
done
```

`sh parse_chrom_all-sites_VCF.sh`

Ready conda environment for pixy analysis:

```
conda deactivate
conda activate pixy
```

The script below will run `pixy` on each scaffold VCF in `processing_files/Hirundo_rustica_Barn2Flycatcher_ChromAssigned.list` for populations in `processing_files/hirundo.popmap`.

pixyloop.sh:

```
list=$1
for chrom in `cat $list`; do
	pixy --stats pi --vcf ../../vcf/vcf_all-sites_chrom/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.$chrom.vcf.gz --zarr_path pixy_zarr --window_size 100000 --populations hirundo.popmap --variant_filter_expression 'DP>=3' --invariant_filter_expression 'DP>=3' --outfile_prefix ./pixy_results/pixy_100kb_$chrom
	rm -r pixy_zarr/$chrom/*
done
```

`sh pixyloop.sh ./processing_files/Hirundo_rustica_Barn2Flycatcher_ChromAssigned.list`

#### Concatenate results

First, make ordered list of populations in output:

`tail -n +2  ./pixy_results/pixy_100kb_QRBI01000187.1_pi.txt | awk '{print $1}' | uniq > parselist.pi`

The script below will concatenate pixy results.

concatenate_pixy_pi.sh:

```
parselist=$1
chromlist=$2
echo "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy_100kb_all_pi.txt
for pop in `cat $parselist`; do
	for chrom in `cat $chromlist`; do
		grep -w $pop ./pixy_results/pixy_100kb_${chrom}_pi.txt >> pixy_100kb_all_pi.txt
	done
done
```

`sh concatenate_pixy_pi.sh ./processing_files/parselist.pi Hirundo_rustica_Barn2Flycatcher_ChromAssigned.list`

The concatenated results of this script are in `pixy_results`.

### Tajima's *D* analysis

We'll use `vcftools` to estimate Tajima's *D* in sliding windows across the autosomes and Z chromosome.

#### Set up environment

```
cd ./analysis/
mkdir popgen_stats
cd popgen_stats
mkdir tajimaD
cd tajimaD
```

#### Make directory with sample lists for analysis

`mkdir sample_lists`

The directory contains lists of individuals per subspecies and hybrid zone (this subdirectory is in `processing_files`)

#### Estimate Tajima's *D* in sliding windows

The script below will calculate *D* in windows of a specified size. Give it a window-size abbreviation and missing data descriptor too.

window_tajimaD.sh:

```
vcf=$1
window=$2
abbrev=$3
miss=$4
for list in ./sample_lists/*.list; do
	name=`echo $list | cut -d'.' -f3`
	echo $name
	vcftools --gzvcf $vcf --keep $list --TajimaD $window --out $name.$miss.$abbrev
done
```

```
sh window_tajimaD.sh ../../../vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.vcf.gz 100000 100kb miss04
```

Results are in `tajima_d_results`.

### *Fst* and *PBS* analysis

We'll estimate relative population differentiation (*Fst*) using `vcftools`.

#### Set up environment

```
cd ./analysis/popgen_stats
mkdir fst
mkdir pbs
```

#### Estimate *Fst* in sliding windows

The script below will calculate *Fst* between pairs of populations in `processing_files/pairwise_fst.list`.
Feed it a window size and abbreviation too.

window_fst.sh:

```
list=$1
window=$2
abbrev=$3
for pair in `cat $list`; do
	pair1=`echo $pair | cut -d',' -f1`
	pair2=`echo $pair | cut -d',' -f2`
	name1=`echo $pair1 | cut -d'.' -f2`
	name2=`echo $pair2 | cut -d'.' -f2`
	echo calculating windowed Fst between $name1 and $name2
	vcftools --gzvcf ../../../vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.vcf.gz --weir-fst-pop ../../sample_lists/$pair1 --weir-fst-pop ../../sample_lists/$pair2 --fst-window-size $window --fst-window-step $window --out ./${name1}_${name2}.$abbrev
done
```

`sh window_fst.sh ./processing_files/pairwise_fst.list 100000 100kb`

Results are `fst_results`.

#### Calculate *PBS* from *Fst* results

We'll calculate *PBS* in sliding windows in the `R` script `Fst_PBS.R` (see [Analysis in R](#analysis-in-r) section below).

After calculating *PBS*, the scripts below will order scans by chromosome.

`cd ../pbs`

order_scans_pbs.py:

```
import sys

out = open(sys.argv[3],'w')
out.write('CHROM'+'\t'+'SCAFF'+'\t'+'BIN_START'+'\t'+'BIN_END'+'\t'+'pbs'+'\n')

chrom_name = sys.argv[1].split('list.')[1]
chrom_name = chrom_name.split('.txt')[0]

for line in open(sys.argv[1],'r'):
	chrom = line.split()[0]
	matches = []
	order = []
	with open(sys.argv[2],'r') as pbs:
		next(pbs)
		for p in pbs:
			scaff = p.split()[0]
			start = p.split()[1]
			if str(scaff) == str(chrom):
				matches.append(p)
				order.append(int(start))
	sort_order = sorted(order)
	for o in sort_order:
		for m in matches:
			if int(m.split()[1]) == int(o):
				scaff = m.split()[0]
				start = m.split()[1]
				end = m.split()[2]
				pbs = m.split()[14]
				out.write(str(chrom_name)+'\t'+str(scaff)+'\t'+str(start)+'\t'+str(end)+'\t'+str(pbs)+'\n')
```

order_chrom_pbs.sh (wrapper to run `order_scans_pbs.py` per subspecies):

```
for chrom in `cat ../chrom.list`; do
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.er.data.txt pbs.er.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.gu.data.txt pbs.gu.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.ru.data.txt pbs.ru.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.sa.data.txt pbs.sa.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.tr.data.txt pbs.tr.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.ty.data.txt pbs.ty.$chrom.txt
done
```

`sh order_chrom_pbs.sh`


### ADMIXTURE analysis

We'll run `ADMIXTURE` for a series of *K* genetic clusters on Z-linked and autosomal SNPs.

#### Set up environment

```
cd ./analysis/
mkdir admixture
cd admixture
mkdir analysis_chrZ
mkdir anlaysis_auto
mkdir input
```

#### Convert SNPs in VCF to .ped format with `plink`

```
cd input
plink --vcf ../../../vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ.vcf --make-bed --out hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ --allow-extra-chr --recode12
plink --vcf ../../../vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto.vcf --make-bed --out hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto --allow-extra-chr --recode12
```

#### Fix scaffold names in .map and .bim files so `ADMIXTURE` doesn't barf

`ADMIXTURE` wants scaffolds to be represented as integers, and will fail if bootstrapping is used on non-integer .map and .bim input.

```
sed -i.bak -e 's/QRBI010000//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ.bim 
sed -i.bak -e 's/QRBI01000//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ.bim 
sed -i.bak -e 's/\.1//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ.bim 

sed -i.bak -e 's/QRBI010000//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ.map 
sed -i.bak -e 's/QRBI01000//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ.map
sed -i.bak -e 's/\.1//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ.map

sed -i.bak -e 's/QRBI010000//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto.bim 
sed -i.bak -e 's/QRBI01000//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto.bim 
sed -i.bak -e 's/\.1//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto.bim 

sed -i.bak -e 's/QRBI010000//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto.map 
sed -i.bak -e 's/QRBI01000//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto.map
sed -i.bak -e 's/\.1//g' hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto.map
```

#### Run ADMIXTURE

The script below will run `ADMIXTURE` for a series of *K* values:

run_admixture.sh:

```
ped=$1
bootstraps=$2
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
	admixture --cv --B$bootstraps $ped $K | tee log${K}.out
done
```

Copy script to analysis directories:

```
cp run_admixture.sh analysis_chrZ
cp run_admixture.sh analysis_auto
```

Perform analysis:

```
cd analysis_auto
sh run_admixture.sh ../input/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.auto.ped 200
cd ../analysis_chrZ
sh run_admixture.sh ../input/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.snps.miss04.maf05.thin100bp.ingroup.chrZ.ped 200
cd ..
```

#### Evaluate CV error

```
grep -h CV analysis_chrZ/log*.out
grep -h CV analysis_auto/log*.out
```

### Topology weighting analysis

We'll use the topology weighting procedure in `TWISST` to measure support for alternative triplet topologies across the genome.

#### The triplets

Each analysis will use *Hirundo smithii* as an outgroup. There are two 'shallow' and two 'deep' timescale analyses.

Shallow triplets:
1. *rustica*, *savignii*, and *transitiva*
2. *erythrogaster*, *gutturalis*, and *tytleri*

Deep triplets:
1. *savignii*, *transitiva*, and *erythrogaster*
2. *erythrogaster*, *tytleri*, and *savignii*

Various sample lists and popmaps for topology weighting analyses are in `processing_files/twisst`.

