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
