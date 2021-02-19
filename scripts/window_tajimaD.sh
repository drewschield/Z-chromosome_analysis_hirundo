vcf=$1
window=$2
abbrev=$3
miss=$4
for list in ./sample_lists/*.list; do
	name=`echo $list | cut -d'.' -f3`
	echo $name
	vcftools --gzvcf $vcf --keep $list --TajimaD $window --out $name.$miss.$abbrev
done
