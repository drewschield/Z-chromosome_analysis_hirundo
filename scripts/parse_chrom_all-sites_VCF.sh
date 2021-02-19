for scaff in `cat Hirundo_rustica_Barn2Flycatcher_ChromAssigned.list`; do
	echo parsing $scaff VCF
	bcftools view --threads 8 -r $scaff -O z -o ./vcf_all-sites_chrom/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.$scaff.vcf.gz hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.vcf.gz
	tabix -p vcf ./vcf_all-sites_chrom/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.$scaff.vcf.gz
	echo indexed $scaff VCF
done
