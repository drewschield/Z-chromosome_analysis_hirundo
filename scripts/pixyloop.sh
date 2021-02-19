list=$1
for chrom in `cat $list`; do
	pixy --stats pi --vcf ../../vcf/vcf_all-sites_chrom/hirundo_rustica+smithii.allsites.HardFilter.recode.depth.chrom.final.$chrom.vcf.gz --zarr_path pixy_zarr --window_size 100000 --populations hirundo.popmap --variant_filter_expression 'DP>=3' --invariant_filter_expression 'DP>=3' --outfile_prefix ./pixy_results/pixy_100kb_$chrom
	rm -r pixy_zarr/$chrom/*
done
