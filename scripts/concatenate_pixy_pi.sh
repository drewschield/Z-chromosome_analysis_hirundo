parselist=$1
chromlist=$2
echo "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing" > pixy_100kb_all_pi.txt
for pop in `cat $parselist`; do
	for chrom in `cat $chromlist`; do
		grep -w $pop ./pixy_results/pixy_100kb_${chrom}_pi.txt >> pixy_100kb_all_pi.txt
	done
done
