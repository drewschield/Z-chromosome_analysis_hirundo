for chrom in `cat ../chrom.list`; do
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.er.data.txt pbs.er.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.gu.data.txt pbs.gu.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.ru.data.txt pbs.ru.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.sa.data.txt pbs.sa.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.tr.data.txt pbs.tr.$chrom.txt
	python order_scans_pbs.py ../hirundo_rustica_scaffold_list.$chrom.txt pbs.ty.data.txt pbs.ty.$chrom.txt
done
