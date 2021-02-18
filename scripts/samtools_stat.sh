for i in *.bam; do
	echo calculating mapping statistics for $i
	samtools stat -@ 8 $i > ./stats/$i.stat.txt
done
