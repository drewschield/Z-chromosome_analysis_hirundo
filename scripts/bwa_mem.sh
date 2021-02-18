for line in `cat sample.list`; do
	name=$line
	out=`echo $name | sed 's/L_/HR/'`
	echo "Mapping filtered $out data to reference"
	bwa mem -t 16 -R "@RG\tID:$out\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:$out" Hirundo_rustica_Chelidonia.fasta ./fastq_filtered/${name}_1_P.trim.fq.gz ./fastq_filtered/${name}_2_P.trim.fq.gz | samtools sort -@ 16 -O bam -T temp -o ./bam/$out.bam -
done
