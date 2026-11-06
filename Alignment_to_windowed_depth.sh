#!/bin/sh

samtools faidx genome.fasta
bedtools makewindows -g genome.fasta.fai -w 50000 | sort -k1,1 -k2,2n > genome_windows_sorted_50k.bed

### Pass the bam files
while [ $# -gt 0 ]
do
  	INNAME=`echo $1 | sed "s/.bam//"`

	bedtools bamtobed -i ${INNAME}.bam > ${INNAME}.bed
	sort ${INNAME}.bed -k1,1 -k2,2n > ${INNAME}_sorted.bed

	awk '$5 == 60' ${INNAME}_sorted.bed > ${INNAME}_mapq60_sorted.bed
        bedtools coverage -sorted -a genome_windows_sorted_50k.bed -b ${INNAME}_mapq60_sorted.bed > ${INNAME}_mapq60_sorted_coverage.txt
	awk '{ if (!($1 in mapping)) {
		count++; mapping[$1] = count;
        	} print mapping[$1], $0;}' ${INNAME}_mapq60_sorted_coverage.txt > ${INNAME}_mapq60_sorted_coverage_clean.txt
        [ -f ${INNAME}_mapq60_sorted_coverage_clean.txt ] && rm ${INNAME}_mapq60_sorted_coverage.txt

        shift
done
