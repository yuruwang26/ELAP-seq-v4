#!/bin/bash


samples="HEK-sictrl-input-IV-rep1 HEK-sictrl-input-IV-rep2 HEK-siDKC1-input-IV-rep1 HEK-siDKC1-input-IV-rep2 HEK-sictrl-IP-IV-rep1 HEK-sictrl-IP-IV-rep2 HEK-siDKC1-IP-IV-rep1 HEK-siDKC1-IP-IV-rep2" 
 
for s in $samples

do
~/Tools/anaconda3/bin/clumpify.sh in=$s.fastq out=$s-dedupe.fastq >>duplicates_removal_1.log

cutadapt -u 8 -o $s-tmp.trimmed.fastq $s-dedupe.fastq
cutadapt -u -5 -q 10,10 -m 19 -o $s-trim.fastq $s-tmp.trimmed.fastq
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file align_summary -p 4 -U $s-trim.fastq |samtools view -bS |samtools sort -o .$s.bam
samtools index $s.bam $s.bai

done
