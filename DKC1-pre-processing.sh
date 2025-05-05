#!/bin/bash


samples="Input-sictrlD-rep1 Input-sictrlD-rep2 Input-siDKC1-rep1 Input-siDKC1-rep2 Input-sictrlT-rep1 Input-sictrlT-rep2 Input-siTruB1-rep1 Input-siTruB1-rep2 adp-sictrlD-rep1 adp-sictrlD-rep2 adp-siDKC1-rep1 adp-siDKC1-rep2 adp-sictrlT-rep1 adp-sictrlT-rep2 adp-siTruB1-rep1 adp-siTruB1-rep2 IP-sictrlD-rep1 IP-sictrlD-rep2 IP-siDKC1-rep1 IP-siDKC1-rep2 IP-sictrlT-rep1 IP-sictrlT-rep2 IP-siTruB1-rep1 IP-siTruB1-rep2" 
 
for s in $samples

do
~/Tools/anaconda3/bin/clumpify.sh in=./sorted/$s.fastq out=./sorted/$s-dedupe.fastq >>duplicates_removal_1.log

cutadapt -u 8 -o ./sorted/$s-tmp.trimmed.fastq ./sorted/$s-dedupe.fastq
cutadapt -u -5 -q 10,10 -m 19 -o ./sorted/$s-trim.fastq ./sorted/$s-tmp.trimmed.fastq
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file align_summary -p 4 -U ./sorted/$s-trim.fastq |samtools view -bS |samtools sort -o ./sorted/$s.bam
samtools index ./sorted/$s.bam ./sorted/$s.bai

done
