#!/bin/bash

###
### ELAP-seq.sh 
###    - input bam files for IP and input samples, and the regions covered by IP peaks
### 
###  
### Usage:
###   ELAP-seq.sh <III-IP.bam> <III-IV-IP.bam> <III-input.bam> <III-IV-input.bam> <peaks.bed> <name> 
###
### Options:
###   <input.bam>         Input bam file to read. 
###   <IP.bam>         IP bam file to read
###   <peaks.bed>      IP peaks corresonding to the bam
###   <name>            given to specifiy the replicate
###
### Output:
###  chr start	end     strand    base	rep1_III/IV_In_stop	rep1_III/IV_In_sum  	rep1_III/IV_IP_sum	rep1_In_stop	rep1_IP_stop	rep1_In_sum rep1_IP_sum	peak	samp	rep2_III/IV_In_stop	rep2_III/IV_In_sum	rep2_III/IV_IP_sum	rep2_In_stop	rep2_IP_stop	rep2_In_sum	rep2_IP_sum	peak	samp
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###
mkdir $6
bash arrest.sh $3 $1 $5 ./$6/$6-III-inside.out ./$6/$6-III-outside.out
bash arrest.sh $4 $2 $5 ./$6/$6-III-IV-inside.out ./$6/$6-III-IV-outside.out

bash calculate1.sh $1 ./$6/$6-III-inside.out ./$6/$6-III-outside.out ./$6/$6-III-inside-unfiltered.bed ./$6/$6-III-outside-unfiltered.bed ./$6/$6-III-in ./$6/$6-III-out $6-III-unfiltered.bed
bash calculate2.sh $2 ./$6/$6-III-IV-inside.out ./$6/$6-III-IV-outside.out ./$6/$6-III-IV-inside-unfiltered.bed ./$6/$6-III-IV-outside-unfiltered.bed ./$6/$6-III-IV-in $6-III-IV-out ./$6/$6-III-IV-unfiltered.bed
python3 Rm_bg_1.py ./$6/$6-III-inside-unfiltered.bed | tr ' ' '\t' > ./$6/$6-III-inside-unfiltered-ab.bed
python3 Rm_bg_1.py ./$6/$6-III-outside-unfiltered.bed | tr ' ' '\t' > ./$6/$6-III-outside-unfiltered-ab.bed
python3 Rm_bg_1.py ./$6/$6-III-IV-inside-unfiltered.bed | tr ' ' '\t' > ./$6/$6-III-IV-inside-unfiltered-ab.bed
python3 Rm_bg_1.py ./$6/$6-III-IV-outside-unfiltered.bed | tr ' ' '\t' > ./$6/$6-III-IV-outside-unfiltered-ab.bed


python3 Rm_bg_2.py ./$6/$6-III-inside-unfiltered-ab.bed | tr ' ' '\t' > ./$6/$6-III-inside-block.bed
python3 Rm_bg_2.py ./$6/$6-III-outside-unfiltered-ab.bed | tr ' ' '\t' > ./$6/$6-III-outside-block.bed
python3 Rm_bg_2.py ./$6/$6-III-IV-inside-unfiltered-ab.bed | tr ' ' '\t' > ./$6/$6-III-IV-inside-block.bed
python3 Rm_bg_2.py ./$6/$6-III-IV-outside-unfiltered-ab.bed | tr ' ' '\t' > ./$6/$6-III-IV-outside-block.bed


bedtools subtract -a ./$6/$6-III-inside-unfiltered-ab.bed -b ./$6/$6-III-inside-block.bed > ./$6/$6-III-inside-unfiltered-2.bed
bedtools subtract -a ./$6/$6-III-outside-unfiltered-ab.bed -b ./$6/$6-III-outside-block.bed > ./$6/$6-III-outside-unfiltered-2.bed
bedtools subtract -a ./$6/$6-III-IV-inside-unfiltered-ab.bed -b ./$6/$6-III-IV-inside-block.bed > ./$6/$6-III-IV-inside-unfiltered-2.bed
bedtools subtract -a ./$6/$6-III-IV-outside-unfiltered-ab.bed -b ./$6/$6-III-IV-outside-block.bed > ./$6/$6-III-IV-outside-unfiltered-2.bed


python3 Rm_bg_3.py ./$6/$6-III-inside-unfiltered-2.bed > ./$6/$6-III-inside-unfiltered-low.bed
python3 Rm_bg_3.py ./$6/$6-III-outside-unfiltered-2.bed > ./$6/$6-III-outside-unfiltered-low.bed
python3 Rm_bg_3.py ./$6/$6-III-IV-inside-unfiltered-2.bed > ./$6/$6-III-IV-inside-unfiltered-low.bed
python3 Rg_bg_3.py ./$6/$6-III-IV-outside-unfiltered-2.bed > ./$6/$6-III-IV-outside-unfiltered-low.bed


awk '!visited[$0]++' ./$6/$6-III-inside-unfiltered-low.bed > ./$6/$6-III-inside-unfiltered-low-1.bed
awk '!visited[$0]++' ./$6/$6-III-outside-unfiltered-low.bed > ./$6/$6-III-outside-unfiltered-low-1.bed
awk '!visited[$0]++' ./$6/$6-III-IV-inside-unfiltered-low.bed > ./$6/$6-III-IV-inside-unfiltered-low-1.bed
awk '!visited[$0]++' ./$6/$6-III-IV-outside-unfiltered-low.bed > ./$6/$6-III-IV-outside-unfiltered-low-1.bed


cat ./$6/$6-III-inside-unfiltered-low-1.bed | tr ' ' '\t' > ./$6/$6-III-inside-unfiltered-low-2.bed
cat ./$6/$6-III-outside-unfiltered-low-1.bed | tr ' ' '\t' > ./$6/$6-III-outside-unfiltered-low-2.bed
cat ./$6/$6-III-IV-inside-unfiltered-low-1.bed | tr ' ' '\t' > ./$6/$6-III-IV-inside-unfiltered-low-2.bed
cat ./$6/$6-III-IV-outside-unfiltered-low-1.bed | tr ' ' '\t' > ./$6/$6-III-IV-outside-unfiltered-low-2.bed

bedtools subtract -a ./$6/$6-III-inside-unfiltered-2.bed -b ./$6/$6-III-inside-unfiltered-low-2.bed > ./$6/$6-III-inside-unfiltered-3.bed
bedtools subtract -a ./$6/$6-III-outside-unfiltered-2.bed -b ./$6/$6-III-outside-unfiltered-low-2.bed > ./$6/$6-III-outside-unfiltered-3.bed
bedtools subtract -a ./$6/$6-III-IV-inside-unfiltered-2.bed -b ./$6/$6-III-IV-inside-unfiltered-low-2.bed > ./$6/$6-III-IV-inside-unfiltered-3.bed
bedtools subtract -a ./$6/$6-III-IV-outside-unfiltered-2.bed -b ./$6/$6-III-IV-outside-unfiltered-low-2.bed > ./$6/$6-III-IV-outside-unfiltered-3.bed

cat ./$6/$6-III-IV-inside-unfiltered-2.bed ./$6/$6-III-IV-outside-unfiltered-2.bed > ./$6/$6-III-IV-unfiltered-2.bed

bash filter.sh $1 ./$6/$6-III-inside-unfiltered-3.bed ./$6/$6-III-outside-unfiltered-3.bed ./$6/$6-III-unique.bed
bash filter.sh $2 ./$6/$6-III-IV-inside-unfiltered-3.bed ./$6/$6-III-IV-outside-unfiltered-3.bed ./$6/$6-III-IV-unique.bed

sort -k1,1 -k2,2n ./$6/$6-III-unique.bed > ./$6/$6-III-unique-1.bed
sort -k1,1 -k2,2n ./$6/$6-III-IV-unique.bed > ./$6/$6-III-IV-unique-1.bed


python3 Stutter1.py ./$6/$6-III-unique-1.bed | tr ' ' '\t' > ./$6/$6-III-stutter-filter.bed
python3 Stutter1.py ./$6/$6-III-IV-unique-1.bed | tr ' ' '\t' > ./$6/$6-III-IV-stutter-filter.bed


python3 Stutter2.py ./$6/$6-III-stutter-filter.bed | tr ' ' '\t' > ./$6/$6-III-remove.bed
python3 Stutter2.py ./$6/$6-III-IV-stutter-filter.bed | tr ' ' '\t' > ./$6/$6-III-IV-remove.bed


bedtools subtract -a ./$6/$6-III-stutter-filter.bed -b ./$6/$6-III-remove.bed > ./$6/$6-III-stutter-filter-2.bed
bedtools subtract -a ./$6/$6-III-IV-stutter-filter.bed -b ./$6/$6-III-IV-remove.bed > ./$6/$6-III-IV-stutter-filter-2.bed
bedtools subtract -a ./$6/$6-III-IV-stutter-filter-2.bed -b ./$6/$6-III-stutter-filter-2.bed > ./$6/$6-new.bed

cat ./$6/$6-III-stutter-filter-2.bed ./$6/$6-new.bed | sort -k1,1 > ./$6/$6-combined.bed

bedtools intersect -wa -wb -a ./$6/$6-III-IV-unfiltered-2.bed -b ./$6/$6-combined.bed > ./$6/$6-combined-1.bed
awk '{ if($10 >2 && $13 > 9) print $0;}' HeLa-rep1-combined-1.bed > HeLa-rep1-combined-filtered.bed
awk '{print $1,$2,$3,$5,$7,$12,$13,$14,$29,$30,$31,$32,$33,$34}' ./$6/$6-combined-filtered.bed | awk -v OFS="\t" '$1=$1' > ./$6/$6-combined-2.bed
