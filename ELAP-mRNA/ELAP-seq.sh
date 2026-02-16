#!/bin/bash

###
### ELAP-seq.sh 
###    - input bam files for IP and input samples, and the regions covered by IP peaks
### 
###  
### Usage:
###   ELAP-seq.sh <rep1-III-IP.bam> <rep1-III-IV-IP.bam> <rep1-III-input.bam> <rep1-III-IV-input.bam> <rep1-III-peaks.bed> <rep1-III-IV-peaks.bed> <rep1-name> 
###   ELAP-seq.sh <rep2-III-IP.bam> <rep2-III-IV-IP.bam> <rep2-III-input.bam> <rep2-III-IV-input.bam> <rep2-III-peaks.bed> <rep2-III-IV-peaks.bed> <rep2-name> 
### Options:
###   <input.bam>         Input bam file to read. 
###   <IP.bam>         IP bam file to read
###   <peaks.bed>      IP peaks corresonding to the bam
###   <name>            given to specifiy the replicate
###
### Output:
###  chr start	end     strand    base	rep1_III/IV_In_stop	rep1_III/IV_In_sum  	rep1_III/IV_IP_sum	rep1_In_stop	rep1_IP_stop	rep1_In_sum rep1_IP_sum	peak	samp
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###
mkdir $7
bash arrest.sh $3 $1 $5 ./$7/$7-III-inside.out ./$7/$7-III-outside.out
bash arrest.sh $4 $2 $6 ./$7/$7-III-IV-inside.out ./$7/$7-III-IV-outside.out

bash calculate_III.sh $1 ./$7/$7-III-inside.out ./$7/$7-III-outside.out ./$7/$7-III-inside-unfiltered.bed ./$7/$7-III-outside-unfiltered.bed ./$7/$7-III-in ./$7/$7-III-out $7-III-unfiltered.bed
bash calculate_III_IV.sh $2 ./$7/$7-III-IV-inside.out ./$7/$7-III-IV-outside.out ./$7/$7-III-IV-inside-unfiltered.bed ./$7/$7-III-IV-outside-unfiltered.bed ./$7/$7-III-IV-in ./$7/$7-III-IV-out ./$7/$7-III-IV-unfiltered.bed

sort -k1,1 -k2,2n -k13,13 ./$7/$7-III-inside-unfiltered.bed > ./$7/$7-III-inside-unfiltered-sort.bed
sort -k1,1 -k2,2n -k13,13 ./$7/$7-III-outside-unfiltered.bed > ./$7/$7-III-outside-unfiltered-sort.bed
sort -k1,1 -k2,2n -k13,13 ./$7/$7-III-IV-inside-unfiltered.bed > ./$7/$7-III-IV-inside-unfiltered-sort.bed
sort -k1,1 -k2,2n -k13,13 ./$7/$7-III-IV-outside-unfiltered.bed > ./$7/$7-III-IV-outside-unfiltered-sort.bed

python3 Rm_bg_1.py ./$7/$7-III-inside-unfiltered-sort.bed | tr ' ' '\t' > ./$7/$7-III-inside-unfiltered-ab.bed
python3 Rm_bg_1.py ./$7/$7-III-outside-unfiltered-sort.bed | tr ' ' '\t' > ./$7/$7-III-outside-unfiltered-ab.bed
python3 Rm_bg_1.py ./$7/$7-III-IV-inside-unfiltered-sort.bed | tr ' ' '\t' > ./$7/$7-III-IV-inside-unfiltered-ab.bed
python3 Rm_bg_1.py ./$7/$7-III-IV-outside-unfiltered-sort.bed | tr ' ' '\t' > ./$7/$7-III-IV-outside-unfiltered-ab.bed


python3 Rm_bg_2.py ./$7/$7-III-inside-unfiltered-ab.bed | tr ' ' '\t' > ./$7/$7-III-inside-block.bed
python3 Rm_bg_2.py ./$7/$7-III-outside-unfiltered-ab.bed | tr ' ' '\t' > ./$7/$7-III-outside-block.bed
python3 Rm_bg_2.py ./$7/$7-III-IV-inside-unfiltered-ab.bed | tr ' ' '\t' > ./$7/$7-III-IV-inside-block.bed
python3 Rm_bg_2.py ./$7/$7-III-IV-outside-unfiltered-ab.bed | tr ' ' '\t' > ./$7/$7-III-IV-outside-block.bed


bedtools subtract -a ./$7/$7-III-inside-unfiltered-ab.bed -b ./$7/$7-III-inside-block.bed > ./$7/$7-III-inside-unfiltered-2.bed
bedtools subtract -a ./$7/$7-III-outside-unfiltered-ab.bed -b ./$7/$7-III-outside-block.bed > ./$7/$7-III-outside-unfiltered-2.bed
bedtools subtract -a ./$7/$7-III-IV-inside-unfiltered-ab.bed -b ./$7/$7-III-IV-inside-block.bed > ./$7/$7-III-IV-inside-unfiltered-2.bed
bedtools subtract -a ./$7/$7-III-IV-outside-unfiltered-ab.bed -b ./$7/$7-III-IV-outside-block.bed > ./$7/$7-III-IV-outside-unfiltered-2.bed


python3 Rm_bg_3.py ./$7/$7-III-inside-unfiltered-2.bed > ./$7/$7-III-inside-unfiltered-low.bed
python3 Rm_bg_3.py ./$7/$7-III-outside-unfiltered-2.bed > ./$7/$7-III-outside-unfiltered-low.bed
python3 Rm_bg_3.py ./$7/$7-III-IV-inside-unfiltered-2.bed > ./$7/$7-III-IV-inside-unfiltered-low.bed
python3 Rg_bg_3.py ./$7/$7-III-IV-outside-unfiltered-2.bed > ./$7/$7-III-IV-outside-unfiltered-low.bed


awk '!visited[$0]++' ./$7/$7-III-inside-unfiltered-low.bed > ./$7/$7-III-inside-unfiltered-low-1.bed
awk '!visited[$0]++' ./$7/$7-III-outside-unfiltered-low.bed > ./$7/$7-III-outside-unfiltered-low-1.bed
awk '!visited[$0]++' ./$7/$7-III-IV-inside-unfiltered-low.bed > ./$7/$7-III-IV-inside-unfiltered-low-1.bed
awk '!visited[$0]++' ./$7/$7-III-IV-outside-unfiltered-low.bed > ./$7/$7-III-IV-outside-unfiltered-low-1.bed


cat ./$7/$7-III-inside-unfiltered-low-1.bed | tr ' ' '\t' > ./$7/$7-III-inside-unfiltered-low-2.bed
cat ./$7/$7-III-outside-unfiltered-low-1.bed | tr ' ' '\t' > ./$7/$7-III-outside-unfiltered-low-2.bed
cat ./$7/$7-III-IV-inside-unfiltered-low-1.bed | tr ' ' '\t' > ./$7/$7-III-IV-inside-unfiltered-low-2.bed
cat ./$7/$7-III-IV-outside-unfiltered-low-1.bed | tr ' ' '\t' > ./$7/$7-III-IV-outside-unfiltered-low-2.bed

bedtools subtract -a ./$7/$7-III-inside-unfiltered-2.bed -b ./$7/$7-III-inside-unfiltered-low-2.bed > ./$7/$7-III-inside-unfiltered-3.bed
bedtools subtract -a ./$7/$7-III-outside-unfiltered-2.bed -b ./$7/$7-III-outside-unfiltered-low-2.bed > ./$7/$7-III-outside-unfiltered-3.bed
bedtools subtract -a ./$7/$7-III-IV-inside-unfiltered-2.bed -b ./$7/$7-III-IV-inside-unfiltered-low-2.bed > ./$7/$7-III-IV-inside-unfiltered-3.bed
bedtools subtract -a ./$7/$7-III-IV-outside-unfiltered-2.bed -b ./$7/$7-III-IV-outside-unfiltered-low-2.bed > ./$7/$7-III-IV-outside-unfiltered-3.bed

bash filter_III.sh $1 ./$7/$7-III-inside-unfiltered-3.bed ./$7/$7-III-outside-unfiltered-3.bed ./$7/$7-III-unique.bed
bash filter_III_IV.sh $2 ./$7/$7-III-IV-inside-unfiltered-3.bed ./$7/$7-III-IV-outside-unfiltered-3.bed ./$7/$7-III-IV-unique.bed

sort -k1,1 -k2,2n ./$7/$7-III-unique.bed > ./$7/$7-III-unique-1.bed
sort -k1,1 -k2,2n ./$7/$7-III-IV-unique.bed > ./$7/$7-III-IV-unique-1.bed


python3 Stutter1.py ./$7/$7-III-unique-1.bed | tr ' ' '\t' > ./$7/$7-III-stutter-filter.bed
python3 Stutter1.py ./$7/$7-III-IV-unique-1.bed | tr ' ' '\t' > ./$7/$7-III-IV-stutter-filter.bed


python3 Stutter2.py ./$7/$7-III-stutter-filter.bed | tr ' ' '\t' > ./$7/$7-III-remove.bed
python3 Stutter2.py ./$7/$7-III-IV-stutter-filter.bed | tr ' ' '\t' > ./$7/$7-III-IV-remove.bed


bedtools subtract -a ./$7/$7-III-stutter-filter.bed -b ./$7/$7-III-remove.bed > ./$7/$7-III-stutter-filter-2.bed
bedtools subtract -a ./$7/$7-III-IV-stutter-filter.bed -b ./$7/$7-III-IV-remove.bed > ./$7/$7-III-IV-stutter-filter-2.bed

awk '{ if($10 >4) print $0;}' ./$7/$7-III-stutter-filter-2.bed > ./$7/$7-III-filter1.bed
awk '{ if($10 >4) print $0;}' ./$7/$7-III-IV-stutter-filter-2.bed > ./$7/$7-III-IV-filter1.bed

awk '($14 <= 0.1) || ($8 < 3) || ($15 / $14 >= 3)' ./$7/$7-III-filter1.bed > ./$7/$7-III-filter2.bed
awk '($14 <= 0.1) || ($8 < 3) || ($15 / $14 >= 3)' ./$7/$7-III-IV-filter1.bed > ./$7/$7-III-IV-filter2.bed


bedtools subtract -a ./$7/$7-III-IV-filter-2.bed -b ./$7/$7-III-filter-2.bed > ./$7/new.bed
cat ./$7/$7-III-filter-2.bed ./$7/new.bed | sort -k1,1 > ./$7/$7-combined.bed
