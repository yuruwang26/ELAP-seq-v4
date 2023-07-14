#!/bin/bash

###
### filter1.sh 
###    -  filter sites from library III based on stop ratio and coverage, and select for sites covered by at least 1 uniquely mapped reads
### 
###  
### Usage:
###   filter1.sh <IP.bam>  <inside-unfiltered.out> <outside-unfiltered.out> <unique.bed>
###
### Options:
###   
###   <IP.bam>         IP bam file to read
###   <inside-unfiltered.out>      sites inside IP peaks after removing background
###   <outside-unfiltered.out>     sites outside IP peaks after removing background
###   <unique.bed>           sites filtered based on stop ratio and coverage and covered by at least 1 uniquely mapped reads
###
### Output:
###  chr start end p-value strand arrest_score base arrest_input readthrough_input arrest_IP readthrough_IP input_sum IP_sum arrest_ratio_input arrest_ratio_IP peak sample
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###

awk '{ if(($15 >= 0.8 && $13 > 5 ) || ($15 > 0.2 && $15 < 0.8 && $13 > 7) || ($10 > 10 && $15 > 0.10 && $15 < 0.2 )) print $0;}' $2 > 6-18.out

awk '{ if(($15 >= 0.8 && $13 > 4 ) || ($15 >= 0.3 && $15 < 0.8 && $13 > 6) || ($10 > 1 && $15 >= 0.1 && $15 < 0.3)) print $0;}' $3 > 7-18.out

awk -v OFS="\t" '$1=$1' 6-18.out > 6-18.bed
awk -v OFS="\t" '$1=$1' 7-18.out > 7-18.bed
cat 6-18.bed 7-18.bed > all.bed
sort -k1,1 all.bed > all-1.bed
rm 6-18.out 6-18.bed 7-18.out 7-18.bed

bedtools bamtobed -split -i $1 > coverage.bed
awk '{ if($5 == 60) print $0;}' coverage.bed | awk '{print $1,$2,$3}' OFS="\t" | sort -k1,1 -k2,2n | bedtools merge -i > MQ-60-final.bed

bedtools intersect -a all-1.bed -b MQ-60-final.bed > $4

rm all.bed all-1.bed MQ-60-final.bed
