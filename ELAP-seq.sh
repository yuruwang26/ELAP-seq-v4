#!/bin/bash

###
### ELAP-seq.sh 
###    - input bam files for IP and input samples, and the regions covered by IP peaks
### 
###  
### Usage:
###   ELAP-seq.sh <rep1-III-IP.bam> <rep1-III-IV-IP.bam> <rep1-III-input.bam> <rep1-III-IV-input.bam> <rep1-peaks.bed> <rep1-name> <rep2-III-IP.bam> <rep2-III-IV-IP.bam> <rep2-III-input.bam> <rep2-III-IV-input.bam> <rep2-peaks.bed> <rep2-name> 
###
### Options:
###   <input.bam>         Input bam file to read. 
###   <IP.bam>         IP bam file to read
###   <peaks.bed>      IP peaks corresonding to the bam
###   <name>            given to specifiy the replicate
###   
###
### Output:
###  chr start	end     strand    base	rep1_III/IV_In_stop	rep1_III/IV_In_sum  	rep1_III/IV_IP_sum	rep1_In_stop	rep1_IP_stop	rep1_In_sum rep1_IP_sum	peak	samp	rep2_III/IV_In_stop	rep2_III/IV_In_sum	rep2_III/IV_IP_sum	rep2_In_stop	rep2_IP_stop	rep2_In_sum	rep2_IP_sum	peak	samp
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###

bash ELAP-seq-pre.sh $1 $2 $3 $4 $5 $6 
bash ELAP-seq-pre.sh $7 $8 $9 $10 $11 $12 
bedtools intersect -wa -wb -a ./$6/$6-combined-2.bed -b ./$12/$12-combined-2.bed > HeLa.bed
awk '!visited[$0]++' HeLa.bed | awk '{print $1,$2,$3,$4,$5,$8,$6,$7,$11,$12,$9,$10,$13,$14,$22,$20,$21,$25,$26,$23,$24,$27,$28}' | awk -v OFS="\t" '{$1=$1; print}' | tr ' ' '\t' | sort -k1,1 -k2,2n > HeLa-sort.bed
awk '{print $0"\t"($9+$18)/2}' HeLa-sort.bed > HeLa-input-avg.bed
awk '{print $0"\t"($10+$19)/2}' HeLa-input-avg.bed > HeLa-IP-avg.bed
python3 In_stop.py HeLa-IP-avg.bed > HeLa-filter-2.bed
awk '{ if($5 == "T") print $0;}' HeLa-filter-2.bed > HeLa-filter-T.bed
