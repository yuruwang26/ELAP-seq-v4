#!/bin/bash

###
### ELAP-seq.sh 
###    - input bam files for IP and input samples, and the regions covered by IP peaks
### 
###  
### Usage:
###   ELAP-seq.sh <rep1-III-IP.bam> <rep1-III-IV-IP.bam> <rep1-III-input.bam> <rep1-III-IV-input.bam> <rep1-peaks.bed> <rep1-name> 
###   ELAP-seq.sh <rep2-III-IP.bam> <rep2-III-IV-IP.bam> <rep2-III-input.bam> <rep2-III-IV-input.bam> <rep2-peaks.bed> <rep2-name> 
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
bedtools intersect -wa -wb -a ./$6/$6-III-IV-unfiltered-2.bed -b ./$6/$6-combined.bed > ./$6/$6-combined-1.bed
awk '{print $1,$2,$3,$5,$7,$12,$13,$14,$29,$30,$31,$32,$33,$34}' ./$6/$6-combined-1.bed | awk -v OFS="\t" '$1=$1' > ./$6/$6-combined-2.bed


