#!/bin/bash

###
### rRNA_filter_p2.sh 
###    - remove sites whose stop ration in the IP sample is > 0.05
### 
###  
### Usage:
###   rRNA_filter_p2.sh <input.bed> <output.bed> 
###
### Options:
###   <input.bed>         Input file after remiving stutter sites
###   <outpust.bed>         output file with sites whose stop ratio in the IP sample is > 0.05
###
### Output:
###  chr start end p-value strand arrest_score base arrest_input readthrough_input arrest_IP readthrough_IP input_sum IP_sum input_stop_ratio IP_stop_ratio 
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###


awk '{ if($14 >= 0.05) print $0;}' $1 > 1.out

awk -v OFS="\t" '$1=$1' 1.out > 2.out

awk '!visited[$0]++' 2.out > 3.out

cat 3.out | tr ' ' '\t' > $2

rm 1.out 2.out 3.out 



