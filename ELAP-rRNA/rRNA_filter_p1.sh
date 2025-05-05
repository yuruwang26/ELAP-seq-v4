#!/bin/bash

###
### rRNA_filter_p1.sh 
###    - Calculate stop ratios and coverage
### 
###  
### Usage:
###   rRNA_filter_p1.sh <input.out> <output.bed> 
###
### Options:
###   <input.out>         arrested sites called by JACUSA. 
###   <output.bed>         sites with stop ratios and reads coverage calculated
###
### Output:
###  chr start end p-value strand arrest_score base arrest_input readthrough_input arrest_IP readthrough_IP input_sum IP_sum input_stop_ratio IP_stop_ratio 
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###


awk '{ if($9 != "*") print $0;}1' OFS="\t" $1 > 6-1.out
awk '{ if($10 == "*") $10="0,0,0,0";}1' OFS="\t" 6-1.out > 6-3.out
awk '{ if($8 == "*") $8="0,0,0,0";}1' OFS="\t" 6-3.out > 6-4.out
awk '{ if($7 == "*") $7="0,0,0,0";}1' OFS="\t" 6-4.out > 6-5.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$10)}1' 6-5.out > 6-6.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$9)}1' 6-6.out > 6-7.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$8)}1' 6-7.out > 6-8.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$7)}1' 6-8.out > 6-9.out
awk '{print $0"\t"($7+$8+$9+$10)}' 6-9.out >  6-10.out
awk '{print $0"\t"($11+$12+$13+$14)}' 6-10.out >  6-11.out
awk '{print $0"\t"($15+$16+$17+$18)}' 6-11.out >  6-12.out
awk '{print $0"\t"($19+$20+$21+$22)}' 6-12.out >  6-13.out
awk '{print $0"\t"($26+$27)}' 6-13.out >  6-14.out
awk '{print $0"\t"($28+$29)}' 6-14.out >  6-15.out

awk '{if ($31 > 0) {print $0"\t"$28/$31} }' 6-15.out >  6-16.out
awk '{if ($30 > 0) {print $0"\t"$26/$30} }' 6-16.out >  6-17.out
awk '{$4=$7=$8=$9=$10=$11=$12=$13=$14=$15=$16=$17=$18=$19=$20=$21=$22=$24=""; print $0}' 6-17.out > 6-18.out
awk -v OFS="\t" '$1=$1' 6-18.out > 6-19.out

awk '!visited[$0]++' 6-19.out > 6-20.out

cat 6-20.out | tr ' ' '\t' > $2

