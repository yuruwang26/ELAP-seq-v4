#!/bin/bash

###
### calculate1.sh 
###    - calculate sum total reads for input and IP, and calculate stop ratios in input and IP
###  
### Usage:
###   calculate1.sh  <output.out> <calculated.out> 
###
### Options:
###   <output.out>     the ourput file generated from DKC1-arrest.sh
###   <calculated.out>    Sites calculated
### Output:
###  chr start end p-value strand arrest_score base arrest_input readthrough_input arrest_IP readthrough_IP input_sum IP_sum arrest_ratio_input arrest_ratio_IP peak sample
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###
awk '{ if($10 == "*") $10="0,0,0,0";}1' OFS="\t" $1 > $1-2.out
awk '{ if($9 == "*") $9="0,0,0,0";}1' OFS="\t" $1-2.out > $1-3.out
awk '{ if($8 == "*") $8="0,0,0,0";}1' OFS="\t" $1-3.out > $1-4.out
awk '{ if($7 == "*") $7="0,0,0,0";}1' OFS="\t" $1-4.out > $1-5.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$10)}1' $1-5.out > $1-6.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$9)}1' $1-6.out > $1-7.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$8)}1' $1-7.out > $1-8.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$7)}1' $1-8.out > $1-9.out
awk '{print $0"\t"($7+$8+$9+$10)}' $1-9.out >  $1-10.out
awk '{print $0"\t"($11+$12+$13+$14)}' $1-10.out >  $1-11.out
awk '{print $0"\t"($15+$16+$17+$18)}' $1-11.out >  $1-12.out
awk '{print $0"\t"($19+$20+$21+$22)}' $1-12.out >  $1-13.out
awk '{print $0"\t"($26+$27)}' $1-13.out >  $1-14.out
awk '{print $0"\t"($28+$29)}' $1-14.out >  $1-15.out


awk '{$4=$7=$8=$9=$10=$11=$12=$13=$14=$15=$16=$17=$18=$19=$20=$21=$22=""; print $0}' $1-15.out > $2

rm $1-2.out $1-3.out $1-4.out $1-5.out $1-6.out $1-7.out $1-8.out $1-9.out $1-10.out $1-11.out $1-12.out $1-13.out $1-14.out $1-15.out 

