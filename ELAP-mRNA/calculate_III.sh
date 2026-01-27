#!/bin/bash

###
### calculate_III.sh 
###    - calculate sum total reads for input and IP, calculate stop ratios in input and IP, and assign the originality of the site
### 
###  
### Usage:
###   calculate_III.sh <IP.bam>  <inside.out> <outside.out> <inside-unfiltered.bed> <outside-unfiltered.bed> <name-in> <name-out> <unfiltered>
###
### Options:
###   <IP.bam>         IP bam file to read
###   <inside.out>     Sites inside IP peaks
###   <outside.out>    Sites outside IP peaks
###   <inside-unfilter.out>            sites inside the IP peaks and stop ratio calculated
###   <outside-unfilter.out>           sites inside the IP peaks and stop ratio calculated
###   <name-in>       name for temporary storage   
###   <name-out>      name for temporary storage
###   <unfiltered>    all sites with stop ratio calculated
###
### Output:
###  chr start end p-value strand arrest_score base arrest_input readthrough_input arrest_IP readthrough_IP input_sum IP_sum arrest_ratio_input arrest_ratio_IP peak sample
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###

awk '{ if($9 != "*") print $0;}1' OFS="\t" $2 > $6-1.out
awk '{ if($10 == "*") $10="0,0,0,0";}1' OFS="\t" $6-1.out > $6-3.out
awk '{ if($8 == "*") $8="0,0,0,0";}1' OFS="\t" $6-3.out > $6-4.out
awk '{ if($7 == "*") $7="0,0,0,0";}1' OFS="\t" $6-4.out > $6-5.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$10)}1' $6-5.out > $6-6.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$9)}1' $6-6.out > $6-7.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$8)}1' $6-7.out > $6-8.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$7)}1' $6-8.out > $6-9.out
awk '{print $0"\t"($7+$8+$9+$10)}' $6-9.out >  $6-10.out
awk '{print $0"\t"($11+$12+$13+$14)}' $6-10.out >  $6-11.out
awk '{print $0"\t"($15+$16+$17+$18)}' $6-11.out >  $6-12.out
awk '{print $0"\t"($19+$20+$21+$22)}' $6-12.out >  $6-13.out
awk '{print $0"\t"($26+$27)}' $6-13.out >  $6-14.out
awk '{print $0"\t"($28+$29)}' $6-14.out >  $6-15.out

awk '{if ($30 > 0) {print $0"\t"$26/$30} }' $6-15.out >  $6-16.out
awk '{if ($31 > 0) {print $0"\t"$28/$31} }' $6-16.out >  $6-17.out
awk '{$4=$7=$8=$9=$10=$11=$12=$13=$14=$15=$16=$17=$18=$19=$20=$21=$22=$24=""; print $0}' $6-17.out > $6-18.out
awk '{ $16="in";}1' OFS="\t" $6-18.out >  $6-19.out
awk '{ $17="III";}1' OFS="\t" $6-19.out >  $6-20.out
rm $6-1.out $6-3.out $6-4.out $6-5.out $6-6.out $6-7.out $6-8.out $6-9.out $6-10.out $6-11.out $6-12.out $6-13.out $6-14.out $6-15.out $6-16.out $6-17.out $6-18.out $6-19.out

awk '{ if($9 != "*") print $0;}1' OFS="\t" $3 > $7-1.out
awk '{ if($10 == "*") $10="0,0,0,0";}1' OFS="\t" $7-1.out > $7-3.out
awk '{ if($8 == "*") $8="0,0,0,0";}1' OFS="\t" $7-3.out > $7-4.out
awk '{ if($7 == "*") $7="0,0,0,0";}1' OFS="\t" $7-4.out > $7-5.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$10)}1' $7-5.out > $7-6.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$9)}1' $7-6.out > $7-7.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$8)}1' $7-7.out > $7-8.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$7)}1' $7-8.out > $7-9.out
awk '{print $0"\t"($7+$8+$9+$10)}' $7-9.out >  $7-10.out
awk '{print $0"\t"($11+$12+$13+$14)}' $7-10.out >  $7-11.out
awk '{print $0"\t"($15+$16+$17+$18)}' $7-11.out >  $7-12.out
awk '{print $0"\t"($19+$20+$21+$22)}' $7-12.out >  $7-13.out
awk '{print $0"\t"($26+$27)}' $7-13.out >  $7-14.out
awk '{print $0"\t"($28+$29)}' $7-14.out >  $7-15.out
awk '{if ($30 > 0) {print $0"\t"$26/$30} }' $7-15.out >  $7-16.out
awk '{if ($31 > 0) {print $0"\t"$28/$31} }' $7-16.out >  $7-17.out
awk '{$4=$7=$8=$9=$10=$11=$12=$13=$14=$15=$16=$17=$18=$19=$20=$21=$22=$24=""; print $0}' $7-17.out > $7-18.out

awk '{ $16="out";}1' $7-18.out >  $7-19.out
awk '{ $17="III";}1' $7-19.out >  $7-20.out
rm $7-1.out $7-3.out $7-4.out $7-5.out $7-6.out $7-7.out $7-8.out $7-9.out $7-10.out $7-11.out $7-12.out $7-13.out $7-14.out $7-15.out $7-16.out $7-17.out $7-18.out $7-19.out

awk -v OFS="\t" '$1=$1' $6-20.out > $4
awk -v OFS="\t" '$1=$1' $7-20.out > $5
cat $4 $5 > $8
