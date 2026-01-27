#!/bin/bash

###
### arrest.sh 
###    - get coverage regions and run JACUSA under the regions to determine RT arrest
### 
###  
### Usage:
###   arrest.sh <input.bam> <IP.bam> <IP-peaks.bed> <inside.out> <outside.out>
###
### Options:
###   <input.bam>         Input bam file to read. 
###   <IP.bam>         IP bam file to read
###   <IP-peaks.bed>      IP peaks corresonding to the bam
###   <inside-filter.out>            sites inside the IP peaks and filtered
###   <outside-filter.out>           sites outside the IP peaks and filtered
###
### Output:
###  chr start   end     name    pvalue  strand  arrest_input        through_input       arrest_IP        through_IP       arrest_score    filter  ref
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###

bedtools bamtobed -split -i $2 > 1.bed
awk '{print $1,$2,$3}' 1.bed OFS="\t" > 2.bed
cat 2.bed | tr ' ' '\t' > 3.bed
sort -k1,1 -k2,2n 3.bed > 4.bed
bedtools merge -i 4.bed > 5.bed

bedtools merge -i $3 > 6.bed

bedtools subtract -a 5.bed -b 6.bed > 7.bed

rm 1.bed 2.bed 3.bed 4.bed 5.bed
 

java -jar JACUSA_v2.0.1.jar rt-arrest -m 0 -p 2 -c 1 -P FR_SECONDSTRAND -R /home/yuruwang/Database/genome/hg38/hg38_UCSC.fa -b 6.bed -r $4 $1 $2
java -jar JACUSA_v2.0.1.jar rt-arrest -m 0 -p 2 -c 1 -P FR_SECONDSTRAND -R /home/yuruwang/Database/genome/hg38/hg38_UCSC.fa -b 7.bed -r $5 $1 $2


