#!/bin/bash

###
### arrest.sh 
###    - get coverage regions and run JACUSA under the regions to determine RT arrest
### 
###  
### Usage:
###   combine1.sh <input.bam> <IP.bam> <IP-peaks.bed> <inside.out> <outside.out>
###
### Options:
###   <input.bam>         Input bam file to read. 
###   <IP.bam>         IP bam file to read
###   <IP-peaks.bed>      IP peaks corresonding to the bam
###   <filter.out>           sites filtered

###
### Output:
###  chr start   end     name    pvalue  strand  arrest_input        through_input       arrest_IP        through_IP       arrest_score    filter  ref
###
### Contact:
###  Yuru Wang (yuruwang26@hotmail.com)
###


bedtools merge -i $3 > 6.bed

java -jar /home/yuruwang/Pseudouridine/4-21-22-HeLa-HEK-mRNA/cutadapt/dedupe/trim/sam_files/JACUSA_v2.0.1.jar rt-arrest -m 0 -p 2 -c 1 -P FR_SECONDSTRAND -R /home/yuruwang/Database/genome/hg38/hg38_UCSC.fa -b 6.bed -r $4 $1 $2