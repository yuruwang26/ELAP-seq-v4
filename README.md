# ELAP-seq-v4

Each replicate should contain one input library and one IP library. The libraries are built with superscript IV or III.

## 1. processing of Fastq reads

Only reads R2 is used. After trimming UMI, the begnning of R2 will be the RT stop site. Can otherwise design the library construction workflow in a way so that R1 from single end sequencing is sufficient for analysis.
pre-processing.sh is used. This includes five steps:

### 1) trim adapter

```bash
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-input-III-rep1-cutadapt.fq.gz HeLa-input-III-rep1_R2.fq.gz >> adaptorTrim.log
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-input-IV-rep1-cutadapt.fq.gz HeLa-input-IV-rep1_R2.fq.gz >> adaptorTrim.log
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-IP-III-rep1-cutadapt.fq.gz HeLa-IP-III-rep1_R2.fq.gz >> adaptorTrim.log
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-IP-IV-rep1-cutadapt.fq.gz HeLa-IP-IV-rep1_R2.fq.gz >> adaptorTrim.log
```

### 2) remove duplicates

```bash
~/Tools/bbmap/clumpify.sh in=HeLa-input-III-rep1-cutadapt.fq.gz out=HeLa-input-III-rep1-dedupe.fq.gz dedupe >> duplicates_removal_1.log
~/Tools/bbmap/clumpify.sh in=HeLa-input-IV-rep1-cutadapt.fq.gz out=HeLa-input-IV-rep1-dedupe.fq.gz dedupe >> duplicates_removal_1.log
~/Tools/bbmap/clumpify.sh in=HeLa-IP-III-rep1-cutadapt.fq.gz out=HeLa-IP-III-rep1-dedupe.fq.gz dedupe >> duplicates_removal_1.log
~/Tools/bbmap/clumpify.sh in=HeLa-IP-IV-rep1-cutadapt.fq.gz out=HeLa-IP-IV-rep1-dedupe.fq.gz dedupe >> duplicates_removal_1.log
```

### 3) trim UMI

```bash
conda activate cutadaptenv
```

```bash
cutadapt -u 6 -o tmp.trimmed.fastq HeLa-input-III-rep1-dedupe.fq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-input-III-rep1-trim.fq.gz tmp.trimmed.fastq
```

```bash
cutadapt -u 7 -o tmp.trimmed.fastq HeLa-input-IV-rep1-dedupe.fq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-input-IV-rep1-trim.fq.gz tmp.trimmed.fastq
```
```bash
cutadapt -u 6 -o tmp.trimmed.fastq HeLa-IP-III-rep1-dedupe.fq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-IP-III-rep1-trim.fq.gz tmp.trimmed.fastq
```

```bash
cutadapt -u 7 -o tmp.trimmed.fastq HeLa-IP-IV-rep1-dedupe.fq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-IP-IV-rep1-trim.fq.gz tmp.trimmed.fastq
```
```bash
conda deactivate
```
### 4) map reads to the genome

```bash
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file YW_align_summary -p 4 -U HeLa-input-III-rep1-trim.fq.gz |samtools view -bS |samtools sort -o HeLa-input-III-rep1.bam
samtools index HeLa-input-III-rep1.bam HeLa-input-III-rep1.bai
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file YW_align_summary -p 4 -U HeLa-input-IV-rep1-trim.fq.gz |samtools view -bS |samtools sort -o HeLa-input-IV-rep1.bam
samtools index HeLa-input-IV-rep1.bam HeLa-input-IV-rep1.bai
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file YW_align_summary -p 4 -U HeLa-IP-III-rep1-trim.fq.gz |samtools view -bS |samtools sort -o HeLa-IP-III-rep1.bam
samtools index HeLa-IP-III-rep1.bam HeLa-IP-III-rep1.bai
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file YW_align_summary -p 4 -U HeLa-IP-IV-rep1-trim.fq.gz |samtools view -bS |samtools sort -o HeLa-IP-IV-rep1.bam
samtools index HeLa-IP-IV-rep1.bam HeLa-IP-IV-rep1.bai
```

### 5) combine .bam files from III and IV
```bash
samtools merge HeLa-input-III-IV-rep1.bam HeLa-input-III-rep1.bam HeLa-input-IV-rep1.bam
samtools merge HeLa-IP-III-IV-rep1.bam HeLa-IP-III-rep1.bam HeLa-IP-IV-rep1.bam
samtools index HeLa-input-III-IV-rep1.bam HeLa-input-III-IV-rep1.bai
samtools index HeLa-IP-III-IV-rep1.bam HeLa-IP-III-IV-rep1.bai
```

## 2. Call IP peaks
This step uses the script call.sh

```bash
macs2 callpeak -t HeLa-IP-III-rep1.bam -c HeLa-input-III-rep1.bam -n test_t2 -f BAM -g 994080837 -q 0.01 --slocal 1000 --extsize 150 --nomodel --keep-dup all --call-summits --outdir HeLa-III-peadDir
macs2 callpeak -t HeLa-IP-III-IV-rep1.bam -c HeLa-input-III-IV-rep1.bam -n test_t2 -f BAM -g 994080837 -q 0.01 --slocal 1000 --extsize 150 --nomodel --keep-dup all --call-summits --outdir HeLa-III-IV-peadDir
```

## 3. Call arrested sites inside and outside of the IP peaks. 
This process includes three steps:
### 1) obtain regions covered by IP peaks 
Extend the IP peaks by -5 and +5 of the start and end positions respectively.
Remove unknown chromosomes or random chromosomes manually
Save the resulting .bed file as HeLa-peaks.bed

### 2) Call all stop sites inside and outside of IP peaks
This step uses the script Arrest.sh
```bash
bash Arrest.sh HeLa-input-III-rep1.bam HeLa-IP-III-rep1.bam HeLa-peaks.bed HeLa-III-rep1-inside.out HeLa-III-rep1-outside.out
bash Arrest.sh HeLa-input-III-IV-rep1.bam HeLa-IP-III-IV-rep1.bam HeLa-peaks.bed HeLa-III-IV-rep1-inside.out HeLa-III-IV-rep1-outside.out
```
### 3) Calculate arrest rate of each site and assign the originality of the site (whether it is from the sample built with superscript III or combined samples built with superscript III and IV, whether it is inside of the IP peak or outside of the IP peak)
This step uses scripts calculate1.sh and calculate2.sh
```bash
bash calculate1.sh HeLa-IP-III-rep1.bam HeLa-III-rep1-inside.out HeLa-III-rep1-outside.out HeLa-III-rep1-inside-unfiltered.bed HeLa-III-rep1-outside-unfiltered.bed HeLa-III-rep1-unfiltered.bed
bash calculate2.sh HeLa-IP-III-IV-rep1.bam HeLa-III-IV-rep1-inside.out HeLa-III-IV-rep1-outside.out HeLa-III-IV-rep1-inside-unfiltered.bed HeLa-III-IV-rep1-outside-unfiltered.bed HeLa-III-rep1-unfiltered.bed
```

## 4. Remove backgrounds 
This process removes background originated from multiple mapping and stutter effect 
### 1) select the most abundant nucleoside
```bash
python3 R1.py HeLa-III-rep1-inside-unfiltered.bed > HeLa-III-rep1-inside-unfiltered-ab.bed
python3 R1.py HeLa-III-rep1-outside-unfiltered.bed > HeLa-III-rep1-outside-unfiltered-ab.bed
cat HeLa-III-rep1-inside-unfiltered-ab.bed | tr ' ' '\t' > HeLa-III-rep1-inside-unfiltered-1.bed
cat HeLa-III-rep1-outside-unfiltered-ab.bed | tr ' ' '\t' > HeLa-III-rep1-outside-unfiltered-1.bed
```
### 2) remove sites that are covered by reads that all share the same start and end mapping position
```bash
python3 R2.py HeLa-III-rep1-inside-unfiltered-1.bed > HeLa-III-rep1-inside-block.bed
python3 R2.py HeLa-III-rep1-outside-unfiltered-1.bed > HeLa-III-rep1-outside-block.bed
cat HeLa-III-rep1-inside-block.bed | tr ' ' '\t' > HeLa-III-rep1-inside-block-1.bed
cat HeLa-III-rep1-outside-block.bed | tr ' ' '\t' > HeLa-III-rep1-outside-block-1.bed
bedtools subtract -a HeLa-III-rep1-inside-unfiltered-1.bed -b HeLa-III-rep1-inside-block-1.bed > HeLa-III-rep1-inside-unfiltered-2.bed
bedtools subtract -a HeLa-III-rep1-outside-unfiltered-1.bed -b HeLa-III-rep1-outside-block-1.bed > HeLa-III-rep1-outside-unfiltered-2.bed
```

### 3) Filter sites due to peak tails (sites within 30 nt has > 5 fold coverage than this site)
```bash
python3 R3.py HeLa-III-rep1-inside-unfiltered-2.bed > HeLa-III-rep1-inside-unfiltered-low.bed
python3 R3.py HeLa-III-rep1-outside-unfiltered-2.bed > HeLa-III-rep1-outside-unfiltered-low.bed
awk '!visited[$0]++' HeLa-III-rep1-inside-unfiltered-low.bed > HeLa-III-rep1-inside-unfiltered-low-1.bed
awk '!visited[$0]++' HeLa-III-rep1-outside-unfiltered-low.bed > HeLa-III-rep1-outside-unfiltered-low-1.bed
cat HeLa-III-rep1-inside-unfiltered-low-1.bed | tr ' ' '\t' > HeLa-III-rep1-inside-unfiltered-low-2.bed
cat HeLa-III-rep1-outside-unfiltered-low-1.bed | tr ' ' '\t' > HeLa-III-rep1-outside-unfiltered-low-2.bed
bedtools subtract -a HeLa-III-rep1-inside-unfiltered-2.bed -b HeLa-III-rep1-inside-unfiltered-low-2.bed > HeLa-III-rep1-inside-unfiltered-3.bed
bedtools subtract -a HeLa-III-rep1-outside-unfiltered-2.bed -b HeLa-III-rep1-outside-unfiltered-low-2.bed > HeLa-III-rep1-outside-unfiltered-3.bed
```

## 5. filter sites based on stop ratio and coverage

```bash
bash filter1.sh HeLa-IP-III-rep1.bam HeLa-III-rep1-inside-unfiltered-3.bed HeLa-III-rep1-outside-unfiltered-3.bed HeLa-III-rep1-unique.bed
bash filter2.sh HeLa-IP-III-IV-rep1.bam HeLa-III-IV-rep1-inside-unfiltered-3.bed HeLa-III-IV-rep1-outside-unfiltered-3.bed HeLa-III-IV-rep1-unique.bed
sort -k1,1 -k2,2n HeLa-III-rep1-unique.bed > HeLa-III-rep1-unique-1.bed 
sort -k1,1 -k2,2n HeLa-III-IV-rep1-unique.bed > HeLa-III-IV-rep1-unique-1.bed
```
## 6. remove stutter sites
```bash
python3 stutter1_III.py HeLa-III-rep1-unique-1.bed > HeLa-III-rep1-stutter-filter.bed
python3 stutter1_III_IV.py HeLa-III-IV-rep1-unique-1.bed > HeLa-III-IV-rep1-stutter-filter.bed
cat HeLa-III-rep1-stutter-filter.bed | tr ' ' '\t' > HeLa-III-rep1-stutter-filter-1.bed
cat HeLa-III-IV-rep1-stutter-filter.bed | tr ' ' '\t' > HeLa-III-IV-rep1-stutter-filter-1.bed
python3 stutter2_III.py HeLa-III-rep1-stutter-filter-1.bed > HeLa-III-rep1-remove.bed
python3 stutter2_III_IV.py HeLa-III-IV-rep1-stutter-filter-1.bed > HeLa-III-IV-rep1-remove.bed
cat HeLa-III-rep1-remove.bed | tr ' ' '\t' > HeLa-III-rep1-remove-1.bed
cat HeLa-III-IV-rep1-remove.bed | tr ' ' '\t' > HeLa-III-IV-rep1-remove-1.bed
bedtools subtract -a HeLa-III-rep1-stutter-filter-1.bed -b HeLa-III-rep1-remove-1.bed > HeLa-III-rep1-stutter-filter-2.bed
bedtools subtract -a HeLa-III-IV-rep1-stutter-filter-1.bed -b HeLa-III-IV-rep1-remove-1.bed > HeLa-III-IV-rep1-stutter-filter-2.bed
```
## 7. merge all sites identified from III and III+IV and replicate 1 and replicate 2.
### 1) 
```bash
bedtools subtract -a HeLa-III-IV-rep1-stutter-filter-2.bed -b HeLa-III-rep1-stutter-filter-2.bed > new.bed
cat HeLa-III-rep1-stutter-filter-2.bed new.bed > HeLa-rep1-combined.bed
sort -k1,1 HeLa-rep1-combined.bed > HeLa-rep1-combined-1.bed
```

### 2) overlap with the III-IV-unfiltered.bed to obtain reads in the input file and IP files of III+IV which are used for quantification later
need to maually change the unfiltered file and combined-1 file in order to select for the same base as well. save as unfiltered-mn.bed and combined-1-1.bed respectively.
```bash
#intersect with input file
bedtools intersect -wa -wb -a ../../$s-III-IV-unfiltered-ab.bed -b HeLa-rep1-combined-1.bed > HeLa-rep1-combined-2.bed
awk '{print $1,$2,$3,$5,$7,$12,$13,$14,$28,$29,$30,$31,$32,$33}' HeLa-rep1-combined-2.bed > HeLa-rep1-combined-3.bed
awk -v OFS="\t" '$1=$1' HeLa-rep1-combined-3.bed > HeLa-rep1-combined-4.bed
```

### 3) Intersect two replicates
```bash
bedtools intersect -wa -wb -a HeLa-rep1-combined-4.bed -b HeLa-rep2-combined-4.bed > HeLa.bed
awk '!visited[$0]++' HeLa.bed > HeLa-1.bed
awk '{print $1,$2,$3,$4,$5,$8,$6,$7,$11,$12,$9,$10,$13,$14,$22,$20,$21,$25,$26,$23,$24,$27,$28}' HeLa-1.bed > HeLa-2.bed
awk -v OFS="\t" '$1=$1' HeLa-2.bed > HeLa-3.bed
cat HeLa-3.bed | tr ' ' '\t' > HeLa-4.bed
sort -k1,1 -k2,2n HeLa-4.bed > HeLa-sort.bed
awk '{print $0"\t"($9+$18)/2}' HeLa-sort.bed > HeLa-input-avg.bed
awk '{print $0"\t"($10+$19)/2}' HeLa-input-avg.bed > HeLa-IP-avg.bed
```

## 8. Final filter
### 1) remove stutter site again
```bash
python3 stutter-final.py HeLa-IP-avg.bed > HeLa-stop-filter.bed
```
### 2) select sites fulfiling cutoffs for average stop ratios

## 9. post-processing

then manually split the chr and base to the orignal columns

```bash
awk '{print $1,$2,$3,$5,$7,$12,$13,$14,$26,$27,$28}' ./K1/K1-combined-2.bed > ./K1/K1-combined-3.bed
awk -v OFS="\t" '$1=$1' ./K1/K1-combined-3.bed > ./K1/K1-combined-4.bed
```
need to do another stutter filter here

#### Intersect two replicates
```bash
bedtools intersect -wa -wb -a ./K1/K1-combined-4.bed -b ./K3/K3-combined-4.bed > K1-K3-v13.bed
awk '{print $1,$2,$3,$4,$5,$8,$6,$7,$11,$9,$10,$19,$17,$18,$22,$20,$21}' K1-K3-v13.bed > K1-K3-v13-1.bed
awk -v OFS="\t" '$1=$1' K1-K3-v13-1.bed > K1-K3-v13-2.bed
sort -k1,1 -k2,2n L2-L3-v13-2.bed > L2-L3-v13-3.bed
cat HeLa.csv | tr ' ' '\t' > HeLa.bed
```
#### filter again the stutter sites
```bash
python3 stop_v14_second_filter.py K1-K3-v13-2.bed > HEK-v14.csv
cat HeLa.csv | tr ' ' '\t' > HeLa.bed
awk '{OFS=" "; print $1,$2,$3,$4,$5,$6,($7/3.289),($8/3.55),$9,$10,$11,$12,($13/3.133),($14/3.87),$15,$16,$17}' OFS="\t" HEK-v14.bed > HEK-v14-RPM.bed
```

combine1.sh includes three steps:

### 1) prepare chromosome locations inside and outside the IP peaks

Extend the IP peaks (.xls file) by -5 and +5 of the start and end positions respectively.
Remove unknown chromosomes and random chromosomes manually
Save the resulting .bed file as L2-peaks-end.bed
merge chromosome locations:

```bash
bedtools merge -i L2-peaks-end.bed > 6.bed
```

Obtain all regions covered by reads in the IP sample:

```bash
bedtools bamtobed -split -i L2-III-mRNA-IP.bam > 1.bed
awk '{print $1,$2,$3}' 1.bed OFS="\t" > 2.bed
cat 2.bed | tr ' ' '\t' > 3.bed
sort -k1,1 -k2,2n 3.bed > 4.bed
bedtools merge -i 4.bed > 5.bed
```

Get regions outside of IP peaks:

```bash
bedtools subtract -a 5.bed -b 6.bed > 7.bed
```


### 2) call RT-arrest ratio for sites inside and outside IP peaks using JACUSA

```bash
java -jar /home/Wang_yuru/Pseudouridine/4-21-22-HeLa-HEK-mRNA/cutadapt/dedupe/trim/sam_files/JACUSA_v2.0.1.jar rt-arrest -m 0 -p 2 -c 1 -P FR_SECONDSTRAND -R /home/Wang_yuru/Database/genome/hg38/hg38_UCSC.fa -b 6.bed -r 6.out $1 $2
java -jar /home/Wang_yuru/Pseudouridine/4-21-22-HeLa-HEK-mRNA/cutadapt/dedupe/trim/sam_files/JACUSA_v2.0.1.jar rt-arrest -m 0 -p 2 -c 1 -P FR_SECONDSTRAND -R /home/Wang_yuru/Database/genome/hg38/hg38_UCSC.fa -b 7.bed -r 7.out $1 $2
```

combine2-p1.sh includes three steps:

### 3) filter sites

sites inside IP peaks:
select for sites whose arrest in the IP sample is not 0;
select for sites whose identity is T in sequencing;
assign the value to 0,0,0,0 if reat-through in IP sample is \*;
calculate stop ratio and filter


```bash
awk '{ if($9 != "*") print $0;}' 6.out > 6-1.out
awk '{ if($13 == "T") print $0;}' 6-1.out > 6-2.out
awk '{ if($10 == "*") $10="0,0,0,0";}1' OFS="\t" 6-2.out > 6-3.out
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
awk '{print $0"\t"($23+$24)}' 6-13.out >  6-14.out
awk '{print $0"\t"($25+$26)}' 6-14.out >  6-15.out
awk '{print $0"\t"$25/$28}' 6-15.out >  6-16.out
awk '{$7=$8=$9=$10=$11=$12=$13=$14=$15=$16=$17=$18=$19=$20=$21=$22=""; print $0}' 6-17.out > 6-18.out
rm 6-1.out 6-2.out 6-3.out 6-4.out 6-5.out 6-6.out 6-7.out 6-8.out 6-9.out 6-10.out 6-11.out 6-12.out 6-13.out 6-14.out 6-15.out 6-16.out 6-17.out 
```

sites outside IP peaks:

```bash
awk '{ if($9 != "*") print $0;}' 7.out > 7-1.out
awk '{ if($13 == "T") print $0;}' 7-1.out > 7-2.out
awk '{ if($10 == "*") $10="0,0,0,0";}1' OFS="\t" 7-2.out > 7-3.out
awk '{ if($8 == "*") $8="0,0,0,0";}1' OFS="\t" 7-3.out > 7-4.out
awk '{ if($7 == "*") $7="0,0,0,0";}1' OFS="\t" 7-4.out > 7-5.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$10)}1' 7-5.out > 7-6.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$9)}1' 7-6.out > 7-7.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$8)}1' 7-7.out > 7-8.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$7)}1' 7-8.out > 7-9.out
awk '{print $0"\t"($7+$8+$9+$10)}' 7-9.out >  7-10.out
awk '{print $0"\t"($11+$12+$13+$14)}' 7-10.out >  7-11.out
awk '{print $0"\t"($15+$16+$17+$18)}' 7-11.out >  7-12.out
awk '{print $0"\t"($19+$20+$21+$22)}' 7-12.out >  7-13.out
awk '{print $0"\t"($23+$24)}' 7-13.out >  7-14.out
awk '{print $0"\t"($25+$26)}' 7-14.out >  7-15.out
awk '{print $0"\t"$25/$28}' 7-15.out >  7-16.out
awk '{$7=$8=$9=$10=$11=$12=$13=$14=$15=$16=$17=$18=$19=$20=$21=$22=""; print $0}' 7-17.out > 7-18.out
rm 7-1.out 7-2.out 7-3.out 7-4.out 7-5.out 7-6.out 7-7.out 7-8.out 7-9.out 7-10.out 7-11.out 7-12.out 7-13.out 7-14.out 7-15.out 7-16.out 7-17.out
```
combine2-p2.sh
```bash 
awk '{ if(($29 > 0.2) || ($25 > 15 && $29 > 0.10 )) print $0;}' 6-16.out > 6-17.out
awk '{ if(($29 > 0.4 && $28 > 10) || ($25 > 10 && $29 > 0.2 )) print $0;}' 7-16.out > 7-17.out
```
### 4. Intersect between pU sites and locations that are covered by uniquely mapped reads in the IP sample

need to change the format of the files to tab delimited

```bash
awk -v OFS="\t" '$1=$1' 6-18.out > 6-18.bed
awk -v OFS="\t" '$1=$1' 7-18.out > 7-18.bed
```

merge sites from inside and outside IP peaks:
```bash
cat 6-18.bed 7-18.bed > all.bed
```

obtain all locations that are covered by uniquely mapped reads

```bash
bedtools bamtobed -split -i $2 > coverage.bed
awk '{ if($5 == 60) print $0;}' coverage.bed > MQ-60.bed
awk '{print $1,$2,$3}' MQ-60.bed OFS="\t" > MQ-60-1.bed
cat MQ-60-1.bed | tr ' ' '\t' > MQ-60-2.bed
sort -k1,1 -k2,2n MQ-60-2.bed > MQ-60-3.bed
bedtools merge -i MQ-60-3.bed > MQ-60-final.bed
rm coverage.bed MQ-60.bed MQ-60-1.bed MQ-60-2.bed MQ-60-3.bed
```

Intersect pU sites list with MQ-60-final.bed

```bash
bedtools intersect -a all.bed -b MQ-60-final.bed > $4
```

## 6. Overlap sites from two replicates
```bash
bedtools intersect -a L2-unique.bed -b L3-unique.bed > L2-L3.bed
```

## 7. remove sites due to stutter
```bash
awk '{print $1,$2,$3,$6,$10,$11,$12,$13,$14,$16,$15}' L2-L3.bed > L2-L3-1.bed
awk -v OFS="\t" '$1=$1' L2-L3-1.bed > L2-L3-2.bed
sort -k1,1 -k2,2n L2-L3-v13-2.bed > L2-L3-v13-3.bed
python3 L_stop.py > HeLa.csv
cat HeLa.csv | tr ' ' '\t' > HeLa.bed
```

## 8. Notes

Combine.sh

```bash
#!/bin/bash

###
### combine.sh
###    - get coverage regions, run JACUSA under the regions, and filter sites
###
###
### Usage:
###   combine.sh <input1.bam> <input2.bam> <IP-peaks.bed> <filter.bed>
###
### Options:
###   <input1.bam>         Input bam file to read.
###   <input2.bam>         IP bam file to read
###   <IP-peaks.bed>      IP peaks corresonding to the bam
###   <inside-filter.out>            sites in uniquely mapped regions
###
###
### Output:
###  chr start end type p-value strand arrest_A arrest_C arrest_G arrest_T readthrough_A readthrough_C readthrough_G readthrough_T read arrest_ratio
###
### Contact:
###  Yuru Wang (ywang26@uchicago.edu)
###

bedtools bamtobed -split -i $2 > 1.bed
awk '{print $1,$2,$3}' 1.bed OFS="\t" > 2.bed
cat 2.bed | tr ' ' '\t' > 3.bed
sort -k1,1 -k2,2n 3.bed > 4.bed
bedtools merge -i 4.bed > 5.bed

bedtools merge -i $3 > 6.bed

bedtools subtract -a 5.bed -b 6.bed > 7.bed

rm 1.bed 2.bed 3.bed 4.bed 5.bed
 

java -jar /home/Wang_yuru/Pseudouridine/4-21-22-HeLa-HEK-mRNA/cutadapt/dedupe/trim/sam_files/JACUSA_v2.0.1.jar rt-arrest -m 0 -p 2 -P FR_SECONDSTRAND -R /home/Wang_yuru/Database/genome/hg38/hg38_UCSC.fa -b 6.bed -r 6.out $1 $2
java -jar /home/Wang_yuru/Pseudouridine/4-21-22-HeLa-HEK-mRNA/cutadapt/dedupe/trim/sam_files/JACUSA_v2.0.1.jar rt-arrest -m 0 -p 2 -P FR_SECONDSTRAND -R /home/Wang_yuru/Database/genome/hg38/hg38_UCSC.fa -b 7.bed -r 7.out $1 $2

awk '{ if($9 != "*") print $0;}' 6.out > 6-1.out
awk '{ if($13 == "T") print $0;}' 6-1.out > 6-2.out
awk '{ if($10 == "*") $10="0,0,0,0";}1' OFS="\t" 6-2.out > 6-3.out
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

awk '{print $0"\t"$28/$31}' 6-15.out >  6-16.out
awk '{ if(($32 > 0.2) || ($28 > 15 && $32 > 0.10 )) print $0;}' 6-16.out > 6-17.out
awk '{$7=$8=$9=$10=$11=$12=$13=$14=$15=$16=$17=$18=$19=$20=$21=$22=""; print $0}' 6-17.out > 6-18.out
rm 6-1.out 6-2.out 6-3.out 6-4.out 6-5.out 6-6.out 6-7.out 6-8.out 6-9.out 6-10.out 6-11.out 6-12.out 6-13.out 6-14.out 6-15.out 6-16.out 6-17.out

awk '{ if($9 != "*") print $0;}' 7.out > 7-1.out
awk '{ if($13 == "T") print $0;}' 7-1.out > 7-2.out
awk '{ if($10 == "*") $10="0,0,0,0";}1' OFS="\t" 7-2.out > 7-3.out
awk '{ if($8 == "*") $8="0,0,0,0";}1' OFS="\t" 7-3.out > 7-4.out
awk '{ if($7 == "*") $7="0,0,0,0";}1' OFS="\t" 7-4.out > 7-5.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$10)}1' 7-5.out > 7-6.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$9)}1' 7-6.out > 7-7.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$8)}1' 7-7.out > 7-8.out
awk 'BEGIN{FS=OFS="\t"}{gsub(",","\t",$7)}1' 7-8.out > 7-9.out
awk '{print $0"\t"($7+$8+$9+$10)}' 7-9.out >  7-10.out
awk '{print $0"\t"($11+$12+$13+$14)}' 7-10.out >  7-11.out
awk '{print $0"\t"($15+$16+$17+$18)}' 7-11.out >  7-12.out
awk '{print $0"\t"($19+$20+$21+$22)}' 7-12.out >  7-13.out
awk '{print $0"\t"($26+$27)}' 7-13.out >  7-14.out
awk '{print $0"\t"($28+$29)}' 7-14.out >  7-15.out

awk '{print $0"\t"$28/$31}' 7-15.out >  7-16.out
awk '{ if(($32 > 0.4 && $31 > 10) || ($28 > 10 && $32 > 0.2 )) print $0;}' 7-16.out > 7-17.out
awk '{$7=$8=$9=$10=$11=$12=$13=$14=$15=$16=$17=$18=$19=$20=$21=$22=""; print $0}' 7-17.out > 7-18.out
rm 7-1.out 7-2.out 7-3.out 7-4.out 7-5.out 7-6.out 7-7.out 7-8.out 7-9.out 7-10.out 7-11.out 7-12.out 7-13.out 7-14.out 7-15.out 7-16.out 7-17.out 

awk -v OFS="\t" '$1=$1' 6-18.out > 6-18.bed
awk -v OFS="\t" '$1=$1' 7-18.out > 7-18.bed
cat 6-18.bed 7-18.bed > all.bed

bedtools bamtobed -split -i $2 > coverage.bed
awk '{ if($5 == 60) print $0;}' coverage.bed > MQ-60.bed
awk '{print $1,$2,$3}' MQ-60.bed OFS="\t" > MQ-60-1.bed
cat MQ-60-1.bed | tr ' ' '\t' > MQ-60-2.bed
sort -k1,1 -k2,2n MQ-60-2.bed > MQ-60-3.bed
bedtools merge -i MQ-60-3.bed > MQ-60-final.bed
rm coverage.bed MQ-60.bed MQ-60-1.bed MQ-60-2.bed MQ-60-3.bed
bedtools intersect -a all.bed -b MQ-60-final.bed > $4
```

L2-III-stop.py
```bash
import pandas as pd
import csv

reader = pd.read_csv('L2-L3-2.bed', delimiter="\t", names=["chr","p", "position","strand", "in_arrest", "in_readthrough", "IP_arrest","IP_readthrough","in_sum","IP_stop","IP_sum"], chunksize=1000000)
#reader.ncolumns = ["chr","p", "position","strand", "in_arrest", "in_readthrough", "IP_arrest","IP_readthrough","in_sum","IP_stop","IP_sum"]
for df in reader:
    #print (r + len(df))
    df["position"] = pd.to_numeric(df["position"])
    for i in range(len(df)):
        pos1 = df["position"].iloc[i] - 1
        pos2 = df["position"].iloc[i] + 1
        if df["strand"].iloc[i] == "-":
            if i-1 > 0 and df["position"].iloc[i-1] == pos1:
                stop_decrease = (df["IP_stop"].iloc[i]-df["IP_stop"].iloc[i-1])/df["IP_stop"].iloc[i]
                read_increase = (df["IP_sum"].iloc[i-1]-df["IP_sum"].iloc[i])/df["IP_sum"].iloc[i-1]
                if i+1 < len(df) and df["position"].iloc[i+1] == pos2:
                    read_decrease = (df["IP_sum"].iloc[i]-df["IP_sum"].iloc[i+1])/df["IP_sum"].iloc[i]
                    if stop_decrease > 0.2 and read_increase < 0.2 and read_decrease > 0.2:
                        print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["in_arrest"].iloc[i],df["in_readthrough"].iloc[i],df["IP_arrest"].iloc[i],df["IP_readthrough"].iloc[i],df["in_sum"].iloc[i],df["IP_sum"].iloc[i],df["IP_stop"].iloc[i])
                elif i+1 < len(df) and df["position"].iloc[i+1] != pos2:
                    if stop_decrease > 0.2 and read_increase < 0.2:
                        print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["in_arrest"].iloc[i],df["in_readthrough"].iloc[i],df["IP_arrest"].iloc[i],df["IP_readthrough"].iloc[i],df["in_sum"].iloc[i],df["IP_sum"].iloc[i],df["IP_stop"].iloc[i])
                
            elif i-1 > 0 and df["position"].iloc[i-1] != pos1:     
                if i+1 < len(df) and df["position"].iloc[i+1] == pos2:
                    read_decrease = (df["IP_sum"].iloc[i]-df["IP_sum"].iloc[i+1])/df["IP_sum"].iloc[i]
                    if read_decrease > 0.2:
                        print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["in_arrest"].iloc[i],df["in_readthrough"].iloc[i],df["IP_arrest"].iloc[i],df["IP_readthrough"].iloc[i],df["in_sum"].iloc[i],df["IP_sum"].iloc[i],df["IP_stop"].iloc[i])
                elif i+1 < len(df) and df["position"].iloc[i+1] != pos2:
                    print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["in_arrest"].iloc[i],df["in_readthrough"].iloc[i],df["IP_arrest"].iloc[i],df["IP_readthrough"].iloc[i],df["in_sum"].iloc[i],df["IP_sum"].iloc[i],df["IP_stop"].iloc[i])
        elif df["strand"].iloc[i] == "+":
            if i-1 > 0 and df["position"].iloc[i-1] == pos1:
                read_decrease = (df["IP_sum"].iloc[i]-df["IP_sum"].iloc[i-1])/df["IP_sum"].iloc[i]
                if i+1 < len(df) and df["position"].iloc[i+1] == pos2:
                    stop_decrease = (df["IP_stop"].iloc[i]-df["IP_stop"].iloc[i+1])/df["IP_stop"].iloc[i]
                    read_increase = (df["IP_sum"].iloc[i+1]-df["IP_sum"].iloc[i])/df["IP_sum"].iloc[i+1]
                    if read_decrease > 0.2 and read_increase < 0.2 and stop_decrease > 0.2:
                        print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["in_arrest"].iloc[i],df["in_readthrough"].iloc[i],df["IP_arrest"].iloc[i],df["IP_readthrough"].iloc[i],df["in_sum"].iloc[i],df["IP_sum"].iloc[i],df["IP_stop"].iloc[i])
                elif i+1 < len(df) and df["position"].iloc[i+1] != pos2:
                    if read_decrease > 0.2:
                        print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["in_arrest"].iloc[i],df["in_readthrough"].iloc[i],df["IP_arrest"].iloc[i],df["IP_readthrough"].iloc[i],df["in_sum"].iloc[i],df["IP_sum"].iloc[i],df["IP_stop"].iloc[i])
                
            elif i-1 > 0 and df["position"].iloc[i-1] != pos1:     
                if i+1 < len(df) and df["position"].iloc[i+1] == pos2:
                    stop_decrease = (df["IP_stop"].iloc[i]-df["IP_stop"].iloc[i+1])/df["IP_stop"].iloc[i]
                    read_increase = (df["IP_sum"].iloc[i+1]-df["IP_sum"].iloc[i])/df["IP_sum"].iloc[i+1]
                    if stop_decrease > 0.2 and read_increase < 0.2:
                        print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["in_arrest"].iloc[i],df["in_readthrough"].iloc[i],df["IP_arrest"].iloc[i],df["IP_readthrough"].iloc[i],df["in_sum"].iloc[i],df["IP_sum"].iloc[i],df["IP_stop"].iloc[i])
                elif i+1 < len(df) and df["position"].iloc[i+1] != pos2:
                    print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["in_arrest"].iloc[i],df["in_readthrough"].iloc[i],df["IP_arrest"].iloc[i],df["IP_readthrough"].iloc[i],df["in_sum"].iloc[i],df["IP_sum"].iloc[i],df["IP_stop"].iloc[i])
```

