# ELAP-seq-rRNA-workflow

Each replicate should contain an input library and an IP library. The input library and the IP library are built with superscript III.

## 1. processing of Fastq reads

Only R2 is used. After trimming UMI, the begnning of R2 will be the RT stop site. Can otherwise redesign the ligation adapters to use R1 from single end sequencing.

### 1) trim adapter

```bash
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-IP-III-rep1-cutadapt.fastq.gz HeLa-IP-III-rep1_R2.fq.gz >> adaptorTrim.log
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-input-III-rep1-cutadapt.fastq.gz HeLa-input-III-rep1_R2.fq.gz >> adaptorTrim.log
```

### 2) remove duplicates

```bash
~/Tools/bbmap/clumpify.sh in=HeLa-IP-III-rep1-cutadapt.fastq.gz out=HeLa-IP-III-rep1-dedupe.fastq.gz dedupe >> duplicates_removal_1.log
~/Tools/bbmap/clumpify.sh in=HeLa-input-III-rep1-cutadapt.fastq.gz out=HeLa-input-III-rep1-dedupe.fastq.gz dedupe >> duplicates_removal_1.log
```

### 3) trim UMI

```bash
conda activate cutadaptenv
```

```bash
cutadapt -u 6 -o tmp.trimmed.fastq HeLa-IP-III-rep1-dedupe.fastq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-IP-III-rep1-trim.fastq.gz tmp.trimmed.fastq
```

```bash
cutadapt -u 6 -o tmp.trimmed.fastq HeLa-input-III-rep1-dedupe.fastq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-input-III-rep1-trim.fastq.gz tmp.trimmed.fastq
```

## 2. map reads to the genome

```bash
hisat2 -x /home/Wang_yuru/rRNA_genome/Human_rRNA_1 -k 1 --no-softclip --mp 1000,999 --summary-file align_summary -p 4 HeLa-IP-III-rep1-trim.fastq.gz |samtools view -bS |samtools sort -o HeLa-IP-III-rRNA-rep1.bam
hisat2 -x /home/Wang_yuru/rRNA_genome/Human_rRNA_1 -k 1 --no-softclip --mp 1000,999 --summary-file align_summary -p 4 HeLa-input-III-rep1-trim.fastq.gz |samtools view -bS |samtools sort -o HeLa-input-III-rRNA-rep1.bam
```

## 3. call RT-arrest ratio for all sites 

```bash
java -jar /home/Wang_yuru/Pseudouridine/4-21-22-HeLa-HEK-mRNA/cutadapt/dedupe/trim/sam_files/JACUSA_v2.0.1.jar rt-arrest -m 0 -p 2 -c 1 -P FR_SECONDSTRAND -R /home/Wang_yuru/rRNA_genome/Human_rRNA_1.fa -b rRNA.bed -r rRNA-rep1.out HeLa-input-III-rRNA-rep1.bam HeLa-IP-III-rRNA-rep1.bam
```

## 4. filter sites

```bash
bash rRNA_filter_p1.sh rRNA-rep1.out rRNA-rep1.bed
```

```bash
python3 Stutter_removal_rRNA.py rRNA-rep1.bed > rRNA-rep1-stutter.bed
```
```bash
bash rRNA_filter_p2.sh rRNA-rep1-stutter.bed rRNA-rep1-filter.bed
```
```bash
bedtools intersect -wa -wb -a rRNA-rep1-filter.bed rRNA-rep2-filter.bed > rRNA-out.bed
```


# ELAP-seq-poly-A-RNA-workflow

Each replicate should contain two input libraries and an two IP libraries. The libraries are built with superscript IV or III.

## 1. processing of Fastq reads

Only reads R2 is used. After trimming UMI, the begnning of R2 is the RT stop site. Can otherwise design the library construction workflow in a way that single end sequencing is sufficient.
This process includes five steps:

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
### 1) MACS2 is used to call IP peaks.

```bash
macs2 callpeak -t HeLa-IP-III-rep1.bam -c HeLa-input-III-rep1.bam -n test_t2 -f BAM -g 994080837 -q 0.01 --slocal 1000 --extsize 150 --nomodel --keep-dup all --call-summits --outdir HeLa-III-peadDir
macs2 callpeak -t HeLa-IP-III-IV-rep1.bam -c HeLa-input-III-IV-rep1.bam -n test_t2 -f BAM -g 994080837 -q 0.01 --slocal 1000 --extsize 150 --nomodel --keep-dup all --call-summits --outdir HeLa-III-IV-peadDir
```
### 2) obtain regions covered by IP peaks 
Extend the IP peaks by -5 and +5 of the start and end positions respectively.
Remove unknown chromosomes or random chromosomes manually
Save the resulting .bed file as HeLa-peaks.bed

## 3. the downstream analysis uses a single command : ELAP-seq.sh
This command include the following steps:
### 1. Call arrested sites inside and outside of the IP peaks. 
This process includes two steps:

#### 1) Call all stop sites inside and outside of IP peaks
This step uses the script arrest.sh
```bash
bash arrest.sh HeLa-input-III-rep1.bam HeLa-IP-III-rep1.bam HeLa-peaks.bed HeLa-III-rep1-inside.out HeLa-III-rep1-outside.out
bash arrest.sh HeLa-input-III-IV-rep1.bam HeLa-IP-III-IV-rep1.bam HeLa-peaks.bed HeLa-III-IV-rep1-inside.out HeLa-III-IV-rep1-outside.out
```
#### 2) Calculate arrest rate of each site and assign the originality of the site (i.e., whether it is from the library built with superscript III or combined libraries built with superscript III and IV, and whether it is inside or outside of the IP peak)
This step uses scripts calculate1.sh for the III librairy and calculate2.sh for the III+IV library
```bash
bash calculate1.sh HeLa-IP-III-rep1.bam HeLa-III-rep1-inside.out HeLa-III-rep1-outside.out HeLa-III-rep1-inside-unfiltered.bed HeLa-III-rep1-outside-unfiltered.bed HeLa-III-in HeLa-III-out HeLa-III-rep1-unfiltered.bed
bash calculate2.sh HeLa-IP-III-IV-rep1.bam HeLa-III-IV-rep1-inside.out HeLa-III-IV-rep1-outside.out HeLa-III-IV-rep1-inside-unfiltered.bed HeLa-III-IV-rep1-outside-unfiltered.bed HeLa-III-IV-in HeLa-III-IV-out HeLa-III-rep1-unfiltered.bed
```

### 2. Remove backgrounds 
This process removes background originated from multiple mapping and low processivity of the RT enzyme
#### 1) determine the main nucleoside for sites mapped with multiple identities
```bash
python3 R1.py HeLa-III-rep1-inside-unfiltered.bed | tr ' ' '\t' > HeLa-III-rep1-inside-unfiltered-ab.bed
python3 R1.py HeLa-III-rep1-outside-unfiltered.bed | tr ' ' '\t' > HeLa-III-rep1-outside-unfiltered-ab.bed
```
#### 2) remove regions that are covered by reads that all share the same start and end mapping position
```bash
python3 R2.py HeLa-III-rep1-inside-unfiltered-ab.bed | tr ' ' '\t' > HeLa-III-rep1-inside-block.bed
python3 R2.py HeLa-III-rep1-outside-unfiltered-ab.bed | tr ' ' '\t' > HeLa-III-rep1-outside-block.bed
bedtools subtract -a HeLa-III-rep1-inside-unfiltered-ab.bed -b HeLa-III-rep1-inside-block.bed > HeLa-III-rep1-inside-unfiltered-2.bed
bedtools subtract -a HeLa-III-rep1-outside-unfiltered-ab.bed -b HeLa-III-rep1-outside-block.bed > HeLa-III-rep1-outside-unfiltered-2.bed
```

#### 3) Filter away low-coverage stop sites inside a peak (coverage at the site is < 1/3.5 of any sites within 30 nt of the current site)
```bash
python3 R3.py HeLa-III-rep1-inside-unfiltered-2.bed | awk '!visited[$0]++' | tr ' ' '\t' > HeLa-III-rep1-inside-unfiltered-low.bed
python3 R3.py HeLa-III-rep1-outside-unfiltered-2.bed | awk '!visited[$0]++' | tr ' ' '\t' > HeLa-III-rep1-outside-unfiltered-low.bed
bedtools subtract -a HeLa-III-rep1-inside-unfiltered-2.bed -b HeLa-III-rep1-inside-unfiltered-low.bed > HeLa-III-rep1-inside-unfiltered-3.bed
bedtools subtract -a HeLa-III-rep1-outside-unfiltered-2.bed -b HeLa-III-rep1-outside-unfiltered-low.bed > HeLa-III-rep1-outside-unfiltered-3.bed
```

### 3. filter sites based on stop ratio and coverage
use script filter1.sh for libraries from superscript III and filter2.sh for combined libraries from superscript III and IV
```bash
bash filter1.sh HeLa-IP-III-rep1.bam HeLa-III-rep1-inside-unfiltered-3.bed HeLa-III-rep1-outside-unfiltered-3.bed HeLa-III-rep1-unique.bed
bash filter2.sh HeLa-IP-III-IV-rep1.bam HeLa-III-IV-rep1-inside-unfiltered-3.bed HeLa-III-IV-rep1-outside-unfiltered-3.bed HeLa-III-IV-rep1-unique.bed
sort -k1,1 -k2,2n HeLa-III-rep1-unique.bed > HeLa-III-rep1-unique-1.bed 
sort -k1,1 -k2,2n HeLa-III-IV-rep1-unique.bed > HeLa-III-IV-rep1-unique-1.bed
```
### 4. remove stutter sites
remove  sites within 1 nt upstream and downstream of the current site whose arrested reads are at least 15% lower than the current site.
```bash
python3 Stutter1.py HeLa-III-rep1-unique-1.bed | tr ' ' '\t' >  HeLa-III-rep1-stutter-filter.bed
python3 Stutter1.py HeLa-III-IV-rep1-unique-1.bed | tr ' ' '\t' > HeLa-III-IV-rep1-stutter-filter.bed
```

remove sites within 6 nt downstream of the current site whose stop ratios are lower than the current site.
```bash
python3 Stutter2.py HeLa-III-rep1-stutter-filter.bed | tr ' ' '\t' > HeLa-III-rep1-remove.bed
python3 Stutter2.py HeLa-III-IV-rep1-stutter-filter.bed | tr ' ' '\t' > HeLa-III-IV-rep1-remove.bed
bedtools subtract -a HeLa-III-rep1-stutter-filter.bed -b HeLa-III-rep1-remove.bed > HeLa-III-rep1-stutter-filter-2.bed
bedtools subtract -a HeLa-III-IV-rep1-stutter-filter.bed -b HeLa-III-IV-rep1-remove.bed > HeLa-III-IV-rep1-stutter-filter-2.bed
```
### 5. merge all sites 
#### 1) merge sites identified from III and III+IV
```bash
bedtools subtract -a HeLa-III-IV-rep1-stutter-filter-2.bed -b HeLa-III-rep1-stutter-filter-2.bed > new.bed
cat HeLa-III-rep1-stutter-filter-2.bed new.bed | sort -k1,1 HeLa-rep1-combined.bed
```

#### 2) for quantification purpose later, obtain input reads and IP reads in libraries combining III and IV data.
```bash
# prepare the file containing the information of input reads and IP reads in libraries combining the III and IV data
cat HeLa-III-IV-rep1-inside-unfiltered-2.bed HeLa-III-IV-rep1-outside-unfiltered-2.bed > HeLa-III-IV-rep1-unfiltered-2.bed
```
```bash
#intersect with the file containing sites passing all filter steps so far. 
bedtools intersect -wa -wb -a HeLa-III-IV-rep1-unfiltered-2.bed -b HeLa-rep1-combined.bed > HeLa-rep1-combined-1.bed
awk '{print $1,$2,$3,$5,$7,$12,$13,$14,$29,$30,$31,$32,$33,$34}' HeLa-rep1-combined-1.bed | awk -v OFS="\t" '$1=$1' > HeLa-rep1-combined-2.bed
```
#### 3) Intersect two biological replicates
```bash
bedtools intersect -wa -wb -a HeLa-rep1-combined-2.bed -b HeLa-rep2-combined-2.bed > HeLa.bed
awk '!visited[$0]++' HeLa.bed | awk '{print $1,$2,$3,$4,$5,$8,$6,$7,$11,$12,$9,$10,$13,$14,$22,$20,$21,$25,$26,$23,$24,$27,$28}' | awk -v OFS="\t" '{$1=$1; print}' | tr ' ' '\t' | sort -k1,1 -k2,2n > HeLa-sort.bed
```

### 6. Further filtering based on stop ratios in the pull-down and input samples
select sites fulfiling cutoffs for average stop ratios
```bash
awk '{print $0"\t"($9+$18)/2}' HeLa-sort.bed > HeLa-input-avg.bed
awk '{print $0"\t"($10+$19)/2}' HeLa-input-avg.bed > HeLa-IP-avg.bed
python3 stop.py HeLa-IP-avg.bed > HeLa-stop-filter.bed
```

### 7. determing confidence levels and modification levels
#### 1) calculate RPM
```bash
awk '{OFS=" "; print $1,$2,$3,$4,$5,$6,($7/10.3128),($8/4.9938),$9,$10,$11,$12,$13,$14,$15,($16/11.3053),($17/5.67979),$18,$19,$20,$21,$22,$23,$24,$25}' OFS="\t" HeLa-stop-filter.bed > HeLa-RPM.bed
```
#### 2) combine with sequence context preference determined by synthetic oligos
```bash
awk '{ if($5 == "T") print $0;}' HeLa-RPM.bed > HeLa-RPM-T.bed
awk '{OFS=" "; print $1,($2-2),($3+2),$5,$6,$4,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' OFS="\t" HeLa-RPM-T.bed | tr ' ' '\t' > HeLa-extend.bed
bedtools getfasta -fi /home/Wang_yuru/Database/genome/hg38/hg38_UCSC.fa -bed HeLa-extend.bed -s -tab > HeLa-seq.bed
```
combine HeLa-seq.bed and HeLa-RPM-T.bed manually to make HeLa-seq-full.bed and sort
```bash
sort -k1,1 HeLa-seq-full.bed > HeLa-sort.bed
sort -k1,1 oligo.bed > oligo-sort.bed
```
change sequence to upper case and merge with synthetic oligo data
```bash
sed 's/^[a-z]*/\U&/' HeLa-sort.bed > HeLa-upper.bed
awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,h[$1]}' oligo-sort.bed HeLa-upper.bed > HeLa-oligo-combine.bed
```


### 10. Notes
