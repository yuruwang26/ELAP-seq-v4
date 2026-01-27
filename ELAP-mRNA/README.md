# ELAP-seq-poly-A-RNA-workflow

Each replicate should contain two input libraries and an two IP libraries. The libraries are built with superscript III or IV.

## 1. processing of Fastq data (ELAP-seq-pre.sh)

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
Remove unknown chromosomes or random chromosomes manually
Save the resulting .bed file as HeLa-peaks.bed

## 3. Downstream analysis to detect pseudouridine (ELAP-seq.sh)
From step 3 to step 4, can use the command ELAP-seq.sh. This command include the following steps:

### 1. Call arrested sites inside and outside of the IP peaks. 

#### 1) Call all stop sites inside and outside of IP peaks
This step uses the script arrest.sh
```bash
bash arrest.sh HeLa-input-III-rep1.bam HeLa-IP-III-rep1.bam HeLa-peaks.bed HeLa-III-rep1-inside.out HeLa-III-rep1-outside.out
bash arrest.sh HeLa-input-III-IV-rep1.bam HeLa-IP-III-IV-rep1.bam HeLa-peaks.bed HeLa-III-IV-rep1-inside.out HeLa-III-IV-rep1-outside.out
```
#### 2) Calculate arrest rate of each site and assign the originality of the site (i.e., whether it is from the library built with superscript III or combined libraries built with superscript III and IV, and whether it is inside or outside of the IP peak)
This step uses scripts calculate_III.sh for the III librairy and calculate_III_IV.sh for the III+IV library
```bash
bash calculate_III.sh HeLa-IP-III-rep1.bam HeLa-III-rep1-inside.out HeLa-III-rep1-outside.out HeLa-III-rep1-inside-unfiltered.bed HeLa-III-rep1-outside-unfiltered.bed HeLa-III-in HeLa-III-out HeLa-III-rep1-unfiltered.bed
bash calculate_III_IV.sh HeLa-IP-III-IV-rep1.bam HeLa-III-IV-rep1-inside.out HeLa-III-IV-rep1-outside.out HeLa-III-IV-rep1-inside-unfiltered.bed HeLa-III-IV-rep1-outside-unfiltered.bed HeLa-III-IV-in HeLa-III-IV-out HeLa-III-rep1-unfiltered.bed
```

### 2. Remove backgrounds 
This process removes background originated from multiple mapping and diminished processivity of the RT enzyme past a major modification site
#### 1) determine the main nucleoside for sites mapped with multiple identities
```bash
python3 Rm_bg_1.py HeLa-III-rep1-inside-unfiltered.bed | tr ' ' '\t' > HeLa-III-rep1-inside-unfiltered-ab.bed
python3 Rm_bg_1.py HeLa-III-rep1-outside-unfiltered.bed | tr ' ' '\t' > HeLa-III-rep1-outside-unfiltered-ab.bed
```
#### 2) remove regions that are covered by reads that all share the same start and end mapping position
```bash
python3 Rm_bg_2.py HeLa-III-rep1-inside-unfiltered-ab.bed | tr ' ' '\t' > HeLa-III-rep1-inside-block.bed
python3 Rm_bg_2.py HeLa-III-rep1-outside-unfiltered-ab.bed | tr ' ' '\t' > HeLa-III-rep1-outside-block.bed
bedtools subtract -a HeLa-III-rep1-inside-unfiltered-ab.bed -b HeLa-III-rep1-inside-block.bed > HeLa-III-rep1-inside-unfiltered-2.bed
bedtools subtract -a HeLa-III-rep1-outside-unfiltered-ab.bed -b HeLa-III-rep1-outside-block.bed > HeLa-III-rep1-outside-unfiltered-2.bed
```
#### 3) Filter away low-coverage stop sites within 50 nt downstream a major stop site
```bash
python3 Rm_bg_3.py HeLa-III-rep1-inside-unfiltered-2.bed | awk '!visited[$0]++' | tr ' ' '\t' > HeLa-III-rep1-inside-unfiltered-low.bed
python3 Rm_bg_3.py HeLa-III-rep1-outside-unfiltered-2.bed | awk '!visited[$0]++' | tr ' ' '\t' > HeLa-III-rep1-outside-unfiltered-low.bed
bedtools subtract -a HeLa-III-rep1-inside-unfiltered-2.bed -b HeLa-III-rep1-inside-unfiltered-low.bed > HeLa-III-rep1-inside-unfiltered-3.bed
bedtools subtract -a HeLa-III-rep1-outside-unfiltered-2.bed -b HeLa-III-rep1-outside-unfiltered-low.bed > HeLa-III-rep1-outside-unfiltered-3.bed
```

### 3. Filter sites based on stop ratios, stopped reads

#### 1) select for sites that have stop ratio >=0.1 in the pull-down library, stop ratio (pull-down) -stop ratio (input) >= 0.05 (for data from SSIII alone) or 0.1 (for combined data of SSIII and SSIV),  and are covered by at least 1 uniquely mapped reads
use script filter.sh 
```bash
bash filter_III.sh HeLa-IP-III-rep1.bam HeLa-III-rep1-inside-unfiltered-3.bed HeLa-III-rep1-outside-unfiltered-3.bed HeLa-III-rep1-unique.bed && sort -k1,1 -k2,2n HeLa-III-rep1-unique.bed > HeLa-III-rep1-unique-1.bed
bash filter_III_IV.sh HeLa-IP-III-IV-rep1.bam HeLa-III-IV-rep1-inside-unfiltered-3.bed HeLa-III-IV-rep1-outside-unfiltered-3.bed HeLa-III-IV-rep1-unique.bed && sort -k1,1 -k2,2n HeLa-III-IV-rep1-unique.bed > HeLa-III-IV-rep1-unique-1.bed
```
#### 2). remove stutter sites
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
#### 3). select for sites with at least 5 stopped reads
```bash
awk '{ if($10 >4) print $0;}' HeLa-III-rep1-stutter-filter-2.bed > HeLa-III-rep1-filter1.bed
awk '{ if($10 >4) print $0;}' HeLa-III-rep2-stutter-filter-2.bed > HeLa-III-rep2-filter1.bed
awk '{ if($10 >4) print $0;}' HeLa-III-IV-rep1-stutter-filter-2.bed > HeLa-III-IV-rep1-filter1.bed
awk '{ if($10 >4) print $0;}' HeLa-III-IV-rep2-stutter-filter-2.bed > HeLa-III-IV-rep2-filter1.bed
```

#### 4). Remove sites whose stop ratios are >= 0.1 in the input, (stop ratio in pull-down)/(stop ratio in input) are < 3, and stopped reads in the input are >=3 
```bash
awk '($14 <= 0.1) || ($8 < 3) || ($15 / $14 >= 3)' HeLa-III-rep1-filter1.bed > HeLa-III-rep1-filter2.bed
```
#### 5). Optional: for evaluainge reproducibility, focus on sites that are covered by at least five reads in one other replicate and require that stop ratio * stopped reads is >=1.5 before intersecting replicates
```bash
awk '($13 >=5 && $10*$15 > = 1.5)' HeLa-III-rep1-filter2.bed > HeLa-III-rep1-filter3.bed
bedtools intersect -a HeLa-III-rep1-filter3.bed HeLa-III-rep2-filter3.bed > HeLa-III-rep1-rep2.bed
bedtools intersect -a HeLa-III-rep1-filter3.bed HeLa-III-rep3-filter3.bed > HeLa-III-rep1-rep3.bed
bedtools subtract -a HeLa-III-rep1-rep2.bed -b HeLa-III-rep1-rep3.bed > tmp.bed
cat HeLa-III-rep1-rep3.bed tmp.bed > HeLa-III-rep1-filter4.bed
```
## 4 Intersect two biological replicates and further filter (if using superscript III data alone)


### 1) Intersect two biological replicates
### 2) Select for sites whose average value of stop ratio * stopped reads between the two pull-down replicates is >=1.5.


## 4 Intersect two biological replicates and further filter (if using superscript III and IV data)

### 1) combine sites identified from III and new sites identified from III+IV
```bash
bedtools subtract -a HeLa-III-IV-rep1-filter2.bed -b HeLa-III-rep1-filter2.bed > new.bed
cat HeLa-III-rep1-filter2.bed new.bed | sort -k1,1 > HeLa-rep1-combined.bed
```
### 2) for quantification purpose later, obtain input reads and IP reads in libraries combining III and IV data & cleanup the table
```bash
# prepare the file containing the information of input reads and IP reads in libraries combining the III and IV data
cat HeLa-III-IV-rep1-inside-unfiltered-2.bed HeLa-III-IV-rep1-outside-unfiltered-2.bed > HeLa-III-IV-rep1-unfiltered-2.bed
```
```bash
#intersect with the file containing sites passing all filter steps so far. 
bedtools intersect -wa -wb -a HeLa-III-IV-rep1-unfiltered-2.bed -b HeLa-rep1-combined.bed > HeLa-rep1-combined-1.bed
```

```bash
awk '{print $1,$2,$3,$5,$7,$12,$13,$14,$29,$30,$31,$32,$33,$34}' HeLa-rep1-combined-filtered.bed | awk -v OFS="\t" '$1=$1' > HeLa-rep1-combined-2.bed
```

### 3) Intersect two biological replicates
```bash
bedtools intersect -wa -wb -a HeLa-rep1-combined-2.bed -b HeLa-rep2-combined-2.bed > HeLa.bed
awk '!visited[$0]++' HeLa.bed | awk '{print $1,$2,$3,$4,$5,$8,$6,$7,$11,$12,$9,$10,$13,$14,$22,$20,$21,$25,$26,$23,$24,$27,$28}' | awk -v OFS="\t" '{$1=$1; print}' | tr ' ' '\t' | sort -k1,1 -k2,2n > HeLa-sort.bed
```

### 4) Select for sites whose average value of stop ratio * stopped reads between the two pull-down replicates is >=1.5.


## 5. post-processing : determing confidence levels and modification levels
### 1) calculate RPM
```bash
awk '{OFS=" "; print $1,$2,$3,$4,$5,$6,($7/10.3128),($8/4.9938),$9,$10,$11,$12,$13,$14,$15,($16/11.3053),($17/5.67979),$18,$19,$20,$21,$22,$23,$24,$25}' OFS="\t" HeLa-filter-T.bed > HeLa-RPM.bed
```
### 2) combine with sequence context preference determined by synthetic oligos
Obtain sequence context surrounding the modification site and make all upper case
```bash
awk '{OFS=" "; print $1,($2-2),($3+2),$5,$6,$4,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' OFS="\t" HeLa-RPM.bed | tr ' ' '\t' > HeLa-extend.bed
bedtools getfasta -fi /home/Wang_yuru/Database/genome/hg38/hg38_UCSC.fa -bed HeLa-extend.bed -s -tab > HeLa-seq.bed
sed 's/^[a-z]*/\U&/' HeLa-seq.bed > HeLa-upper.bed

```
Join HeLa-upper.bed and HeLa-RPM.bed manually to make HeLa-seq-full.bed and sort
```bash
sort -k1,1 HeLa-seq-full.bed > HeLa-sort.bed
sort -k1,1 oligo.bed > oligo-sort.bed
```
merge with synthetic oligo data using code or using VLOOK function in excel
```bash
awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,h[$1]}' oligo-sort.bed HeLa-sort.bed > HeLa-oligo-combine.bed
```
Confidence level:
The highest-confidence sites are defined as sites having IP stop ratio >= 0.3 in at least replicates and also identified in all three replicates.
The higher-confidence sites are defined as sites either having IP stop ratio >= 0.3 in two replicate or identified in all three replicates.
The lower-confidence sites are defined as sites having IP stop ratio < 0.3 and identified in only two replicates.

Modification level:
Relative modification level = Ln(enrichment at the modification site/enrichment of the oligo with the same sequence context)
Sites are defined as highly modified when relative modification levels are >= 2.3 for sites in ΨU context and >= 3.5 for sites in non-ΨU context.
Sites are defined as moderately modified when relative modification levels are >= 1.6 for sites in ΨU context and >= 1.8 for sites in non-ΨU context.
Sites are defined as lowly modified when relative modification levels are < 1.6 for sites in ΨU context and < 1.8 for sites in non-ΨU context.
