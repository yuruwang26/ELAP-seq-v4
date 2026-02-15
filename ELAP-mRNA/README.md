# ELAP-seq-poly-A-RNA-workflow

Each replicate should contain two input libraries and an two IP libraries. The libraries are built with superscript III or IV.

## 1. Processing of FASTQ data 

Only reads R2 is used. After trimming UMI, the begnning of R2 is the RT stop site, which indicates modification. Can otherwise construct the libraries in a way that single end sequencing is sufficient.
This process includes five steps:

### 1) trim adapter

```bash
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-rep1-III-input-cutadapt.fq.gz HeLa-input-III-rep1_R2.fq.gz >> adaptorTrim.log
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-rep1-IV-input-cutadapt.fq.gz HeLa-input-IV-rep1_R2.fq.gz >> adaptorTrim.log
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-rep1-III-IP-cutadapt.fq.gz HeLa-IP-III-rep1_R2.fq.gz >> adaptorTrim.log
cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  -o HeLa-rep1-IV-IP-cutadapt.fq.gz HeLa-IP-IV-rep1_R2.fq.gz >> adaptorTrim.log
```

### 2) remove duplicates

```bash
~/Tools/bbmap/clumpify.sh in=HeLa-rep1-III-input-cutadapt.fq.gz out=HeLa-rep1-III-input-dedupe.fq.gz dedupe >> duplicates_removal_1.log
~/Tools/bbmap/clumpify.sh in=HeLa-rep1-IV-input-cutadapt.fq.gz out=HeLa-rep1-IV-input-dedupe.fq.gz dedupe >> duplicates_removal_1.log
~/Tools/bbmap/clumpify.sh in=HeLa-rep1-III-IP-cutadapt.fq.gz out=HeLa-rep1-III-IP-dedupe.fq.gz dedupe >> duplicates_removal_1.log
~/Tools/bbmap/clumpify.sh in=HeLa-rep1-IV-IP-cutadapt.fq.gz out=HeLa-rep1-IV-IP-dedupe.fq.gz dedupe >> duplicates_removal_1.log
```

### 3) trim UMI
There is a 6-nt UMI sequence at the 5′ end and a 5-nt UMI sequence at the 3′ end. For Superscript IV data, trim one additional nucleotide from the 5′ end, as Superscript IV tends to elongate by one extra nucleotide past the modification.
```bash
conda activate cutadaptenv
```

```bash
cutadapt -u 6 -o tmp.trimmed.fastq HeLa-rep1-III-input-dedupe.fq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-rep1-III-input-trim.fq.gz tmp.trimmed.fastq
```

```bash
cutadapt -u 7 -o tmp.trimmed.fastq HeLa-rep1-IV-input-dedupe.fq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-rep1-IV-input-trim.fq.gz tmp.trimmed.fastq
```
```bash
cutadapt -u 6 -o tmp.trimmed.fastq HeLa-rep1-III-IP-dedupe.fq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-rep1-III-IP-trim.fq.gz tmp.trimmed.fastq
```

```bash
cutadapt -u 7 -o tmp.trimmed.fastq HeLa-rep1-IV-IP-dedupe.fq.gz
cutadapt -u -5 -q 10,10 -m 19 -o HeLa-rep1-IV-IP-trim.fq.gz tmp.trimmed.fastq
```
```bash
conda deactivate
```
### 4) map reads to the genome

```bash
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file YW_align_summary -p 4 -U HeLa-rep1-III-input-trim.fq.gz |samtools view -bS |samtools sort -o HeLa-rep1-III-input.bam
samtools index HeLa-rep1-III-input.bam HeLa-rep1-III-input.bai
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file YW_align_summary -p 4 -U HeLa-rep1-IV-input-trim.fq.gz |samtools view -bS |samtools sort -o HeLa-rep1-IV-input.bam
samtools index HeLa-rep1-IV-input.bam HeLa-rep1-IV-input.bai
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file YW_align_summary -p 4 -U HeLa-rep1-III-IP-trim.fq.gz |samtools view -bS |samtools sort -o HeLa-rep1-III-IP.bam
samtools index HeLa-rep1-III-IP.bam HeLa-rep1-III-IP.bai
hisat2 -x /home/yuruwang/Database/genome/hg38/hg38_UCSC --known-splicesite-infile /home/yuruwang/Database/genome/hg38/hisat2_splice_sites.txt --rna-strandness F --no-softclip --summary-file YW_align_summary -p 4 -U HeLa-rep1-IV-IP-trim.fq.gz |samtools view -bS |samtools sort -o HeLa-rep1-IV-IP.bam
samtools index HeLa-rep1-IV-IP.bam HeLa-rep1-IV-IP.bai
```

### 5) combine .bam files from superscript III and superscript IV in order to rescue additional sites that are missed by the superscript III library due to low read coverage
```bash
samtools merge HeLa-rep1-III-IV-input.bam HeLa-rep1-III-input.bam HeLa-rep1-IV-input.bam
samtools merge HeLa-rep1-III-IV-IP.bam HeLa-rep1-III-IP.bam HeLa-rep1-IV-IP.bam
samtools index HeLa-rep1-III-IV-input.bam HeLa-rep1-III-IV-input.bai
samtools index HeLa-rep1-III-IV-IP.bam HeLa-rep1-III-IV-IP.bai
```



## 2. Call IP peaks
### 1) MACS2 is used to call IP peaks. 
In the final list, we want to report whether a site is inside an enriched IP peak or not

```bash
macs2 callpeak -t HeLa-rep1-III-IP.bam -c HeLa-rep1-III-input.bam -n test_t2 -f BAM -g 994080837 -q 0.01 --slocal 1000 --extsize 150 --nomodel --keep-dup all --call-summits --outdir HeLa-III-peadDir
macs2 callpeak -t HeLa-rep1-III-IV-IP.bam -c HeLa-rep1-III-IV-input.bam -n test_t2 -f BAM -g 994080837 -q 0.01 --slocal 1000 --extsize 150 --nomodel --keep-dup all --call-summits --outdir HeLa-III-IV-peadDir
```
### 2) obtain regions covered by IP peaks 
Remove unknown chromosomes or random chromosomes manually
Save the resulting .bed file as HeLa-peaks.bed




## 3. Downstream analysis to detect pseudouridine (ELAP-seq.sh)
This step can use the command ELAP-seq.sh. 

Usage: 
```bash
ELAP-seq.sh <rep1-III-IP.bam> <rep1-III-IV-IP.bam> <rep1-III-input.bam> <rep1-III-IV-input.bam> <rep1-peaks.bed> <rep1-name> 
ELAP-seq.sh <rep2-III-IP.bam> <rep2-III-IV-IP.bam> <rep2-III-input.bam> <rep2-III-IV-input.bam> <rep2-peaks.bed> <rep2-name>
```

This command include the following steps:

### 1. Call arrested sites inside and outside of the IP peaks. 

#### 1) Call all stop sites inside and outside of IP peaks
This step uses the script arrest.sh
```bash
bash arrest.sh HeLa-rep1-III-input.bam HeLa-rep1-III-IP.bam HeLa-peaks.bed HeLa-rep1-III-inside.out HeLa-rep1-III-outside.out
bash arrest.sh HeLa-rep1-III-IV-input.bam HeLa-rep1-III-IV-IP.bam HeLa-peaks.bed HeLa-rep1-III-IV-inside.out HeLa-rep1-III-IV-outside.out
```
#### 2) Calculate arrest rate of each site and assign the originality of the site (i.e., whether it is inside or outside of an IP peak, and whether it is identified from the library built with superscript III or combined data of libraries built with superscript III and IV)
This step uses scripts calculate_III.sh for the III librairy and calculate_III_IV.sh for the III+IV library. In the output, "in" means inside an IP peak, "out" means outside an IP peak. "III" means the site was called from the library built with superscript III, and "III_IV" means the site was called from combined libraries built with superscript III and IV.
```bash
bash calculate_III.sh HeLa-rep1-III-IP.bam HeLa-rep1-III-inside.out HeLa-rep1-III-outside.out HeLa-rep1-III-inside-unfiltered.bed HeLa-rep1-III-outside-unfiltered.bed HeLa-III-in HeLa-III-out HeLa-rep1-III-unfiltered.bed
bash calculate_III_IV.sh HeLa-rep1-III-IV-IP.bam HeLa-rep1-III-IV-inside.out HeLa-rep1-III-IV-outside.out HeLa-rep1-III-IV-inside-unfiltered.bed HeLa-rep1-III-IV-outside-unfiltered.bed HeLa-III-IV-in HeLa-III-IV-out HeLa-rep1-III-unfiltered.bed
```

### 2. Remove false positives
This process removes noises originated from multiple mapping and diminished processivity of the RT enzyme past a major modification site
#### 1) determine the main nucleoside for sites mapped with multiple identities
```bash
python3 Rm_bg_1.py HeLa-rep1-III-inside-unfiltered.bed | tr ' ' '\t' > HeLa-rep1-III-inside-unfiltered-ab.bed
python3 Rm_bg_1.py HeLa-rep1-III-outside-unfiltered.bed | tr ' ' '\t' > HeLa-rep1-III-outside-unfiltered-ab.bed
```
#### 2) remove regions that are covered by reads that all share the same start and end mapping position
This is likely caused by off-target mapping of reads whose originality is another region in the genome containing a same sequence as the current region 
```bash
python3 Rm_bg_2.py HeLa-rep1-III-inside-unfiltered-ab.bed | tr ' ' '\t' > HeLa-rep1-III-inside-block.bed
python3 Rm_bg_2.py HeLa-rep1-III-outside-unfiltered-ab.bed | tr ' ' '\t' > HeLa-rep1-III-outside-block.bed
bedtools subtract -a HeLa-rep1-III-inside-unfiltered-ab.bed -b HeLa-rep1-III-inside-block.bed > HeLa-rep1-III-inside-unfiltered-2.bed
bedtools subtract -a HeLa-rep1-III-outside-unfiltered-ab.bed -b HeLa-rep1-III-outside-block.bed > HeLa-rep1-III-outside-unfiltered-2.bed
```
#### 3) Filter away low-coverage stop sites within 50 nt downstream a major stop site
These stop sites tend to be false positive sites caused by diminished processivity of the RT enzyme reading beyond the Ψ modification
```bash
python3 Rm_bg_3.py HeLa-rep1-III-inside-unfiltered-2.bed | awk '!visited[$0]++' | tr ' ' '\t' > HeLa-rep1-III-inside-unfiltered-low.bed
python3 Rm_bg_3.py HeLa-rep1-III-outside-unfiltered-2.bed | awk '!visited[$0]++' | tr ' ' '\t' > HeLa-rep1-III-outside-unfiltered-low.bed
awk '!visited[$0]++' HeLa-rep1-III-inside-unfiltered-low.bed > HeLa-rep1-III-inside-unfiltered-low-1.bed
awk '!visited[$0]++' HeLa-rep1-III-outside-unfiltered-low.bed > HeLa-rep1-III-outside-unfiltered-low-1.bed
cat HeLa-rep1-III-inside-unfiltered-low-1.bed | tr ' ' '\t' > HeLa-rep1-III-inside-unfiltered-low-2.bed
cat HeLa-rep1-III-outside-unfiltered-low-1.bed | tr ' ' '\t' > HeLa-rep1-III-outside-unfiltered-low-2.bed
bedtools subtract -a HeLa-rep1-III-inside-unfiltered-2.bed -b HeLa-rep1-III-inside-unfiltered-low-2.bed > HeLa-rep1-III-inside-unfiltered-3.bed
bedtools subtract -a HeLa-rep1-III-outside-unfiltered-2.bed -b HeLa-rep1-III-outside-unfiltered-low-2.bed > HeLa-rep1-III-outside-unfiltered-3.bed
```

### 3. Filter sites based on stop ratios and stopped reads

#### 1) select for sites that have stop ratio >=0.1 in the pull-down library, stop ratio (pull-down) -stop ratio (input) >= 0.05 (for data from SSIII alone) or 0.1 (for combined data of SSIII and SSIV),  and are covered by at least 1 uniquely mapped reads
use script filter.sh 
```bash
bash filter_III.sh HeLa-rep1-III-IP.bam HeLa-rep1-III-inside-unfiltered-3.bed HeLa-rep1-III-outside-unfiltered-3.bed HeLa-rep1-III-unique.bed && sort -k1,1 -k2,2n HeLa-rep1-III-unique.bed > HeLa-rep1-III-unique-1.bed
bash filter_III_IV.sh HeLa-rep1-III-IV-IP.bam HeLa-rep1-III-IV-inside-unfiltered-3.bed HeLa-rep1-III-IV-outside-unfiltered-3.bed HeLa-rep1-III-IV-unique.bed && sort -k1,1 -k2,2n HeLa-rep1-III-IV-unique.bed > HeLa-rep1-III-IV-unique-1.bed
```
#### 2). remove stutter sites
Remove  sites within 1 nt upstream and downstream of the current site whose arrested reads are at least 15% lower than the current site. We select the major arrest sites this way.
```bash
python3 Stutter1.py HeLa-rep1-III-unique-1.bed | tr ' ' '\t' >  HeLa-rep1-III-stutter-filter.bed
python3 Stutter1.py HeLa-rep1-III-IV-unique-1.bed | tr ' ' '\t' > HeLa-rep1-III-IV-stutter-filter.bed
```

Remove sites within 6 nt downstream of the current site whose stop ratios are lower than the current site.This continues to remove noises nearby the major arrest sites.
```bash
python3 Stutter2.py HeLa-rep1-III-stutter-filter.bed | tr ' ' '\t' > HeLa-rep1-III-remove.bed
python3 Stutter2.py HeLa-rep1-III-IV-stutter-filter.bed | tr ' ' '\t' > HeLa-rep1-III-IV-remove.bed
bedtools subtract -a HeLa-rep1-III-stutter-filter.bed -b HeLa-rep1-III-remove.bed > HeLa-rep1-III-stutter-filter-2.bed
bedtools subtract -a HeLa-rep1-III-IV-stutter-filter.bed -b HeLa-rep1-III-IV-remove.bed > HeLa-rep1-III-IV-stutter-filter-2.bed
```
#### 3). select for sites with at least 5 stopped reads
```bash
awk '{ if($10 >4) print $0;}' HeLa-rep1-III-stutter-filter-2.bed > HeLa-rep1-III-filter1.bed
awk '{ if($10 >4) print $0;}' HeLa-rep2-III-stutter-filter-2.bed > HeLa-rep2-III-filter1.bed
awk '{ if($10 >4) print $0;}' HeLa-rep1-III-IV-stutter-filter-2.bed > HeLa-rep1-III-IV-filter1.bed
awk '{ if($10 >4) print $0;}' HeLa-rep2-III-IV-stutter-filter-2.bed > HeLa-rep2-III-IV-filter1.bed
```

#### 4). Remove sites whose stop ratios are > 0.1 in the input, stopped reads in the input are >=3, and (stop ratio in pull-down)/(stop ratio in input) are < 3  
Sites fulfiling these three cutoffs tend to be false positives since the stop signature is apparently present in the input samples.
```bash
awk '($14 <= 0.1) || ($8 < 3) || ($15 / $14 >= 3)' HeLa-rep1-III-filter1.bed > HeLa-rep1-III-filter2.bed
```
#### 5). *Optional: for evaluainge reproducibility, focus on sites that are covered by at least five reads in one other replicate and require that stop ratio * stopped reads is >=1.5 before intersecting replicates
```bash
awk '($13 >=5 && $10*$15 >= 1.5)' HeLa-rep1-III-filter2.bed > HeLa-rep1-III-filter3.bed
bedtools intersect -a HeLa-rep1-III-filter3.bed HeLa-rep2-III-filter3.bed > HeLa-rep1-rep2-III.bed
bedtools intersect -a HeLa-rep1-III-filter3.bed HeLa-rep3-III-filter3.bed > HeLa-rep1-rep3-III.bed
bedtools subtract -a HeLa-rep1-rep2-III.bed -b HeLa-rep1-rep3-III.bed > tmp.bed
cat HeLa-rep1-rep3-III.bed tmp.bed > HeLa-rep1-III-filter4.bed
```

### 4. Combine sites identified from superscript III alone and the new sites identified using the combined data of libraries built with Superscript III and Superscript IV 
```bash
bedtools subtract -a HeLa-rep1-III-IV-filter2.bed -b HeLa-rep1-III-filter2.bed > new.bed
cat HeLa-rep1-III-filter2.bed new.bed | sort -k1,1 > HeLa-rep1-combined.bed
```
The resulting file contains: chr start end pvalue strand arrest_score ref Input_arrest_rep1 Input_readthrough_rep1 IP_arrest_rep1 IP_readthrough_rep1 Input_total_count_rep1 IP_total_count_rep1 Input_stop_ratio_rep1 IP_stop_ratio_rep1 peak_rep1 sample_origin_rep1 

## 4 Intersect two biological replicates and further filter

### 1. If looking at sites identified by superscript III data alone
#### 1) Intersect two biological replicates
```bash
bedtools intersect -wa -wb -a HeLa-rep1-III-filter2.bed -b HeLa-rep2-III-filter2.bed > HeLa.bed
awk '!visited[$0]++' HeLa.bed | awk '{print $1,$2,$3,$5,$7,$10,$12,$13,$14,$15,$16,$17,$27,$29,$30,$31,$32,$33,$34}' | awk -v OFS="\t" '{$1=$1; print}' | tr ' ' '\t' | sort -k1,1 -k2,2n > HeLa-sort.bed
```
The resulting file contains: chr start end strand ref IP_arrest_rep1 Input_total_count_rep1 IP_total_count_rep1 Input_stop_ratio_rep1 IP_stop_ratio_rep1 peak_rep1 sample_origin_rep1 IP_arrest_rep2 Input_total_count_rep2 IP_total_count_rep2 Input_stop_ratio_rep2 IP_stop_ratio_rep2 peak_rep2 sample_origin_rep2
#### 2) Select for sites whose average value of stop ratio * stopped reads between the two pull-down replicates is >=1.5.
This imposes a more stringent requirement on the number of stopped reads when the stop ratio in the pull-down sample is lower than 30%
```bash
awk '($6*$10 + $13*$7)/2 >=1.5 ' HeLa-sort.bed > HeLa-filter.bed
```
#### 3) Select for sites whose stop locate at T
```bash
awk '{ if($5 == "T") print $0;}' HeLa-filter.bed > HeLa-T.bed
```
### 2. If looking at combined sites identified by superscript III-alone data and combined data of superscript III and IV

#### 1) Intersect two biological replicates
```bash
bedtools intersect -wa -wb -a HeLa-rep1-combined.bed -b HeLa-rep2-combined.bed > HeLa.bed
awk '!visited[$0]++' HeLa.bed | awk '{print $1,$2,$3,$5,$7,$10,$12,$13,$14,$15,$16,$17,$27,$29,$30,$31,$32,$33,$34}' | awk -v OFS="\t" '{$1=$1; print}' | tr ' ' '\t' | sort -k1,1 -k2,2n > HeLa-sort.bed
```
The resulting file contains: chr start end strand ref IP_arrest_rep1 Input_total_count_rep1 IP_total_count_rep1 Input_stop_ratio_rep1 IP_stop_ratio_rep1 peak_rep1 sample_origin_rep1 IP_arrest_rep2 Input_total_count_rep2 IP_total_count_rep2 Input_stop_ratio_rep2 IP_stop_ratio_rep2 peak_rep2 sample_origin_rep2
#### 2) Select for sites whose average value of stop ratio * stopped reads between the two pull-down replicates is >=1.5.
```bash
awk '($6*$10 + $13*$7)/2 >=1.5 ' HeLa-sort.bed > HeLa-filter.bed
```
#### 3) Select for sites whose stop locate at T
```bash
awk '{ if($5 == "T") print $0;}' HeLa-filter.bed > HeLa-T.bed
```



## 5. post-processing: determing confidence levels and modification levels
### 1) Determine confidence level for each site:
The highest-confidence sites are defined as sites having IP stop ratio >= 0.3 in at least replicates and also identified in all three replicates.
The higher-confidence sites are defined as sites either having IP stop ratio >= 0.3 in two replicate or identified in all three replicates.
The lower-confidence sites are defined as sites having IP stop ratio < 0.3 and identified in only two replicates.
### 2). For quantification, use the combined data of libraries built with SuperScript III and SuperScript IV to achieve the best read coverage for each site in both input and IP samples 
####1)
```bash
# Prepare the file containing the information of input reads and IP reads obtained from the combined data of libraries built with SuperScript III and SuperScript IV.
cat HeLa-rep1-III-IV-inside-unfiltered-2.bed HeLa-rep1-III-IV-outside-unfiltered-2.bed > HeLa-rep1-III-IV-unfiltered-2.bed
cat HeLa-rep2-III-IV-inside-unfiltered-2.bed HeLa-rep2-III-IV-outside-unfiltered-2.bed > HeLa-rep2-III-IV-unfiltered-2.bed
```
```bash
#intersect with candidate modification sites passing all filter steps so far. 
bedtools intersect -wa -wb -a HeLa-III-IV-rep1-unfiltered-2.bed -b HeLa-rep1-combined.bed > HeLa-rep1-combined-1.bed
```

```bash
# tidy up the table
awk '{print $1,$2,$3,$5,$7,$12,$13,$14,$15,$29,$30,$31,$32,$33,$34}' HeLa-rep1-combined-1.bed | awk -v OFS="\t" '$1=$1' > HeLa-rep1-combined-2.bed
```
The resulting file format: chr start end strand ref reads_input_III_IV_rep1 reads_IP_III_IV_rep1 stop_ratio_input_III_IV_rep1 stop_ratio_IP_III_IV_rep1 reads_input_rep1 reads_IP_rep1 stop_ratio_input_rep1 stop_ratio_IP_rep1 peak_rep1 sample_rep1
#### 2) using data from III+IV data, calculate RPM by dividing the sum total reads at the site in the sample by the total number of mapped reads of the entire sample and calculate enrichment levels by dividing RPM in IP by RPM in input.
#### 3) combine with sequence context preference determined by synthetic oligos
#### 4) Obtain 5-nt sequence context surrounding the modification site 
#### 5) merge with synthetic oligo data with the corresponding sequence context using using VLOOK function in excel
#### 6) calculate adjusted enrichment levels by dividing enrichment levels at the site by the enrichment level of the corresponding sequence context 
#### 7) Calculate relative modification levels using Ln(adjusted enrichment levels)



Modification level:
Relative modification level = Ln(enrichment at the modification site/enrichment of the oligo with the same sequence context)
Sites are defined as highly modified when relative modification levels are >= 2.3 for sites in ΨU context and >= 3.5 for sites in non-ΨU context.
Sites are defined as moderately modified when relative modification levels are >= 1.6 for sites in ΨU context and >= 1.8 for sites in non-ΨU context.
Sites are defined as lowly modified when relative modification levels are < 1.6 for sites in ΨU context and < 1.8 for sites in non-ΨU context.
