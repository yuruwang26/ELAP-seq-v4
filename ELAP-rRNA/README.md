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
bedtools intersect -wa -wb -a rRNA-rep1-filter.bed -b rRNA-rep2-filter.bed > rRNA-out.bed
```
