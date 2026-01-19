# DKC1-Knockdown-analysis-workflow
## 1. Merge paired end reads 
Paired end sequencing were used. Reads R1 and R2 were first merged, then decomplexed based on the internal barcodes. 
```bash
/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-S1_R1.fastq.gz in2=YW-S1_R2.fastq.gz out=YW_S1_merge.fastq.gz

/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-S2_R1.fastq.gz in2=YW-S2_R2.fastq.gz out=YW_S2_merge.fastq.gz

```


## 2. Sort samples based on internal barcodes.
Need to convert .gz to .fastq
```bash
gunzip YW_S1_merge.fastq.gz
gunzip YW_S2_merge.fastq.gz
```

```bash
seqkit grep -s -r -p "^TCT" YW_S1_merge.fastq -o ./sorted/HEK-sictrl-input-IV-rep1.fastq       
seqkit grep -s -r -p "^ACA" YW_S1_merge.fastq -o ./sorted/HEK-sictrl-input-IV-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S1_merge.fastq -o ./sorted/HEK-siDKC1-input-IV-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S1_merge.fastq -o ./sorted/HEK-siDKC1-input-IV-rep2.fastq
seqkit grep -s -r -p "^TCT" YW_S5_merge.fastq -o ./sorted/HEK-sictrl-IP-IV-rep1.fastq
seqkit grep -s -r -p "^ACA" YW_S5_merge.fastq -o ./sorted/HEK-sictrl-IP-IV-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S5_merge.fastq -o ./sorted/HEK-siDKC1-IP-IV-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S5_merge.fastq -o ./sorted/HEK-siDKC1-IP-IV-rep2.fastq
```
## 3. De-duplication, trim adapters and mapping
In DKC1-processing.sh, reads were first de-duplicated using clumpify.sh, and then the internal barcodes and UMI (8 nt in total on the 5' end and 5 nt in total on the 3' end) are further trimmed. The resulting reads are mapped onto human genome hg38 using hisat2.

```bash
bash DKC1-preprocessing.sh
```

## 4. Acquire information of read coverage for known candidate sites identified in HEK293T cells 
Categorize candidate modification sites identified in HEK293T cells according to the mapped strands (strand reversal due to applying applying results from R2-only mapping to the analysis using merged-read mapping)

```bash
bash DKC1-arrest.sh HEK-sictrl-input-IV-rep1.bam HEK-sictrl-IP-IV-rep1.bam ELAP-HEK-pos.bed sictrl-rep1-pos.out
bash DKC1-arrest.sh HEK-sictrl-input-IV-rep1.bam HEK-sictrl-IP-IV-rep1.bam ELAP-HEK-neg.bed sictrl-rep1-neg.out
cat sictrl-rep1-pos.out sictrl-rep1-neg.out > sictrl-rep1.out
bash DKC1-arrest.sh HEK-sictrl-input-IV-rep2.bam HEK-sictrl-IP-IV-rep2.bam ELAP-HEK-pos.bed sictrl-rep2-pos.out
bash DKC1-arrest.sh HEK-sictrl-input-IV-rep2.bam HEK-sictrl-IP-IV-rep2.bam ELAP-HEK-neg.bed sictrl-rep2-neg.out
cat sictrl-rep2-pos.out sictrl-rep2-neg.out > sictrl-rep2.out
bash DKC1-arrest.sh HEK-siDKC1-input-IV-rep1.bam HEK-siDKC1-IP-IV-rep1.bam ELAP-HEK-pos.bed siDKC1-rep1-pos.out
bash DKC1-arrest.sh HEK-siDKC1-input-IV-rep1.bam HEK-siDKC1-IP-IV-rep1.bam ELAP-HEK-neg.bed siDKC1-rep1-neg.out
cat siDKC1-rep1-pos.out siDKC1-rep1-neg.out > siDKC1-rep1.out
bash DKC1-arrest.sh HEK-siDKC1-input-IV-rep2.bam HEK-siDKC1-IP-IV-rep2.bam ELAP-HEK-pos.bed siDKC1-rep2-pos.out
bash DKC1-arrest.sh HEK-siDKC1-input-IV-rep2.bam HEK-siDKC1-IP-IV-rep2.bam ELAP-HEK-neg.bed siDKC1-rep2-neg.out
cat siDKC1-rep2-pos.out siDKC1-rep2-neg.out > siDKC1-rep2.out
bash DKC1-calculate.sh sictrl-rep1 sictrl-rep1-calculate.out
bash DKC1-calculate.sh sictrl-rep2 sictrl-rep2-calculate.out
bash DKC1-calculate.sh siDKC1-rep1 siDKC1-rep1-calculate.out
bash DKC1-calculate.sh siDKC1-rep2 siDKC1-rep2-calculate.out
```
## Process data in exce files
### 1. Calculate local RPM value by dividing read coverage at the modification site by the total number of mapped reads
### 2. Calculate average RPM values for input samples under each condition
### 3. Calculate enrichment level at each site by dividing the RPM value at the site in each IP sample by the average RPM value in input samples
### 4. Calculate fold change of enrichment level between the siDKC1 and sictrl conditions and p-value using student's t-test.

