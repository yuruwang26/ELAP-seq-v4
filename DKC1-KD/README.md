# DKC1-Knockdown-analysis-workflow
## Categorize candidate modification sites identified in HEK293T cells according to the mapped strands (strand reversal due to applying applying results from R2-only mapping to the analysis using merged-read mapping)
## Acquire information of read coverage for known candidate sites identified in HEK293T cells 
```bash
bash DKC1-preprocessing.sh
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
### 1. Calculate RPM value
### 2. Calculate average RPM values for input samples under each condition
### 3. Calculate enrichment level at each site by dividing the RPM value at the site in each IP sample by the average RPM value in input samples
### 4. Calculate fold change of enrichment level between the siDKC1 and sictrl conditions and p-value using one-sided student's t-test.

