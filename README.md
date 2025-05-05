






# Analyzing RNA structure

## 1. install viennaRNA
```bash
conda install -c bioconda viennarna
```
## 2. get chromosome locations 100 nt surrounding the modification sites
e.g. 
chr1	953761	953861	0	0	-


## 3. get RNA fasta sequence
```bash
bedtools getfasta -fi hg38_UCSC.fa -bed 5’UTR.bed -fo 5’UTR-sequence.fa -s -split
```

## 4. run viennaRNA
```bash
RNAfold -p 5’UTR-sequence.fa > 5’UTR-sequence.res
```
## 5. open the .res file in notepad++, keep the topmost prediction and save as 5'UTR-input.bed. 
Use the 'replace' function in notepad++, and replace "\n.*\n.*\n frequency.*\n.*\n.*" with ""
## 6. run python file match.py to obtain the matching status and the matching positions (if matched) for the sites in the middle of the sequence (modification site or a control site)
```bash
python match.py > 5'UTR-output.bed
```
## 7. in excel, paste the 100-nt sequences corresponding to each site. Use the MID function to get the base identity at the matched position. remember to add 1 to the numbering of each position, as python counts the first position as 0.


# DKC1-Knockdown-analysis-workflow
## Acquire information of read coverage for known candidate sites identified in HEK293T cells 
```bash
bash DKC1-preprocessing.sh
bash DKC1-arrest.sh HEK-sictrl-input-IV-rep1.bam HEK-sictrl-IP-IV-rep1.bam ELAP-HEK-all.bed sictrl-rep1.out
bash DKC1-arrest.sh HEK-sictrl-input-IV-rep2.bam HEK-sictrl-IP-IV-rep2.bam ELAP-HEK-all.bed sictrl-rep2.out
bash DKC1-arrest.sh HEK-siDKC1-input-IV-rep1.bam HEK-siDKC1-IP-IV-rep1.bam ELAP-HEK-all.bed siDKC1-rep1.out
bash DKC1-arrest.sh HEK-siDKC1-input-IV-rep2.bam HEK-siDKC1-IP-IV-rep2.bam ELAP-HEK-all.bed siDKC1-rep2.out
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
