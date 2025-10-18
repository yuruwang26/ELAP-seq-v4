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
## 5. open the .res file in notepad++, keep the topmost prediction and save as 5'UTR-input.txt. 
Use the 'replace' function in notepad++, and replace "\r?\n.*\r?\n.*\r?\n frequency.*\r?\n.*\r?\n.*" with ""
## 6. run python file match.py to obtain the matching status and the matching positions (if matched) for the sites in the middle of the sequence (modification site or a control site)
```bash
python match.py > 5'UTR-output.bed
```
## 7. in excel, paste the 100-nt sequences corresponding to each site. Use the MID function to get the base identity at the matched position. remember to add 1 to the numbering of each position, as python counts the first position as 0.
