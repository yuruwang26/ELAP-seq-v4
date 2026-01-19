

/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-S1_R1.fastq.gz in2=YW-S1_R2.fastq.gz out=YW_S1_merge.fastq.gz

/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-S2_R1.fastq.gz in2=YW-S2_R2.fastq.gz out=YW_S2_merge.fastq.gz

gunzip YW_S1_merge.fastq.gz
gunzip YW_S2_merge.fastq.gz

###need to convert .gz to .fastq
seqkit grep -s -r -p "^TCT" YW_S1_merge.fastq -o ./sorted/HEK-sictrl-input-IV-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S1_merge.fastq -o ./sorted/HEK-sictrl-input-IV-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S1_merge.fastq -o ./sorted/HEK-siDKC1-input-IV-rep1.fastq
seqkit grep -s -r -p "^ACA" YW_S1_merge.fastq -o ./sorted/HEK-siDKC1-input-IV-rep2.fastq

seqkit grep -s -r -p "^TCT" YW_S5_merge.fastq -o ./sorted/HEK-sictrl-IP-IV-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S5_merge.fastq -o ./sorted/HEK-sictrl-IP-IV-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S5_merge.fastq -o ./sorted/HEK-siDKC1-IP-IV-rep1.fastq
seqkit grep -s -r -p "^ACA" YW_S5_merge.fastq -o ./sorted/HEK-siDKC1-IP-IV-rep2.fastq



