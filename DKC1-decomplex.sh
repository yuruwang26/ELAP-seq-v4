

/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-6s-1120-S3_R1.fastq.gz in2=YW-6s-1120-S3_R2.fastq.gz out=YW_S3_merge.fastq.gz

/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-6s-1120-S4_R1.fastq.gz in2=YW-6s-1120-S4_R2.fastq.gz out=YW_S4_merge.fastq.gz

/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-6s-1120-S5_R1.fastq.gz in2=YW-6s-1120-S5_R2.fastq.gz out=YW_S5_merge.fastq.gz

/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-6s-1120-S6_R1.fastq.gz in2=YW-6s-1120-S6_R2.fastq.gz out=YW_S6_merge.fastq.gz

/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-6s-1120-S1_R1.fastq.gz in2=YW-6s-1120-S1_R2.fastq.gz out=YW_S1_merge.fastq.gz

/home/yuruwang/Tools/anaconda3/bin/bbmerge.sh in1=YW-6s-1120-S2_R1.fastq.gz in2=YW-6s-1120-S2_R2.fastq.gz out=YW_S2_merge.fastq.gz
gunzip YW_S1_merge.fastq.gz
gunzip YW_S2_merge.fastq.gz
gunzip YW_S3_merge.fastq.gz
gunzip YW_S4_merge.fastq.gz
gunzip YW_S5_merge.fastq.gz
gunzip YW_S6_merge.fastq.gz

###need to convert .gz to .fastq
seqkit grep -s -r -p "^TCT" YW_S1_merge.fastq -o ./sorted/Input-sictrlD-rep1.fastq
seqkit grep -s -r -p "^ACA" YW_S1_merge.fastq -o ./sorted/Input-sictrlD-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S1_merge.fastq -o ./sorted/Input-siDKC1-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S1_merge.fastq -o ./sorted/Input-siDKC1-rep2.fastq
seqkit grep -s -r -p "^TCT" YW_S2_merge.fastq -o ./sorted/Input-sictrlT-rep1.fastq
seqkit grep -s -r -p "^ACA" YW_S2_merge.fastq -o ./sorted/Input-sictrlT-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S2_merge.fastq -o ./sorted/Input-siTruB1-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S2_merge.fastq -o ./sorted/Input-siTruB1-rep2.fastq

seqkit grep -s -r -p "^TCT" YW_S3_merge.fastq -o ./sorted/adp-sictrlD-rep1.fastq
seqkit grep -s -r -p "^ACA" YW_S3_merge.fastq -o ./sorted/adp-sictrlD-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S3_merge.fastq -o ./sorted/adp-siDKC1-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S3_merge.fastq -o ./sorted/adp-siDKC1-rep2.fastq
seqkit grep -s -r -p "^TCT" YW_S4_merge.fastq -o ./sorted/adp-sictrlT-rep1.fastq
seqkit grep -s -r -p "^ACA" YW_S4_merge.fastq -o ./sorted/adp-sictrlT-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S4_merge.fastq -o ./sorted/adp-siTruB1-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S4_merge.fastq -o ./sorted/adp-siTruB1-rep2.fastq

seqkit grep -s -r -p "^TCT" YW_S5_merge.fastq -o ./sorted/IP-sictrlD-rep1.fastq
seqkit grep -s -r -p "^ACA" YW_S5_merge.fastq -o ./sorted/IP-sictrlD-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S5_merge.fastq -o ./sorted/IP-siDKC1-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S5_merge.fastq -o ./sorted/IP-siDKC1-rep2.fastq
seqkit grep -s -r -p "^TCT" YW_S6_merge.fastq -o ./sorted/IP-sictrlT-rep1.fastq
seqkit grep -s -r -p "^ACA" YW_S6_merge.fastq -o ./sorted/IP-sictrlT-rep2.fastq
seqkit grep -s -r -p "^CTG" YW_S6_merge.fastq -o ./sorted/IP-siTruB1-rep1.fastq
seqkit grep -s -r -p "^GAC" YW_S6_merge.fastq -o ./sorted/IP-siTruB1-rep2.fastq

