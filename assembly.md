## Downloading data
```
ssh hughes.annab@xfer.discovery.neu.edu
cd /projects/gatins/2025_HBE_Genome/raw
wget https://s3.amazonaws.com/gtl-public-data/miten/remy/12_1_25_LSK114_AMH_HBE_dorado.1.0.0_sup.bam hughes.annab@xfer.discovery.neu.edu:/projects/gatins/2025_HBE_Genome/raw
```
run time: 03:00:00

## Convert bam to fastq with samtools
```
module load samtools/1.21
samtools bam2fq /projects/gatins/2025_HBE_Genome/raw/12_1_25_LSK114_AMH_HBE_dorado.1.0.0_sup.bam > /projects/gatins/2025_HBE_Genome/assembly/hbe.fastq
```
run time: 01:01:00

## Check fastq stats with SeqKit
```
module load anaconda3/2024.06
source activate /projects/gatins/programs_explorer/SeqKit
seqkit stat hbe.fastq
```
run time: <10 mins

| file   |    format | type   |  num_seqs   |       sum_len | min_len | avg_len   | max_len |
|--------|-----------|--------|-------------|---------------|--------|--------|------|
|hbe.fastq | FASTQ |  DNA |  182,803,626 | 203,672,202,934     |   5 | 1,114.2 | 1,209,147 |

## Trim adapters with Porechop
```
/projects/gatins/programs_explorer/Porechop/porechop-runner.py -i /projects/gatins/2025_HBE_Genome/assembly/hbe.fastq -o hbe_noadapters.fastq --threads 50
```
job id: porechop
run time: 19:34:00

## Filter with SeqKit: min length 10kb, min quality Q5
```
source activate /projects/gatins/programs_explorer/SeqKit/
seqkit seq hbe_noadapters.fastq -m 10000 -Q 5 -j 10 > /projects/gatins/2025_HBE_Genome/assembly/hbe_filtered_10kQ5.fastq
```
### check stats
```
seqkit stat hbe_filtered_10kQ5.fastq
```
| file     |                 format | type | num_seqs    |    sum_len | min_len |  avg_len | max_len |
|--------|-----------|--------|-------------|---------------|--------|--------|------|
| hbe_filtered_10kQ5.fastq | FASTQ |  DNA  |  436,063 | 7,415,328,471  | 10,000 | 17,005.2 | 723,911 |
