## Downloading data
```
ssh hughes.annab@xfer.discovery.neu.edu
cd /projects/gatins/2025_HBE_Genome/raw
wget https://s3.amazonaws.com/gtl-public-data/miten/remy/12_1_25_LSK114_AMH_HBE_dorado.1.0.0_sup.bam hughes.annab@xfer.discovery.neu.edu:/projects/gatins/2025_HBE_Genome/raw
```
run time: started at 13:29

## Convert bam to fastq with samtools
```
module load samtools/1.9
samtools bam2fq /projects/gatins/2025_HBE_Genome/raw/12_1_25_LSK114_AMH_HBE_dorado.1.0.0_sup.bam > /projects/gatins/2025_HBE_Genome/assembly/hbe.fastq
```
run time:

## Check fastq stats with SeqKit
```
module load anaconda3/2024.06
seqkit stat /projects/gatins/2025_HBE_Genome/assembly/*.fastq
```
run time:

## Trim adapters with Porechop
```
/projects/gatins/programs_explorer/Porechop/porechop-runner.py -i /projects/gatins/2025_HBE_Genome/assembly/hbe.fastq -o hbe_noadapters.fastq
```
job id: porechop
run time:
