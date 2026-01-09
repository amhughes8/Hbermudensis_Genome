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
| file     |    format | type | num_seqs    |    sum_len | min_len |  avg_len | max_len |
|--------|-----------|--------|-------------|---------------|--------|--------|------|
| hbe_filtered_10kQ5.fastq | FASTQ |  DNA  |  436,063 | 7,415,328,471  | 10,000 | 17,005.2 | 723,911 |

## Check coverage with jellyfish
```
cd /projects/gatins/programs_explorer/jellyfish_2.2/bin
./jellyfish count -m 21 -s 500M -t 10 -C -o /projects/gatins/2025_HBE_Genome/assembly/hbeQ5_21mer_output /projects/gatins/2025_HBE_Genome/assembly/hbe_filtered_10kQ5.fastq
```
Generate histogram:
```
./jellyfish histo /projects/gatins/2025_HBE_Genome/assembly/hbeQ5_21mer_output > /projects/gatins/2025_HBE_Genome/assembly/hbeQ5_21mer_output.histo
```

The jellyfish output showed that this subset of data only gives about 10x coverage. Let's try filtering to a 5k cutoff and see if this improves.

## Use SeqKit to filter for 5k min length
| file     |    format | type | num_seqs    |    sum_len | min_len |  avg_len | max_len |
|--------|-----------|--------|-------------|---------------|--------|--------|------|
| hbe_filtered_5kQ5.fastq | FASTQ  | DNA  | 3,082,960 | 24,622,611,268 |   5,000 | 7,986.7 | 723,911 |

## check coverage on new dataset with jellyfish
```
cd /projects/gatins/programs_explorer/jellyfish_2.2/bin
./jellyfish count -m 21 -s 500M -t 10 -C -o /projects/gatins/2025_HBE_Genome/assembly/hbe5kQ5_21mer_output /projects/gatins/2025_HBE_Genome/assembly/hbe_filtered_5kQ5.fastq
```
Generate histogram:
```
./jellyfish histo /projects/gatins/2025_HBE_Genome/assembly/hbe5kQ5_21mer_output > /projects/gatins/2025_HBE_Genome/assembly/hbe5kQ5_21mer_output.histo
```
~35x coverage, moving forward

## Remove mtDNA
download mito genome from NCBI (accession: PV742858.1)
```
PV742858.1_HBE_mtdna.fasta
```
Use minimap2 to map filtered fastq reads pre-assembly (hbe_filtered_5kQ5.fastq) to reference mitochondrial sequence (PV742858.1_HBE_mtdna.fasta)
```
/projects/gatins/programs_explorer/minimap2/minimap2 -t 40 -ax map-ont PV742858.1_HBE_mtdna.fasta hbe_filtered_5kQ5.fastq > aln_minimap2.sam
```
Use samtools to extract sequences that mapped to mitochondrial genome
```
module load samtools/1.21
samtools view -Sb -@ 30 aln_minimap2.sam > mtdna_aligned.bam
samtools sort -@ 20 mtdna_aligned.bam -o mito_aln.sorted.bam
samtools index mito_aln.sorted.bam
```
Extract sequences that were unmapped and save them to a new BAM file
```
samtools view -b -f 4 -@ 20 mito_aln.sorted.bam > unmapped.bam
```
Convert BAM to FASTQ
```
samtools fastq unmapped.bam > reads_no_mito.fastq
```

## Assemble with hifiasm
```
/projects/gatins/programs_explorer/hifiasm/hifiasm -o hifiasm_hbe.asm --ont -t32 /projects/gatins/2025_H
BE_Genome/assembly/reads_no_mito.fastq
```
job id: hifiasm_hbe

run time: 4 hrs

## BUSCO
```
module load anaconda3/2022.05
source activate /projects/gatins/programs_explorer/busco
busco -i /projects/gatins/2025_HBE_Genome/jobs/hbe_hifiasm_assembly.fasta --mode genome --lineage_dataset actinopterygii_odb12 --cpu 25 --out hbe_hifiasm_nomito_busco
```

## SeqKit stats on assembly + BUSCO
| file    |        format | type | num_seqs  |  sum_len | min_len |   avg_len  |  max_len |  Q1  | Q2 | Q3 | sum_gap  |  N50 | N50_num | Q20(%) | Q30(%) | AvgQual | GC(%) | sum_n | BUSCO |
|---------|---------------|-------|---------|------------|--------|------------|-----------|-----|-----|----|---------|-------|--------|---------|-------|---------|-------|-------|--------|
|hbe_hifiasm_assembly.fasta | FASTA  | DNA     |   699 | 623,092,146  |  6,162 | 891,405.1 | 8,801,329 | 34,613 | 234,787 | 1,309,041.5  |  0 | 2,363,646   |    80     |  0 |   0  | 0 | 41.38 | 0 | 98.8% |

Despite a very complete assembly, there are a ton of contigs so I'm confused. I used FCS to check for adapter contamination and it looks like Porechop did not completely remove all adapter content, so this could possibly contribute to the issue? Let's try to use Porechop_ABI (ab initio) -- an extension of Porechop.

FCS command:
```
cd /projects/gatins/2025_HCI_Genome/processing/fcs
./run_fcsadaptor.sh --fasta-input /projects/gatins/2025_HBE_Genome/assembly/hifiasm_5kQ5_nomito/hbe_assembly_15kcutoff.fasta --output-dir /projects/gatins/2025_HBE_Genome/assembly/hifiasm_5kQ5_nomito --euk --container-engine singularity --image fcs-adaptor.sif
```
FCS output:
| accession	| length	| action |	range	name |
|-----------|---------|-------|----------|
| ptg000037l |	1956817	| ACTION_TRIM	| 18191..18218,28770..28796,41816..41843,82031..82056,97905..97932	| CONTAMINATION_SOURCE_TYPE_ADAPTOR:NGB02000.1:Oxford Nanopore Technologies Rapid Adapter (RA) Ligation Adapter top (LA) Native Adaptor top (NA) polyT masked |
| ptg000372l	| 65086	| ACTION_TRIM	| 6471..6495	| CONTAMINATION_SOURCE_TYPE_ADAPTOR:NGB02001.1:Oxford Nanopore Technologies Ligation Adapter bottom (LA) |

## Porechop_ABI
running with 100G RAM and 10 threads... we'll see how it goes
```
porechop_abi --ab_initio -i /projects/gatins/2025_HBE_Genome/assembly/hbe.fastq -o /projects/gatins/2025_HBE_Genome/assembly/hbe_noadapters_abi.fastq
```

## Back to filtering, going to follow queen angel protocol exactly
Filtering to 2.5kb min length

| file |  format | type  |  num_seqs | sum_len | min_len | avg_len | max_len  |  Q1  |   Q2   |  Q3 | sum_gap |   N50 | N50_num | Q20(%) | Q30(%) |  AvgQual | GC(%) | sum_n |
|-----|----------|------|----------|--------|----------|---------|---------|-------|------|------|-------|-------|------|-------|-------|------|------|------|
| hbe_filtered2.5kQ5.fastq | FASTQ  | DNA |  15,553,224 | 66,328,006,101  |  2,500 | 4,264.6 | 723,911 | 2,871 | 3,445 | 4,574  | 0 | 4,153 | 49,208  | 87.88 |  81.93  |  17.12 | 41.61  |    0 |

## Re-assembling
First, I'll try without removing the mtDNA but then based on how the assembly looks, ill remove it
