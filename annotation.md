# 1. RepeatModeler to model repetitive elements in the genome
Build database with final genome assembly

*I had to move the dfam-tetools-latest.sif image over to the folder where the fasta was located or else it would say the fasta file did not exist! Super weird... anyway, now everything is located in **/projects/gatins/2025_HBE_Genome/annotation***

```
# running from within /projects/gatins/2025_HBE_Genome/assembly/hifiasm_2.5kQ5/contam_removal then moving all output and .sif ../../../annotation
apptainer exec dfam-tetools-latest.sif BuildDatabase -name hbe_genome_repeats final_assembly_filtered_nocontam.fasta
```

RepeatModeler:

duration: 28:53:01
```
#!/bin/bash
#SBATCH -J repeatmodeler                    # Job name
#SBATCH -p long                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 50                               # Number of tasks/threads
#SBATCH --mem=50G                           # Memory
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=3-00:00:00                     # Maximum run time
cd /projects/gatins/2025_HBE_Genome/annotation
apptainer exec dfam-tetools-latest.sif RepeatModeler -LTRStruct -database hbe_genome_repeats -threads 50
```
output in /projects/gatins/2025_HBE_Genome/annotation:
```
hbe_genome_repeats-families.fa  #Consensus sequences for each family identified.
hbe_genome_repeats-families.stk  #Seed alignments for each family identified.
hbe_genome_repeats-rmod.log  #Execution log.  Useful for reproducing results.
```

# 2. RepeatMasker

```
/projects/gatins/2025_HCI_Genome/annotation/RepeatMasker/RepeatMasker -pa 10 -lib hbe_genome_repeats-families.fa -xsmall -gff /projects/gatins/2025_HBE_Genome/assembly/hifiasm_2.5kQ5/contam_removal/final_assembly_filtered_nocontam.fasta
```

In order to run this, I had to reconfigure the paths for libraries, TRF, and RMBlast manually using perl ./configure within the RepeatMasker folder since the path changed after the change to Explorer UGH. Just making this note so I remember in case this ever happens again.

```
[hughes.annab@explorer-02 contam_removal]$ cat final_assembly_filtered_nocontam.fasta.tbl
==================================================
file name: final_assembly_filtered_nocontam.fasta
sequences:            30
total length:  604362165 bp  (604362165 bp excl N/X-runs)
GC level:         41.38 %
bases masked:  109328236 bp ( 18.09 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements        89617     17913222 bp    2.96 %
   SINEs:            13848      1466567 bp    0.24 %
   Penelope:          3929       500259 bp    0.08 %
   LINEs:            54447     10558231 bp    1.75 %
    CRE/SLACS            0            0 bp    0.00 %
     L2/CR1/Rex      41242      7560202 bp    1.25 %
     R1/LOA/Jockey    1479       438387 bp    0.07 %
     R2/R4/NeSL        554       171398 bp    0.03 %
     RTE/Bov-B        3185       753669 bp    0.12 %
     L1/CIN4          2292       826686 bp    0.14 %
   LTR elements:     17393      5388165 bp    0.89 %
     BEL/Pao          1971       917011 bp    0.15 %
     Ty1/Copia           0            0 bp    0.00 %
     Gypsy/DIRS1      2712      1373191 bp    0.23 %
       Retroviral     5868      1399904 bp    0.23 %

DNA transposons     187704     26133684 bp    4.32 %
   hobo-Activator    59158      7886350 bp    1.30 %
   Tc1-IS630-Pogo    46885      7806147 bp    1.29 %
   En-Spm                0            0 bp    0.00 %
   MULE-MuDR           472        88325 bp    0.01 %
   PiggyBac           5385       441723 bp    0.07 %
   Tourist/Harbinger  8306      1504431 bp    0.25 %
   Other (Mirage,     2967       415841 bp    0.07 %
    P-element, Transib)

Rolling-circles        481        54686 bp    0.01 %

Unclassified:       416794     50085079 bp    8.29 %

Total interspersed repeats:    94131985 bp   15.58 %


Small RNA:            7545       936231 bp    0.15 %

Satellites:           1112       666215 bp    0.11 %
Simple repeats:     282144     12247454 bp    2.03 %
Low complexity:      35967      2029776 bp    0.34 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element


RepeatMasker version 4.1.9 , default mode
run with rmblastn version 2.14.1+
The query was compared to classified sequences in "hbe_genome_repeats-families.fa"
FamDB:
```

# Gene prediction with BRAKER3

First, we need to map RNAseq data to the reference genome assembly. We are going to use RNAseq from the queen angelfish (_H. ciliaris_) and see how well this works.

```
# index the genome
cd /projects/gatins/2025_HBE_Genome/annotation/hci_rnaseq_mapping
/projects/gatins/programs_explorer/hisat2/hisat2-build -p 20 /projects/gatins/2025_HBE_Genome/assembly/hifiasm_2.5kQ5/contam_removal/final_assembly_filtered_nocontam.fasta.masked HBE_masked

# map
# -x indicates the reference genome index. hisat2 looks for the specified index first in the current directory, then in the directory specified in the HISAT2_INDEXES environment variable.
export HISAT2_INDEXES=/projects/gatins/2025_HBE_Genome/annotation/hci_rnaseq_mapping

# map to genome and create SAM file
cd /projects/gatins/2025_HCI_Genome/rnaseq/fastqs/trimmed
for i in `cat files`; do /projects/gatins/programs_explorer/hisat2/hisat2 -x HBE_masked -1 ${i}_1_polyAremoved_val_1.fq.gz -2 ${i}_2_polyAremoved_val_2.fq.gz -S ${i}_HBE.sam -p 20; done

# load modules
module load samtools/1.21

# convert SAM to BAM
for i in `cat files`; do samtools view -u ${i}_HBE.sam | samtools sort -o ${i}_HBE.bam; done

# merge all sample BAM files
samtools merge -@ 32 hbe_mapped_hci_rnaseq.bam ./*HBE.bam
```

```
mv hbe_mapped_hci_rnaseq.bam /projects/gatins/2025_HCI_Genome/annotation/braker/
```

```
cd /projects/gatins/2025_HBE_Genome/assembly/hifiasm_2.5kQ5/contam_removal/
cp final_assembly_filtered_nocontam.fasta.masked /projects/gatins/2025_HCI_Genome/annotation/braker/

cd /projects/gatins/2025_HCI_Genome/annotation/braker

apptainer exec -B /projects/gatins/2025_HCI_Genome/annotation/braker /projects/gatins/2025_HCI_Genome/annotation/braker/braker3.sif braker.pl \
--genome=/projects/gatins/2025_HCI_Genome/annotation/braker/final_assembly_filtered_nocontam.fasta.masked \
--prot_seq=/projects/gatins/2025_HCI_Genome/annotation/braker/Vertebrata.fa \
--bam=/projects/gatins/2025_HCI_Genome/annotation/braker/hbe_mapped_hci_rnaseq.bam \
--threads=30 --species=Hbermudensis --softmasking \
--AUGUSTUS_CONFIG_PATH=/projects/gatins/2025_HCI_Genome/annotation/braker/config &> hbe_nobusco_hci_rnaseq_braker.log
```

number of protein-coding genes predicted:
```
grep -c "^>" braker.aa
```
26811

Running BUSCO on initial output:
```
#!/bin/bash
#SBATCH -J busco_proteins                    # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 10                               # Number of tasks/threads
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

module load anaconda3/2024.06
source activate /projects/gatins/programs_explorer/busco
busco -i braker.aa --mode proteins --lineage_dataset actinopterygii_odb12 --cpu 10 --out hbe_initial_braker_busco
```
| C:94.6%[S:78.7%,D:15.9%],F:0.8%,M:4.6%,n:7207 |
|-----------------------------------------------|
|	6815 Complete BUSCOs (C) |
|	5672 Complete and single-copy BUSCOs (S) |
|	1143 Complete and duplicated BUSCOs (D) |
|	58 Fragmented BUSCOs (F) |
|	334 Missing BUSCOs (M) |
|	7207 Total BUSCO groups searched |
  
Now, we will use TSEBRA to filter.

First, we will filter out single-exon genes because BRAKER3 likely overestimates them:
running this with 5 cores on an interactive short partition node... took ~3 mins
```
apptainer exec braker3.sif tsebra.py \
-g /projects/gatins/2025_HCI_Genome/annotation/braker/hbe_braker/GeneMark-ETP/genemark.gtf \
-k /projects/gatins/2025_HCI_Genome/annotation/braker/hbe_braker/Augustus/augustus.hints.gtf \
-e /projects/gatins/2025_HCI_Genome/annotation/braker/hbe_braker/hintsfile.gff \
--filter_single_exon_genes \
-c /projects/gatins/2025_HCI_Genome/annotation/braker/default.cfg \
-o hbe_braker_nseg.gtf
```

Next, we will filter for the longest isoform:
```
apptainer exec braker3.sif get_longest_isoform.py --gtf hbe_braker_nseg.gtf --out hbe_braker_nseg_li.gtf
```

extract protein sequences:
```
/projects/gatins/programs_explorer/gffread/bin/gffread -w hbe_braker_nseg_li.fa -y hbe_braker_nseg_li.aa -g /projects/gatins/2025_HCI_Genome/annotation/braker/final_assembly_filtered_nocontam.fasta.masked hbe_braker_nseg_li.gtf
```

BUSCO again:
```
module load anaconda3/2024.06
source activate /projects/gatins/programs_explorer/busco
busco -i hbe_braker_nseg_li.aa --mode proteins --lineage_dataset actinopterygii_odb12 --cpu 10 --out hbe_filtered_braker_busco
```
