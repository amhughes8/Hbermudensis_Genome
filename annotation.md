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
