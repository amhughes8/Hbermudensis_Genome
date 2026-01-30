# 1. RepeatModeler to model repetitive elements in the genome
Build database with final genome assembly

*I had to move the dfam-tetools-latest.sif image over to the folder where the fasta was located or else it would say the fasta file did not exist! Super weird... anyway, now everything is located in **/projects/gatins/2025_HBE_Genome/annotation***

```
# running from within /projects/gatins/2025_HBE_Genome/assembly/hifiasm_2.5kQ5/contam_removal then moving all output and .sif ../../../annotation
apptainer exec dfam-tetools-latest.sif BuildDatabase -name hbe_genome_repeats final_assembly_filtered_nocontam.fasta
```

RepeatModeler:

started running at 2:32pm on Jan 28th
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
