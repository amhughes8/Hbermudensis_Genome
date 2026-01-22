Align high quality dataset to ref genome
```
/projects/gatins/programs_explorer/minimap2/minimap2 -t 30 -ax map-ont /projects/gatins/2025_HBE_Genome/assembly/hifiasm_2.5kQ5/contam_removal/filtered_assembly_blobtools.fasta /projects/gatins/2025_HBE_Genome/assembly/hbe_filtered_3kQ10.fastq > HBE_aligned.sam
```

```
module load samtools/1.21
samtools view -Sb -@ 30 -o HBE_aligned.bam HBE_aligned.sam
samtools sort -o HBE_aligned_sorted.bam -O bam -@ 20 HBE_aligned.bam
samtools index -b -@ 20 HBE_aligned_sorted.bam
```

Generate whole-genome consensus file for PSMC and gzip
```
samtools consensus --ambig -f fastq -d 50 HBE_aligned_sorted.bam -o consensus.fq
gzip consensus.fq
```

Run PSMC
```
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/fq2psmcfa -q20 consensus.fq.gz > diploid_HBE.psmcfa
/projects/gatins/2025_HCI_Genome/PSMC/psmc/psmc -N30 -t30 -r5 -p "4+30*2+4+6+10" -o diploid_HBE_final.psmc diploid_HBE.psmcfa
```

Plot
```
source activate /projects/gatins/programs_explorer/gnuplot
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/psmc_plot.pl -u 5.97e-09 -g 5 HBE_t30r5_plot_u597-9g5 diploid_HBE_final.psmc
```
The plot looks so weird? It's just a straight line! I see this has happened to other people before ([see this GitHub issue](https://github.com/lh3/psmc/issues/57)), so going to try changing the r value

```
/projects/gatins/2025_HCI_Genome/PSMC/psmc/psmc -N30 -t30 -r10 -p "4+30*2+4+6+10" -o diploid_HBE_r10.psmc diploid_HBE.psmcfa
source activate /projects/gatins/programs_explorer/gnuplot
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/psmc_plot.pl -u 5.97e-09 -g 5 HBE_t30r10_plot_u597-9g5 diploid_HBE_r10.psmc
```


Hmm.... okay going to run some tests with different r values. Also going to use the depth of 30 consensus file (consensus30.fq.gz)
```
# 10
/projects/gatins/2025_HCI_Genome/PSMC/psmc/psmc -N30 -t30 -r10 -p "4+30*2+4+6+10" -o diploid_HBE30_r10.psmc diploid_HBE30.psmcfa

# 2
/projects/gatins/2025_HCI_Genome/PSMC/psmc/psmc -N30 -t30 -r10 -p "4+30*2+4+6+10" -o diploid_HBE30_r2.psmc diploid_HBE30.psmcfa

#7
/projects/gatins/2025_HCI_Genome/PSMC/psmc/psmc -N30 -t30 -r10 -p "4+30*2+4+6+10" -o diploid_HBE30_r7.psmc diploid_HBE30.psmcfa

# plot
source activate /projects/gatins/programs_explorer/gnuplot
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/psmc_plot.pl -u 5.97e-09 -g 5 HBE30_t30r10_plot_u597-9g5 diploid_HBE30_r10.psmc
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/psmc_plot.pl -u 5.97e-09 -g 5 HBE30_t30r2_plot_u597-9g5 diploid_HBE30_r2.psmc
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/psmc_plot.pl -u 5.97e-09 -g 5 HBE30_t30r7_plot_u597-9g5 diploid_HBE30_r7.psmc
```

Bootstrapping
```
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/splitfa diploid_HBE.psmcfa > diploid_HBE_split.psmcfa
```
```
mkdir bootstrap
cp diploid_HBE_split.psmcfa bootstrap
cp diploid_HBE_final.psmc bootstrap
cd bootstrap
echo split_HBE_{001..100}.psmcfa| xargs -n 1 cp diploid_HBE_split.psmcfa
```

```
#!/bin/bash
#SBATCH -J psmc_array			    # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 2                                # Number of tasks/threads
#SBATCH -o array_%A_%a.out    		    # Name of stdout output file
#SBATCH -e array_%A_%a.err    		    # Name of stdout output file
#SBATCH --array=1-100			    # Array index
#SBATCH --mem=6000MB 			    # Memory to be allocated PER NODE
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#
# ----------------Your Commands------------------- #
#
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# select our filename
N=${SLURM_ARRAY_TASK_ID}
# Comment one of the following two lines, depending on if the file names have leading zeros
#FILENAME=run-${N}.inp # without leading zeros
 FILENAME=split_HCI_$(printf "%03d" ${N}).psmcfa # with leading zeros
# adjust "%03d" to as many digits as are in the numeric part of the file name
echo "My input file is ${FILENAME}"

#
echo $P
#
/projects/gatins/2025_HCI_Genome/PSMC/psmc/psmc -N30 -t30 -r5 -b -p "4+30*2+4+6+10" -o /projects/gatins/2025_HBE_Genome/PSMC/bootstrap/${FILENAME}.psmc /projects/gatins/2025_HBE_Genome/PSMC/bootstrap/${FILENAME}
#

echo "Job finished" `date`
echo "My input file is ${FILENAME}"
```

```
cat *.psmc > HBE_combined.psmc
```

Bootstrapping new params
```
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/splitfa diploid_HBE30.psmcfa > diploid_HBE30_split.psmcfa
```
```
mkdir bootstrap30
cp diploid_HBE30_split.psmcfa bootstrap30
cp diploid_HBE30_final.psmc bootstrap30
cd bootstrap30
echo split_HBE30_{001..100}.psmcfa| xargs -n 1 cp diploid_HBE30_split.psmcfa
```

```
#!/bin/bash
#SBATCH -J psmc_array			    # Job name
#SBATCH -p short                            # Partition
#SBATCH -N 1                                # Number of nodes
#SBATCH -n 2                                # Number of tasks/threads
#SBATCH -o array_%A_%a.out    		    # Name of stdout output file
#SBATCH -e array_%A_%a.err    		    # Name of stdout output file
#SBATCH --array=1-100			    # Array index
#SBATCH --mem=6000MB 			    # Memory to be allocated PER NODE
#SBATCH --mail-user=hughes.annab@northeastern.edu  # Email
#SBATCH --mail-type=END                     # Email notification at job completion
#SBATCH --time=48:00:00                     # Maximum run time

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#
# ----------------Your Commands------------------- #
#
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# select our filename
N=${SLURM_ARRAY_TASK_ID}
# Comment one of the following two lines, depending on if the file names have leading zeros
#FILENAME=run-${N}.inp # without leading zeros
FILENAME=split_HBE30_$(printf "%03d" ${N}).psmcfa # with leading zeros
# adjust "%03d" to as many digits as are in the numeric part of the file name
echo "My input file is ${FILENAME}"

#
echo $P
#
/projects/gatins/2025_HCI_Genome/PSMC/psmc/psmc -N30 -t30 -r10 -b -p "1+1+1+1+30*2+4+6+10" -o /projects/gatins/2025_HBE_Genome/PSMC/bootstrap30/${FILENAME}.psmc /projects/gatins/2025_HBE_Genome/PSMC/bootstrap30/${FILENAME}
#

echo "Job finished" `date`
echo "My input file is ${FILENAME}"
```

```
cat *.psmc > HBE30_combined.psmc
```

```
module load anaconda3/
source activate /projects/gatins/programs_explorer/gnuplot
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/psmc_plot.pl -u 5.97e-09 -g 5 HBE30_t30r10_plot_u597-9g5 HBE_combined.psmc
```
