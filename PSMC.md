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

source activate /projects/gatins/programs_explorer/gnuplot
/projects/gatins/2025_HCI_Genome/PSMC/psmc/utils/psmc_plot.pl -u 5.97e-09 -g 5 HBE_t30r5_plot_u597-9g5 diploid_HBE_final.psmc
