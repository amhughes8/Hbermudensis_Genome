Mapping mtDNA to assembly before filtering and decontamination to see if the mito genome was assembled into a single contig:
```
/projects/gatins/programs_explorer/minimap2/minimap2 -t 40 -a PV742858.1_HBE_mtdna.fasta ./hifiasm_2.5kQ5/hifiasm_hbe_2.5kQ5.fasta > aln_minimap2.sam
module load samtools/1.21
samtools view -Sb -@ 30 aln_minimap2.sam > mtdna_aligned.bam
samtools sort -@ 20 mtdna_aligned.bam -o mito_aln.sorted.bam
samtools index mito_aln.sorted.bam
samtools view -b -F 4 -@ 20 mito_aln.sorted.bam > mapped.bam
samtools fastq mapped.bam > reads_mito.fastq
```

Can't get MitoZ to work with previous singularity image file... so redownloading it with apptainer ugh
```
apptainer pull MitoZ_v3.6.sif docker://guanliangmeng/mitoz:3.6
```
Okay, it turns out the problem was that I just needed to say 'apptainer run' instead of 'apptainer exec' so I'm going to delete the new version of MitoZ I installed. In the future, can just run from source location
```
apptainer run MitoZ_v3.6.sif mitoz all \
--outprefix hbe_mtdna \
--thread_number 10 \
--clade Chordata \
--genetic_code 2 \
--species_name "Holacanthus bermudensis" \
--fq1 reads_mito.fastq \
--skip_filter \
--assembler megahit \
--kmers_megahit 71 99 \
--memory 50 \
--requiring_taxa Chordata
```
Hmm... I think it can't assemble a mitogenome because it's already been assembled by hifiasm? Let me try just playing around with the contig it assembled as the mitogenome. The size is 33k when I extracted it alone, so I think it must be repeated
```
# to extract:
samtools faidx ./hifiasm_2.5kQ5/hifiasm_hbe_2.5kQ5.fasta ptg000086l > mito_genome.fa
```

It's about twice the size so I think it's just doubled. I'm going to do a rough trim with seqkit
```
seqkit subseq -r 1:17000 mito_genome.fa > mito_genome_half.fa
```

Put this fasta into [mitofish](https://mitofish.aori.u-tokyo.ac.jp/annotation/input)'s annotation program and got mito genome with 13 genes, 2 rRNAs, and 23 tRNAs like we would typically expect. 
![plot](photos/circos.png)
