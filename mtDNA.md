```
apptainer exec /projects/gatins/2025_HCI_Genome/processing/mtdna/MitoZ_v3.6.sif mitoz all \
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
