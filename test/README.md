# test data set

Here you will find a small test data set to ensure that kallisto and pizzly are working properly. The contents are:

- `transcripts.fasta.gz`: a gzip compressed transcriptome
- `transcripts.gtf.gz` : a small subset of the GTF file for  Ensembl version 85 corresponding to the transcripts provided
- `reads_X.fastq.gz`: gzip compressed "left" and "right" reads
- `Snakefile`: a sample [`snakemake`](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) file (not required for running kallisto or pizzly)

Running `snakemake` will go through the following pipeline

1. Creates the index for use with `kallisto`
2. Runs `kallisto` with `--fusion` to identify potential fusions
3. Runs `pizzly` on the `kallisto` output to identify potential fusions
4. Creates a new index based on the transcriptome and the fusion transcripts identified by `pizzly`
5. Runs `kallisto` in normal quantification mode on the expanded index to quantify both normal transcripts and fusions.
