# SARS-CoV-2 RNAseq scripts
This README file describes the various python and R scripts used in order to process the human/mouse RNAseq datasets with respiratory virus infections such as SARS-CoV-2, IAV and RSV.

## Table of Content:
1. fasterqdump_downloader.py: Python script used to download all SRR submissions (.fastq files) from NCBI for user entered SRP-ID
2. trimmomatic_qc.py: Python script to trim low quality reads from downloaded .fastq files
3. hisat2_align_extract.py: Python script to align high-quality mRNA reads from .fastq files to host genome using hisat2 and hisat2 host genome indices, followed by extraction of unaligned mRNA reads to be converted to .fastq files and be further aligned to SARS-CoV-2 hisat2 genome indices.
4. bowtie2_align_extract.py: Python script to align high-quality mRNA reads from .fastq files to host genome using bowtie2 and bowtie2 host genome indices, followed by extraction of unaligned mRNA reads to be converted to .fastq files for kraken2 metagenomic analysis downstream.
5. stringtie_bam_to_gtf.py: Python script to align reads from .bam file obtained from step 4 to respective gene annotation (.gtf) files to obtain transcript and gene counts using stringTie
6. stringtie_native_preDE.py: Python script by stringTie developers that allows output .gtf files from step 5 into raw read count matrix for downstream DE analysis.
7. deseq2_man2_human.R: R script that performs DE and downstream analyses and WGCNA on read count matrix obtained from step 6.
