## To run script:
## hisat2_align_extract.py path/to/metadata.csv path/to/trimmed_fastq_directory path/to/host_genome_index path/to/covid19_genome_index path/to/output_directory


import os
import pandas as pd
import subprocess
from multiprocessing import Pool
import argparse

def hisat2_alignment(srr_id, main_dir, host_hisat2_index, cov_idx, out_dir):
    # Set output prefixes for host genome
    host_output_prefix = os.path.join(out_dir, srr_id)
    
    # Check if the cov_aligned BAM file already exists
    cov_aligned_bam = f"{host_output_prefix}_cov_aligned.bam"
    if os.path.exists(cov_aligned_bam):
        print(f"{srr_id} bam file already exists. Skipping the alignment step.")
        return
    
    rna_seq_reads1 = os.path.join(main_dir, f"{srr_id}_1.fastq")
    rna_seq_reads2 = os.path.join(main_dir, f"{srr_id}_2.fastq")

    try:
        # Aligning reads to host genome index
        if os.path.exists(rna_seq_reads1) and os.path.exists(rna_seq_reads2):
            # Paired-end reads
            os.system(f"hisat2 -p 8 --dta -x {host_hisat2_index} -1 {rna_seq_reads1} -2 {rna_seq_reads2} -S {host_output_prefix}_host.sam")
        elif os.path.exists(rna_seq_reads1):
            # Single-end reads
            os.system(f"hisat2 -p 8 --dta -x {host_hisat2_index} -U {rna_seq_reads1} -S {host_output_prefix}_host.sam")
        else:
            print(f"No input files found for {srr_id}")
            return

        # Sort the sam file and convert it to bam
        os.system(f"samtools sort -@ 8 {host_output_prefix}_host.sam -o {host_output_prefix}_host.bam")

        # Filter the bam file to contain only mapped reads
        os.system(f"samtools view -b -F 4 {host_output_prefix}_host.bam > {host_output_prefix}_host_aligned.bam")
        
        # Extracting unmapped reads from {host_output_prefix}_host.sam
        os.system(f"samtools view -b -f 4 {host_output_prefix}_host.bam > {host_output_prefix}_unaligned.bam")
        
        # Converting unmapped reads back to fastq using samtools
        os.system(f"samtools fastq -0 /dev/null -s /dev/null -n {host_output_prefix}_unaligned.bam -1 {host_output_prefix}_unaligned.1.fastq -2 {host_output_prefix}_unaligned.2.fastq")
        
        # Aligning unmapped reads to covid19 genome index
        os.system(f"hisat2 -p 8 --dta -x {cov_idx} -1 {host_output_prefix}_unaligned.1.fastq -2 {host_output_prefix}_unaligned.2.fastq -S {host_output_prefix}_cov_aligned.sam")
        
        # Sorting and converting {host_output_prefix}_cov_aligned.sam to bam file
        os.system(f"samtools view -F 4 -b {host_output_prefix}_cov_aligned.sam | samtools sort -@ 8 -o {host_output_prefix}_cov_aligned.bam")
        
        # Removing intermediates
        os.remove(f"{host_output_prefix}_unaligned.1.fastq")
        os.remove(f"{host_output_prefix}_unaligned.2.fastq")
        
        os.remove(f"{host_output_prefix}_host.sam")
        os.remove(f"{host_output_prefix}_cov_aligned.sam")
        os.remove(f"{host_output_prefix}_unaligned.bam")
        
    except Exception as e:
        print(f"Error: {e}")
        return
    
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process SRR IDs for a given SRP ID.")
    parser.add_argument("csv_file", help="The path to the metadata CSV file")
    parser.add_argument("main_dir", help="The directory containing the trimmed FASTQ files")
    parser.add_argument("host_index", help="The directory of the host genome index")
    parser.add_argument("cov_index", help="The directory of the COVID-19 genome index")
    parser.add_argument("out_dir", help="The directory to save the output BAM files")
    args = parser.parse_args()

    # Read CSV file to get SRR IDs
    phenodata_df = pd.read_csv(args.csv_file)
    srr_ids = phenodata_df['Run'].tolist()

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    num_threads = 3

    # Create a pool of threads
    with Pool(processes=num_threads) as pool:
        pool.starmap(hisat2_alignment, [(srr_id, args.main_dir, args.host_index, args.cov_index, args.out_dir) for srr_id in srr_ids])

    print("All alignment processes completed.")

if __name__ == "__main__":
    main()
