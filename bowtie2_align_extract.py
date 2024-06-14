## To run script:
## bowtie2_align_extract.py SRP123456 path/to/main/directory path/to/reference/genome/index/{genome_idx}


import os
import pandas as pd
import subprocess
from multiprocessing import Pool
import argparse

def process_srr_id(srr_id, input_dir, reference_genome_dir, output_dir):
    # Define input and output paths
    input_fastq1 = os.path.join(input_dir, f"{srr_id}_1.fastq")
    input_fastq2 = os.path.join(input_dir, f"{srr_id}_2.fastq")
    output_unmapped = f"{output_dir}/{srr_id}_non_host.fa"
    output_mapped = f"{output_dir}/{srr_id}_host.sam"
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Check if input files are single-ended or paired-ended
    if os.path.exists(input_fastq1) and os.path.exists(input_fastq2):
        # Paired-end reads
        cmd = [
            "bowtie2",
            "-x", reference_genome_dir,
            "-p", "8",
            "-q",
            "-1", input_fastq1,
            "-2", input_fastq2,
            "--un-conc", output_unmapped,
            "-S", output_mapped
        ]
    elif os.path.exists(input_fastq1):
        # Single-end reads
        cmd = [
            "bowtie2",
            "-x", reference_genome_dir,
            "-p", "8",
            "-q",
            "-U", input_fastq1,
            "--un", output_unmapped,
            "-S", output_mapped
        ]
    else:
        print(f"No input files found for {srr_id}")
        return

    try:
        subprocess.run(cmd, check=True)
        print(f"Processed SRR ID: {srr_id}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {srr_id}: {e}")

def get_srr_ids(csv_file):
    df = pd.read_csv(csv_file)
    return df['Run'].tolist()

def process_srr_ids(csv_file, input_dir, reference_genome_dir, output_dir):
    srr_ids = get_srr_ids(csv_file)
    
    # Use multiprocessing to parallelize the processing
    with Pool(processes=8) as pool:
        pool.starmap(process_srr_id, [(srr_id, input_dir, reference_genome_dir, output_dir) for srr_id in srr_ids])

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process SRR IDs for a given SRP ID.")
    parser.add_argument("srp_id", help="The SRP ID for the project")
    parser.add_argument("base_path", help="The base path for the project directories")
    parser.add_argument("reference_genome_dir", help="The directory of the reference genome")
    args = parser.parse_args()

    # Define paths based on the user input
    csv_file = os.path.join(args.base_path, f'{args.srp_id}_metadata.csv')
    input_dir = os.path.join(args.base_path, 'SRR_trimmed_fastq')
    output_dir = os.path.join(args.base_path, 'bowtie2_output')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Process SRR IDs
    process_srr_ids(csv_file, input_dir, args.reference_genome_dir, output_dir)

if __name__ == "__main__":
    main()
