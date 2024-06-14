## To run script:
## python script_name.py SRP123456 /your/path/to/fastq/directory/

import os
import subprocess
import multiprocessing
from multiprocessing import Pool
import pandas as pd
import argparse

def trim_fastq(srr_id, input_dir, output_dir):
    input_file1 = os.path.join(input_dir, f"{srr_id}_1.fastq")
    input_file2 = os.path.join(input_dir, f"{srr_id}_2.fastq")
    input_file_single = os.path.join(input_dir, f"{srr_id}.fastq")

    output_file1_paired = os.path.join(output_dir, f"{srr_id}_1_trimmed_paired.fastq")
    output_file2_paired = os.path.join(output_dir, f"{srr_id}_2_trimmed_paired.fastq")

    output_file1_unpaired = os.path.join(output_dir, f"{srr_id}_1_trimmed_unpaired.fastq")
    output_file2_unpaired = os.path.join(output_dir, f"{srr_id}_2_trimmed_unpaired.fastq")

    output_file1 = os.path.join(output_dir, f"{srr_id}_1.fastq")
    output_file2 = os.path.join(output_dir, f"{srr_id}_2.fastq")
    output_file_single_trimmed = os.path.join(output_dir, f"{srr_id}.fastq")

    # Check if output files are already present
    if os.path.exists(output_file_single_trimmed) or (os.path.exists(output_file1) and os.path.exists(output_file2)):
        print(f"Files for {srr_id} already present. Moving to the next SRR ID.")
        return

    if os.path.exists(input_file1) and os.path.exists(input_file2):
        # Paired-end reads
        trimmomatic_command_paired = f"trimmomatic PE -phred33 {input_file1} {input_file2} {output_file1_paired} {output_file1_unpaired} {output_file2_paired} {output_file2_unpaired} LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36"
        try:
            # Run the trimmomatic command using subprocess
            subprocess.run(trimmomatic_command_paired, shell=True, check=True)
            # Check if paired files exist before renaming and removing unpaired files
            if os.path.exists(output_file1_paired) and os.path.exists(output_file2_paired):
                # Rename the paired trimmed files to final filenames
                os.rename(output_file1_paired, output_file1)
                os.rename(output_file2_paired, output_file2)
                # Remove the unpaired trimmed files
                os.remove(output_file1_unpaired)
                os.remove(output_file2_unpaired)
                print(f"Trimming complete for paired-end reads {srr_id}")
            else:
                print(f"Trimming failed for {srr_id}: Paired files not found.")
        except Exception as e:
            print(f"Error occurred during trimming for {srr_id}: {e}")

    elif os.path.exists(input_file_single):
        # Single-end reads
        trimmomatic_command_single = f"trimmomatic SE -phred33 {input_file_single} {output_file_single_trimmed} LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36"
        try:
            # Run the trimmomatic command using subprocess
            subprocess.run(trimmomatic_command_single, shell=True, check=True)
            print(f"Trimming complete for single-end reads {srr_id}")
        except Exception as e:
            print(f"Error occurred during trimming for {srr_id}: {e}")
    else:
        print(f"No input files found for {srr_id}")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Trim FASTQ files for a given SRP ID.")
    parser.add_argument("srp_id", help="The SRP ID for the project")
    parser.add_argument("base_path", help="The base path for the project directories")
    args = parser.parse_args()

    # Define the directories based on the user input
    input_dir = os.path.join(args.base_path, f'{args.srp_id}_man3')
    output_dir = os.path.join(args.base_path, f'{args.srp_id}_trimmed_man3')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    phenodata_file = os.path.join(args.base_path, f'{args.srp_id}_metadata.csv')

    # Read the srr_ids from the phenodata file
    phenodata_df = pd.read_csv(phenodata_file)
    srr_ids = phenodata_df['Run'].tolist()

    # Set the number of processes to use (adjust according to your CPU cores)
    num_processes = multiprocessing.cpu_count()

    # Create a multiprocessing pool
    with Pool(processes=num_processes) as pool:
        pool.starmap(trim_fastq, [(srr_id, input_dir, output_dir) for srr_id in srr_ids])

if __name__ == "__main__":
    main()
