## To run the script : 
## fasterqdump_downloader.py SRP123456 /your/path/to/fastq/directory/


import os
import multiprocessing
import pandas as pd
import argparse

def download_sra_run(srr_id, output_dir, check_dir):
    # Check if the file already exists
    fastq_file_path = os.path.join(check_dir, f"{srr_id}_1.fastq")
    
    if not os.path.exists(fastq_file_path):
        # If the file doesn't exist, download it
        os.system(f"fasterq-dump -O {output_dir} {srr_id}")
        print(f"{srr_id} downloaded.")
    else:
        print(f"{srr_id} already exists. Skipping download.")

def download_srr_files(output_dir, phenodata_file, check_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read the phenodata file
    phenodata_df = pd.read_csv(phenodata_file)

    # Get the list of SRR IDs from phenodata_df
    srr_ids = phenodata_df['Run'].tolist()

    # Set the number of processes to use (adjust according to your CPU cores)
    num_processes = multiprocessing.cpu_count()

    # Create a multiprocessing pool
    pool = multiprocessing.Pool(processes=num_processes)

    # Map the download function to the SRR IDs using the multiprocessing pool
    pool.starmap(download_sra_run, [(srr_id, output_dir, check_dir) for srr_id in srr_ids])

    # Close the pool
    pool.close()
    pool.join()

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Download SRR files for a given SRP ID.")
    parser.add_argument("srp_id", help="The SRP ID for the project")
    parser.add_argument("base_path", help="The base path for the project directories")
    args = parser.parse_args()

    # Define the directories based on the user input
    output_dir = os.path.join(args.base_path, f'{args.srp_id}_fastq')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    check_dir = os.path.join(args.base_path, f'{args.srp_id}_trimmed_fastq')
    if not os.path.exists(check_dir):
        os.makedirs(check_dir)

    phenodata_file = os.path.join(args.base_path, f'{args.srp_id}_metadata.csv')

    # Call the function to download the SRR files
    download_srr_files(output_dir, phenodata_file, check_dir)

if __name__ == "__main__":
    main()
