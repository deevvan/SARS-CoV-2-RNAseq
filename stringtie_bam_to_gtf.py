## To run script:
## stringtie_bam_to_gtf.py path/to/metadata.csv path/to/aligned_bam_files/ paht/to/host_gtf_output/ path/to/covid_gtf_output/


import os
import pandas as pd
from multiprocessing import Pool
import argparse

def run_stringtie_host(args):
    srr_id, main_dir, host_gtf_dir = args
    input_file = os.path.join(main_dir, f"{srr_id}_host_aligned.bam")
    output_file = os.path.join(host_gtf_dir, f"{srr_id}_host.gtf")
    annotation_file_human = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/reference_genome/UCSChg38_human.gtf"
  
    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"File {output_file} already exists.")
        return output_file
    
    os.system(f"stringtie -l {srr_id} {input_file} -o {output_file} -G {annotation_file_human} -e -p 10")
    print(f"StringTie analysis completed for {srr_id}. Output saved to {output_file}.")
    return output_file

def run_stringtie_covid(args):
    srr_id, main_dir, covid_gtf_dir = args
    input_file = os.path.join(main_dir, f"{srr_id}_cov_aligned.bam")
    output_file = os.path.join(covid_gtf_dir, f"{srr_id}_covid.gtf")
    annotation_file_covid = "/mmfs1/projects/changhui.yan/DeewanB/Manuscript3/reference_genome/GCF_009858895.2_ASM985889v3_genomic_ANNOTATED.gff"

    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"{output_file} already exists")
        return

    os.system(f"stringtie -l {srr_id} {input_file} -o {output_file} -G {annotation_file_covid} -e -p 10")
    print(f"StringTie analysis completed for {srr_id}. Output saved to {output_file}.")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Run StringTie analysis on host and SARS-CoV-2 aligned BAM files.")
    parser.add_argument("metadata_file", help="Path to the metadata CSV file")
    parser.add_argument("main_dir", help="The directory containing the aligned BAM files")
    parser.add_argument("host_gtf_dir", help="The directory to save host GTF files")
    parser.add_argument("covid_gtf_dir", help="The directory to save SARS-CoV-2 GTF files")
    args = parser.parse_args()

    # Read the metadata CSV file
    phenodata_df = pd.read_csv(args.metadata_file)

    # Create output directories if they don't exist
    os.makedirs(args.host_gtf_dir, exist_ok=True)
    os.makedirs(args.covid_gtf_dir, exist_ok=True)

    srr_ids = phenodata_df['Run'].tolist()

    num_processes = max(1, os.cpu_count() // 2)

    # Create a list of arguments for the run_stringtie_host function
    host_arguments = [(srr_id, args.main_dir, args.host_gtf_dir) for srr_id in srr_ids]

    # Create a multiprocessing pool for running StringTie on host-aligned BAM files
    with Pool(processes=num_processes) as pool:
        results_host = pool.map(run_stringtie_host, host_arguments)

    # Create a list of arguments for the run_stringtie_covid function
    covid_arguments = [(srr_id, args.main_dir, args.covid_gtf_dir) for srr_id in srr_ids]

    # Create a multiprocessing pool for running StringTie on SARS-CoV-2 aligned BAM files
    with Pool(processes=num_processes) as pool:
        results_covid = pool.map(run_stringtie_covid, covid_arguments)

    print("StringTie analysis completed for all samples.")

if __name__ == "__main__":
    main()
