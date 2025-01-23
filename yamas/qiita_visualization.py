from __future__ import annotations

import pandas as pd
import csv
import gzip
import os.path
import pickle
import datetime
from tqdm import tqdm
from metaphlan.utils.merge_metaphlan_tables import merge

from .utilities import run_cmd, ReadsData, check_conda_qiime2
import json
import shutil
import tarfile
import os
import yaml


def combine_tsv_files(file_list, output_file):
    """
    Combines multiple TSV files into one, keeping only the first header.

    Parameters:
    - file_list: List of paths to the TSV files.
    - output_file: Path for the combined output TSV file.

    Returns:
    - None: Writes the combined data to the output file.
    """
    # Initialize an empty list to store DataFrames
    dataframes = []

    for idx, file in enumerate(file_list):
        # Read the file into a DataFrame
        if idx == 0:
            # Include headers for the first file
            df = pd.read_csv(file, sep="\t")
        else:
            # Skip headers for subsequent files
            df = pd.read_csv(file, sep="\t", header=0)
        dataframes.append(df)

    # Combine all DataFrames
    combined_df = pd.concat(dataframes, ignore_index=True)

    # Save the combined DataFrame to a new TSV file
    combined_df.to_csv(output_file, sep="\t", index=False)
    return output_file


def handle_metadata_duplicates(dir_path):
    """
    Checks for duplicate barcodes in metadata.tsv, creates a duplicate_metadata.tsv file,
    and removes all rows with duplicate barcodes from metadata.tsv.

    Parameters:
    dir_path (str): Path to the directory containing metadata.tsv.
    """
    metadata_file = os.path.join(dir_path, "metadata.tsv")
    duplicate_file = os.path.join(dir_path, "duplicate_metadata.tsv")

    # Check if metadata file exists
    if not os.path.exists(metadata_file):
        print(f"No metadata.tsv file found in {dir_path}. Skipping...")
        return

    # Read the metadata file
    metadata = pd.read_csv(metadata_file, sep="\t")

    # Check if the barcode column exists
    if 'barcode' not in metadata.columns:
        print(f"No 'barcode' column found in {metadata_file}. Skipping...")
        return

    # Identify all barcodes that have duplicates
    duplicate_barcodes = metadata['barcode'][metadata['barcode'].duplicated(keep=False)]

    # Separate duplicates and unique rows
    duplicates = metadata[metadata['barcode'].isin(duplicate_barcodes)]
    unique_metadata = metadata[~metadata['barcode'].isin(duplicate_barcodes)]

    # Save duplicates to duplicate_metadata.tsv
    if not duplicates.empty:
        duplicates.to_csv(duplicate_file, sep="\t", index=False)

    # Save the cleaned metadata.tsv
    unique_metadata.to_csv(metadata_file, sep="\t", index=False)

def generate_manifest(input_directory, output_manifest):
    """
    Generate a manifest file for paired-end sequencing data.

    Parameters:
    input_directory (str): Directory containing the input FASTQ files.
    output_manifest (str): Path to the output manifest file.

    Returns:
    str: Path to the generated manifest file.
    """
    manifest_path = os.path.join(input_directory, output_manifest)
    run_cmd([f"touch {manifest_path}"])
    with open(manifest_path, 'w') as manifest:
        # Write the header line
        manifest.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")

        # Iterate over all files matching the forward read pattern
        for forward in sorted(os.listdir(input_directory)):
            if '_R1_' in forward and forward.endswith('.fastq.gz'):
                sample_id = forward.split('_R1_')[0]  # Extract sample ID
                reverse = f"{sample_id}_R2_001.fastq.gz"  # Construct reverse read filename

                forward_path = os.path.join(input_directory, forward)
                reverse_path = os.path.join(input_directory, reverse)

                if os.path.exists(reverse_path):
                    # Add sample ID, forward, and reverse paths to the manifest
                    manifest.write(f"{sample_id}\t{forward_path}\t{reverse_path}\n")
                else:
                    print(f"Warning: Reverse read file not found for {sample_id}")

    return manifest_path


def merge_fastq(main_dir):
    """
    This function merges the fastq fo big project analysis such as the AGP for better user experience.
    :param main_dir: The path to the main directory
    :return: creates a new dir within the main dir that contains the merged items
    """

    os.makedirs(os.path.join(main_dir, "all", "fastq"), exist_ok=True)
    os.makedirs(os.path.join(main_dir, "all", "metadata"), exist_ok=True)
    merged_forward = os.path.join(main_dir, "all", "fastq", "forward.fastq.gz")
    merged_reverse = os.path.join(main_dir, "all", "fastq", "reverse.fastq.gz")
    merged_barcodes = os.path.join(main_dir, "all", "fastq", "barcodes.fastq.gz")
    merged_metadata =  os.path.join(main_dir, "all", "metadata", "metadata.tsv")
    metadata = []
    for directory in os.listdir(main_dir):
        if directory == "all" or not os.path.exists(os.path.join(main_dir, directory, "fastq")):
            continue
        metadata_path = os.path.join(main_dir, directory, "metadata")
        handle_metadata_duplicates(metadata_path)
        metadata.append(os.path.join(metadata_path, "metadata.tsv"))
        fastq_path = os.path.join(main_dir, directory, "fastq")
        forward = os.path.join(fastq_path, "forward.fastq.gz")
        reverse = os.path.join(fastq_path, "reverse.fastq.gz")
        barcodes = os.path.join(fastq_path, "barcodes.fastq.gz")
        if  len(os.listdir(fastq_path)) == 3:
            command = [f"cat {forward} >> {merged_forward} && cat {reverse} >> {merged_reverse} && cat {barcodes} >> {merged_barcodes}"]
            run_cmd(command)

    combine_tsv_files(metadata, merged_metadata)
    return os.path.join(main_dir, "all", "fastq")


def check_paired_flag(fastq_path):
    # Define the expected filenames
    forward_file = os.path.join(fastq_path, "forward.fastq.gz")
    reverse_file = os.path.join(fastq_path, "reverse.fastq.gz")
    barcodes_file = os.path.join(fastq_path, "barcodes.fastq.gz")

    # Check for file presence and set the paired flag
    if os.path.exists(forward_file) and os.path.exists(reverse_file) and os.path.exists(barcodes_file):
        return True
    elif os.path.exists(forward_file):
        return False
    else:
        print("Please check your input, you do not have the files with the names as asked or you entered data that we currently can not process.")


def metaphlan_extraction(reads_data):
    paired = reads_data.rev and reads_data.fwd
    fastq_path = os.path.join(reads_data.dir_path, "fastq")
    export_path = os.path.join(reads_data.dir_path, "export")
    run_cmd([f"mkdir {export_path}"])
    final_output_path = os.path.join(export_path, 'final.txt')
    run_cmd([f"touch {final_output_path}"])
    fastq_files = [a for a in os.listdir(fastq_path) if a.split(".")[-1] == "fastq"]
    args = []
    if paired:
        print("paired")
        for i in tqdm(range(0,len(fastq_files),2)):
            fastq_name = fastq_files[i].split('_')[0]
            fastq_1 = os.path.join(fastq_path, fastq_files[i])
            fastq_2 = os.path.join(fastq_path, fastq_files[i+1])
            output = os.path.join(fastq_path,f"{fastq_name}.bowtie2.bz2")
            command = [f"metaphlan {fastq_1},{fastq_2} --input_type fastq --bowtie2out {output} --nproc 24"]
            run_cmd(command)
            final_output_file = os.path.join(os.path.join(reads_data.dir_path, 'qza'), f'{fastq_name}_profile.txt')
            command = [f"metaphlan {output} --input_type bowtie2out --nproc 24 > {final_output_file}"]
            run_cmd(command)
            args.append(final_output_file)
    else:
        print("not paired")
        for fastq in tqdm(fastq_files):
            output = os.path.join(os.path.join(reads_data.dir_path, 'qza'), f'{fastq}_profile.txt')
            fastq = os.path.join(fastq_path,fastq)
            command = [f"metaphlan {fastq} --input_type fastq --nproc 24 > {output}"]
            run_cmd(command)
            args.append(output)
    merge(args, open(final_output_path, 'w'), gtdb=False)


# noinspection PyTypeChecker
def metaphlan_txt_csv(reads_data, dataset_id):
    export_path = os.path.join(reads_data.dir_path, "export")
    input_file = os.path.join(export_path,f"{dataset_id}_final.txt")
    output_file = os.path.join(export_path,f"{dataset_id}_final_table.csv")
    with open(input_file, 'r') as txt_file:
        lines = txt_file.readlines()

    # Extract data from the text file
    headers = lines[0].strip().split('\t')
    data = [line.strip().split('\t') for line in lines[1:]]

    # Replace "|" with ","
    headers = [header.replace('|', ',') for header in headers]
    data = [[entry.replace('|', ',') for entry in row] for row in data]

    # Transpose the data
    transposed_data = list(map(list, zip(*data)))

    # Write the transposed data to a CSV file
    with open(output_file, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(transposed_data)

    print(f"CSV file '{output_file}' has been created.")

def qiime_import(dir_path,fastq_path, paired):
    qza_path = os.path.join(dir_path, "qza")

    multiplexed_qza_file_path = os.path.join(qza_path, "multiplexed-seqs.qza")

    command = [
        "qiime", "tools", "import",
        "--type", f"{'EMPPairedEndSequences' if paired else 'MultiplexedSingleEndBarcodeInSequence'}",
        "--input-path", fastq_path,
        "--output-path", multiplexed_qza_file_path,
    ]
    run_cmd(command)
    return multiplexed_qza_file_path

def qiime_import_project(dir_path, manifest_path, paired):
    multiplexed_qza_file_path = os.path.join(dir_path, "qza", "merged-demux.qza")

    command = [
        "qiime", "tools", "import",
        "--type", f"{'SampleData[PairedEndSequencesWithQuality]' if paired else 'SampleData[SequencesWithQuality]'}",
        "--input-path", manifest_path,
        "--output-path", multiplexed_qza_file_path,
        "--input-format", "PairedEndFastqManifestPhred33V2"
    ]
    run_cmd(command)
    return multiplexed_qza_file_path


def qiime_demux(dir_path,multiplexed_qza_file_path, metadata_path, paired):

    demux_qza_file_path = os.path.join(dir_path,"qza", f"demux-{'paired' if paired else 'single'}-end.qza")
    untrim_qza_file_path = os.path.join(dir_path, "qza","untrimmed.qza")


    command = [
        "qiime", f"{'demux' if paired else 'cutadapt'}", f"{'emp-paired' if paired else 'demux-single'}",
        "--i-seqs", multiplexed_qza_file_path,
        "--m-barcodes-file", metadata_path,
        "--m-barcodes-column", "barcode",
        f"{'--p-no-golay-error-correction'if paired else '--p-error-rate 0'}",
        "--o-per-sample-sequences", demux_qza_file_path,
        f"{'--o-error-correction-details' if paired else '--o-untrimmed-sequences'}", untrim_qza_file_path,
        "--verbose"
    ]
    run_cmd(command)
    return demux_qza_file_path

def qiime_demux_project(dir_path, multiplexed_qza_file, paired):
    manifest_list = []
    for directory in os.listdir(dir_path):
        print("directory", directory)
        if "directory" not in directory:
            print("continue", directory)
            continue
        if len(os.listdir(os.path.join(dir_path, directory, "fastq"))) == 3:
            demux_qza_file_path = os.path.join(dir_path, directory, "demux-paired-end.qza")
            metadata_path = os.path.join(dir_path, directory, "metadata", "metadata.tsv")
            untrim_qza_file_path = os.path.join(dir_path, directory, "untrimmed.qza")
            command = [
                "qiime", f"{'demux' if paired else 'cutadapt'}", f"{'emp-paired' if paired else 'demux-single'}",
                "--i-seqs", multiplexed_qza_file,
                "--m-barcodes-file", metadata_path,
                "--m-barcodes-column", "barcode",
                f"{'--p-no-golay-error-correction' if paired else '--p-error-rate 0'}",
                "--o-per-sample-sequences", demux_qza_file_path,
                f"{'--o-error-correction-details' if paired else '--o-untrimmed-sequences'}", untrim_qza_file_path,
                "--verbose"
            ]
            run_cmd(command)

            output_path = os.path.join(dir_path, directory, "exported")
            command = ["qiime", "tools", "export", "--input-path", demux_qza_file_path, "--output-path", output_path]
            run_cmd(command)

            manifest_list.append(generate_manifest(output_path, "manifest.tsv"))
        else:
            print("not enough files", directory)
            continue
    return combine_tsv_files(manifest_list, os.path.join(dir_path, "all", "manifest.tsv"))


def check_metadata(metadata_path):

    metadata = pd.read_csv(metadata_path, sep='\t')
    if 'barcode' in metadata.columns:
        return "yes"
    return "no"


def trim_single(dir_path,demux_qza_file_path, paired):
    #Trim adapters from demultiplexed reads
    #If there are sequencing adapters or PCR primers in the reads which you'd like to remove, you can do that next as follows.

    trimmed_seqs_file_path = os.path.join(dir_path,"qza", "trimmed-seqs.qza")
    command = [
        "qiime", "cutadapt", f"trim-{'paired' if paired else 'single'}",
        "--i-demultiplexed-sequences", demux_qza_file_path,
        f"{'--p-front-f' if paired else '--p-front'}", "GCTACGGGGGG",
        "--p-error-rate", "0",
        "--o-trimmed-sequences", trimmed_seqs_file_path,
        "--verbose"
    ]
    run_cmd(command)
    return trimmed_seqs_file_path


def qiime_summarize(dir_path,trimmed_seqs_file_path):
    #Summarize demultiplexed and trimmed reads
    vis_path = os.path.join(dir_path, "vis")
    vis_file_path = os.path.join(vis_path, "trimmed-seqs.qzv")

    command=[
        "qiime", "demux", "summarize",
        "--i-data", trimmed_seqs_file_path,
        "--o-visualization", vis_file_path
    ]
    run_cmd(command)
    return vis_file_path


def get_reads_data(dir_path,demux_qza_file_path, paired):

    output_dir_path = os.path.join(dir_path,"extracted-reads")

    #extract data from demultiplexed-seqs.qza file
    command = [
        "qiime", "tools", "extract",
        "--input-path", demux_qza_file_path,
        "--output-path", output_dir_path
    ]
    run_cmd(command)

    os.chdir(output_dir_path)
    subdirectories = [d for d in os.listdir() if os.path.isdir(d)]
    if not subdirectories:
        print("Something wrong with the data. Its doesn't follow the rules of qiita.")
        return
    else:
        subdirectory_path = subdirectories[0]

    os.chdir(subdirectory_path)

    if paired:
        return ReadsData(dir_path, fwd=True, rev=True)

    else:
        return ReadsData(dir_path, fwd=True, rev=False)



def qiita_visualization(fastq_path,metadata_path,data_type, verbose):
    verbose_print = print if verbose else lambda *a, **k: None

    verbose_print("\n")
    verbose_print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    dir_path = os.path.commonpath([os.path.dirname(fastq_path), os.path.dirname(metadata_path)])

    paired = check_paired_flag(fastq_path)
    verbose_print("\n")
    verbose_print('Checking environment...', end=" ")
    check_conda_qiime2()
    verbose_print('Done.')
    verbose_print("\n")
    verbose_print('Checking metadata...', end= " ")
    if check_metadata(metadata_path) == "no":
        verbose_print("The 'barcode' column does not exist in metadata.tsv. check and try again.")
        return
    verbose_print('Done.')
    verbose_print("\n")
    run_cmd(["mkdir", os.path.join(dir_path, "qza")])
    run_cmd(["mkdir", os.path.join(dir_path, "vis")])


    verbose_print("Find ALL NEW data in the directory you created:", dir_path)

    if data_type == '16S' or data_type == '18S':

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Import the multiplexed sequences' (1/4)")
        multiplexed_qza_file_path= qiime_import(dir_path,fastq_path, paired)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Import the multiplexed sequences' (1/4)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Demultiplex the reads' (2/4)")
        demux_qza_file_path=qiime_demux(dir_path, multiplexed_qza_file_path,metadata_path, paired)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Demultiplex the reads' (2/4)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Trim adapters from demultiplexed reads' (3/4)")
        trimmed_seqs_file_path=trim_single(dir_path,demux_qza_file_path, paired)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'rim adapters from demultiplexed reads' (3/4)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Summarize demultiplexed and trimmed reads' (4/4)")
        vis_file_path= qiime_summarize(dir_path,trimmed_seqs_file_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Summarize demultiplexed and trimmed reads' (4/4)")

        #getting values about fwd and rev
        reads_data= get_reads_data(dir_path,demux_qza_file_path, paired)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating visualization\n")

        pickle.dump(reads_data, open(os.path.join(reads_data.dir_path, "reads_data.pkl"), "wb"))

        print(f"Visualization file is located in {vis_file_path}\n"
              f"Please drag this file to https://view.qiime2.org/ and continue.\n")
        if reads_data.fwd and reads_data.rev:
            print(f"Note: The data has both forward and reverse reads.\n"
                  f"Therefore, you must give the parameters 'trim' and 'trunc' of export() "
                  f"as a tuple of two integers."
                  f"The first place related to the forward read and the second to the reverse.")
        else:
            print(f"Note: The data has only a forward read.\n"
                  f"Therefore, you must give the parameters 'trim' and 'trunc' of export() "
                  f"exactly one integers value which is related to the forward read.")

        return reads_data.dir_path

    else:

        print("YaMAS doesnt support downloading shotgun, yet.")

def qiita_project_visualization(dir_path, data_type, paired, verbose):
    verbose_print = print if verbose else lambda *a, **k: None

    if not paired:
        print("YaMAS currently only supports paired-end reads.")
        return

    verbose_print("\n")
    verbose_print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    os.makedirs(os.path.join(dir_path, "qza"), exist_ok=True)
    os.makedirs(os.path.join(dir_path, "vis"), exist_ok=True)
    verbose_print("\n")
    verbose_print('Checking environment...', end=" ")
    check_conda_qiime2()
    verbose_print('Done.')
    verbose_print("\n")
    run_cmd(["mkdir", os.path.join(dir_path, "qza")])
    run_cmd(["mkdir", os.path.join(dir_path, "vis")])
    verbose_print("Find ALL NEW data in the directory you created:", dir_path)

    if data_type == '16S' or data_type == '18S':
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'merging the fastq.gz files' (1/6)")
        main_path = merge_fastq(dir_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'merging the fastq.gz files' (1/6)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Importing merged fastq' (2/6)")
        multiplexed_qza_file_path = qiime_import(dir_path, main_path, paired)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Importing merged fastq' (2/6)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Analyzing sub-projects' (3/6)")
        combined_manifest = qiime_demux_project(dir_path, multiplexed_qza_file_path, paired)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Analyzing sub-projects' (3/6)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Importing project' (4/6)")
        multiplexed_qza_file_path = qiime_import_project(dir_path, combined_manifest, paired)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Importing project' (4/6)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Trimming the reads' (5/6)")
        trimmed_seqs_file_path = trim_single(dir_path, multiplexed_qza_file_path, paired)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Trimming the reads' (5/6)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Summarizing the reads' (6/6)")
        vis_file_path = qiime_summarize(dir_path, trimmed_seqs_file_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Summarizing the reads' (6/6)")

        print(f"Visualization file is located in {vis_file_path}\n"
              f"Please drag this file to https://view.qiime2.org/ and continue.\n")
        print(f"Note: The data has both forward and reverse reads.\n"
              f"Therefore, you must give the parameters 'trim' and 'trunc' of export() "
              f"as a tuple of two integers."
              f"The first place related to the forward read and the second to the reverse.")