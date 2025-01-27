import rpy2.robjects as robjects
import os
from .extract_file_KEGG import download_file, extract_csv
from .reference_data import reference_data
import multiprocessing
import pandas as pd
from .translate_k0 import translate_k0_to_gene, translate_genes_to_k0s, get_full_hierarchy, generate_pathway_csv_from_sample
from .utilities import run_cmd



def search_species(filename, output_directory, dictionary_path):
    """
    This function using certain flags to find the unique species within the sample
    It won't count any species that don't have files on the NCBI database
    :param filename: metaphlan output (profile.txt)
    :param output_directory: where to save the species list file
    :param dictionary_path: path to the KEGG file
    :return: creates a file containing species in the sample
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    comparison_filename = dictionary_path
    with open(comparison_filename, 'r') as comparison_file:
        comparison_lines = comparison_file.readlines()

    comparison_items = set()
    for line in comparison_lines:
        comparison_items.add(line.strip())

    items = set()
    for line in lines:
        start_index = 0
        while True:
            start_index = line.find('|s__', start_index)
            if start_index == -1:
                break
            end_index = line.find('2|', start_index)
            if end_index == -1:
                break
            item = line[start_index + 3: end_index].strip()
            if '|t__' in item:
                item = item.split('|t__')[0]
            if item.startswith('_'):
                item = item[1:]
            if item in comparison_items:
                items.add(item)
            start_index = end_index + 2

    if items:
        max_length = max(len(item) for item in items)
        formatted_items = [item.ljust(max_length) for item in items]

        output_file = f"{output_directory}/species_list.txt"
        with open(output_file, 'w') as output_file:
            for item in formatted_items:
                output_file.write(item + "\n")

    return list(items)


# This function takes all the FNA files and combines them into one file for alignment
def combine_fna_files(species_list, annotation_directory, output_file_path, verbose_print):
    """
    This function creates a combined fna file for the sample with only the bacteria found in it. This file is used for the analysis and indexing in the R code.
    :param verbose_print: prints the line if wanted
    :param species_list: species list found via the pathways function
    :param annotation_directory:  directory where the data was downloaded from the download_file function
    :param output_file_path: the path to the saved combined fna file
    :return: creates a combined fna file containing all the fna files of the bacteria found in the sample
    """
    try:
        combined_fna_path = os.path.join(output_file_path, "combined_fna.fna")

        with open(combined_fna_path, "w") as output_file:
            for species in species_list:
                species_directory = os.path.join(annotation_directory, species)
                species_fna_file_path = os.path.join(species_directory, f"{species}.fna")

                # Check if the fna file for the species exists
                if os.path.exists(species_fna_file_path):
                    with open(species_fna_file_path, "r") as fna_file:
                        # Read the contents of the fna file and write it to the output file
                        output_file.write(fna_file.read())
                    verbose_print(f"FNA file for {species} added to the output file.")
                else:
                    verbose_print(f"Warning: FNA file for {species} not found.")

        verbose_print(f"All fna files have been combined into {combined_fna_path}.")
        return combined_fna_path

    except Exception as e:
        print(f"Error combining fna files: {e}")
        return None


# The following section will use parallel processing to quickly download the reference species data faster
# This data gives a list of pathways and the genes they contain
def process_species_subset(species_list, annotation_folder_path, dictionary_path, start_idx, end_idx):
    subset_species_list = species_list[start_idx:end_idx]
    reference_data(subset_species_list, annotation_folder_path, dictionary_path)


def create_spreadsheet(bacteria_list, main_folder_path, save_path):
    """
    Creates a csv file with species as rows and genes as columns, with the gene counts in between.
    :param bacteria_list: bacteria list from the pathways code
    :param main_folder_path: the path where the data is stored
    :param save_path: the path where the spreadsheet will be saved
    :return: creates a spreadsheet containing the gene counts of the data
    """
    species_data = {}  # Dictionary to store species data

    # Loop through each species
    for species_name in bacteria_list:
        species_folder = os.path.join(main_folder_path, species_name)
        count_matrix_filename = os.path.join(species_folder, f"{species_name}_gene_counts_translated.txt")

        if os.path.exists(count_matrix_filename):
            # Read the count matrix file for the current species
            with open(count_matrix_filename, 'r') as matrix_file:
                genes_counts = [line.strip().split() for line in matrix_file]

            # Extract genes and counts for the current species
            genes = [gene_count[0] for gene_count in genes_counts]
            counts = [int(gene_count[1]) for gene_count in genes_counts]

            # Add species data to the dictionary
            species_data[species_name] = dict(zip(genes, counts))

    # Create a DataFrame from the collected species data
    main_spreadsheet = pd.DataFrame(species_data)

    # Transpose the DataFrame to have genes as rows and bacteria as columns
    main_spreadsheet = main_spreadsheet.transpose()

    # Fill NaN values with 0
    main_spreadsheet = main_spreadsheet.fillna(0)

    # Save the main spreadsheet to the specified save path
    main_spreadsheet.to_csv(save_path)


def pathways(metaphlan_output, sample_name, reads_data, main_dir, verbose_print):
    """
    This is the main pipeline of pathways. This function does the following:
    1. creates parameters and creates folder of the data
    2. downloads for the species found data using the species2.csv file - all the files will be deleted later
    3. starts the R analysis
    4. final processing of the data
    creates the following files:
    - {species}_unique_genes.txt: Lists of unique genes per species
    - {species}_gene_counts.txt: Gene count data per species
    - {species}_pwRes_output.txt: Pathway enrichment results per species
    - gene_counts.csv: Consolidated gene counts across species
    - pathways_data.csv: Consolidated pathway analysis results
    """
    paired = reads_data.rev and reads_data.fwd
    if paired:
        fastq_file_1 = os.path.join(reads_data.dir_path, "fastq", f"{sample_name}_1.fastq")
        fastq_file_2 = os.path.join(reads_data.dir_path, "fastq", f"{sample_name}_2.fastq")
    else:
        fastq_file_1 = os.path.join(reads_data.dir_path, "fastq", f"{sample_name}")
        fastq_file_2 = ""
    pathway_dir = os.path.join(reads_data.dir_path, "pathways_data", sample_name.split(".")[0])
    data_dir = os.path.join(reads_data.dir_path, "pathways", sample_name.split(".")[0])
    os.makedirs(pathway_dir, exist_ok=True)
    kegg_species = os.path.join(main_dir, "all_KEGG_species.txt")
    main_csv = os.path.join(main_dir, "species2.csv")
    dictionary_path = os.path.join(main_dir, "species_codes.json")
    species_list = search_species(metaphlan_output, pathway_dir, kegg_species)
    if not species_list:
        print(f"for {sample_name} we could not find any species.")
        return
    species_list.sort()
    reference_path = os.path.join(reads_data.dir_path, "Index_files")
    os.makedirs(reference_path, exist_ok=True)
    for species in species_list:
        species_directory = os.path.join(reference_path, species)
        if not os.path.exists(species_directory):
            # Download the fna/gtf file only if the species folder doesn't exist
            download_file(species, "fna", reference_path, main_csv)
            download_file(species, "gtf", reference_path, main_csv)
        else:
            verbose_print(f"Skipping {species}. Folder already exists in the output directory.")

    annotation_file = combine_fna_files(species_list, reference_path, pathway_dir, verbose_print)

    # The following steps of the precessing will be done using R, this is the code
    # this code will do the following: (broad explanation)
    # 1. Align the reads to the reference annotation file
    # 2. Extract gene IDs from the aligned reads
    # 3. Perform pathway analysis using the extracted gene IDs
    # 4. Save the results to a file for each species
    r_code = f"""
    # Function to check and load required packages
    load_required_packages <- function(packages) {{
        for (pkg in packages) {{
            if (!pkg %in% .packages()) {{
                library(pkg, character.only = TRUE)
            }}
        }}
    }}

    # List of required packages
    required_packages <- c("Rsubread", "GenomicRanges", "GenomicAlignments", 
                           "rtracklayer", "clusterProfiler", "jsonlite")

    # Load the required packages
    load_required_packages(required_packages)

    # Set a folder path
    folder_path <- "{pathway_dir}"

    # Set the paths to your input files
    readfile1 <- "{fastq_file_1}"  # Path to the first FASTQ file
    readfile2 <- "{fastq_file_2}"  # Path to the second FASTQ file
    annotation_file <- "{annotation_file}"  # Reference annotation file
    annotation_directory <- "{reference_path}"  # Path to the annotation directory
    dictionary_path <- "{dictionary_path}"  # Path to the JSON dictionary file
    buildindex(basename = "my_index", reference = annotation_file)

    # Align reads only if the BAM file doesn't already exist
    bam_file <- file.path(folder_path, "aligned_reads.bam")
    if (!file.exists(bam_file)) {{
        if (readfile2 != "") {{
            align(index = "my_index", readfile1 = readfile1, readfile2 = readfile2, type = "rna", output_file = bam_file)
        }} else {{
            align(index = "my_index", readfile1 = readfile1, type = "rna", output_file = bam_file)
        }}
    }} else {{
        cat("BAM file already exists. Skipping alignment.\\n")
    }}

    # Load species_list from a .txt file in the folder path
    species_list <- scan(file.path(folder_path, "species_list.txt"), what = "character")

    # Loop through species_list
    for (species_name in species_list) {{
        # Create a subfolder for the current species
        species_folder <- file.path(folder_path, species_name)
        dir.create(species_folder, showWarnings = FALSE)

        # Build the GTF file path for each species
        gtf_file <- file.path(annotation_directory, species_name, paste0(species_name, ".gtf"))
        if (!file.exists(gtf_file)) {{
            cat("GTF file for species", species_name, "does not exist. Skipping to the next species.\\n")
            next
        }}

        # Load the aligned reads and create gene id list
        gr <- import(gtf_file, format = "gtf")
        bam_data <- readGAlignmentPairs(bam_file)
        overlaps <- findOverlaps(bam_data, gr)
        gene_ids <- mcols(gr)$gene_id[subjectHits(overlaps)]
        gene_ids <- gene_ids[!is.na(gene_ids)]

        # Get unique gene_ids
        unique_genes <- unique(gene_ids)

        # Count the occurrences of each gene_id
        gene_counts <- table(gene_ids)

        # Save the unique genes to a file
        unique_genes_output_file <- file.path(species_folder, paste0(species_name, "_unique_genes.txt"))
        write.table(unique_genes, file = unique_genes_output_file, col.names = FALSE, row.names = FALSE, quote = FALSE)

        # Save the gene counts to a file
        gene_counts_output_file <- file.path(species_folder, paste0(species_name, "_gene_counts.txt"))
        write.table(as.data.frame(gene_counts), file = gene_counts_output_file, col.names = TRUE, row.names = FALSE, quote = FALSE)

        # Perform pathway analysis
        json_data <- fromJSON(dictionary_path)
        organism <- json_data[[species_name]]

        # Create a KEGG pathway analysis
        pwRes <- enrichKEGG(gene = unique_genes,  
                            organism = organism,
                            keyType = "kegg",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)

        # Get the enriched pathways and associated genes
        enriched_data <- data.frame(
            Pathway = pwRes$ID,
            Description = pwRes$Description,
            GeneRatio = pwRes$GeneRatio,
            BgRatio = pwRes$BgRatio,
            pvalue = pwRes$pvalue,
            qvalue = pwRes$qvalue,
            geneID = pwRes$geneID
        )

        # Save the enriched pathways to a file
        output_file_pwRes <- file.path(species_folder, paste0(species_name, "_pwRes_output.txt"))
        write.table(enriched_data, file = output_file_pwRes, sep = "\\t", row.names = FALSE, quote = FALSE)
    }}
    """
    robjects.r(r_code)

    num_processes = 3  # Number of parallel processes you want to run

    # Divide the species list into equal parts for each process
    chunk_size = len(species_list) // num_processes
    processes = []

    # Start a process for each subset of the species list
    for i in range(num_processes):
        start_idx = i * chunk_size
        end_idx = (i + 1) * chunk_size if i < num_processes - 1 else len(species_list)

        # Create a process for each subset of the species list
        p = multiprocessing.Process(target=process_species_subset,
                                    args=(species_list, reference_path, main_csv, start_idx, end_idx))
        processes.append(p)
        p.start()

    # Wait for all processes to finish
    for p in processes:
        p.join()

    # Translate the unique gene names of each count matrix to KO codes, which are recognized for all species in the KEGG database
    for species in species_list:
        count_matrix = f"{pathway_dir}/{species}/{species}_gene_counts.txt"
        if not os.path.exists(count_matrix):
            verbose_print(f"skipping {species}")
            continue
        translate_genes_to_k0s(species, count_matrix, main_csv, verbose_print)

    # Create a spreadsheet containing the gene counts of the data and the pathways data
    gene_counts_csv = os.path.join(data_dir, "gene_counts.csv")
    create_spreadsheet(species_list, pathway_dir, gene_counts_csv)
    # translate_k0_to_gene(gene_counts_csv, gene_counts_csv)  #takes too much time currently.
    df = pd.read_csv(gene_counts_csv)
    df["Unnamed: 0"] = df["Unnamed: 0"].apply(get_full_hierarchy)
    df.to_csv(gene_counts_csv)
    generate_pathway_csv_from_sample(pathway_dir)
    if os.path.exists(os.path.join(data_dir, "pathways_data.csv)")):
        df = pd.read_csv(os.path.join(data_dir, "pathways_data.csv)"))
        df["Species"] = df["Species"].apply(get_full_hierarchy)
        df.to_csv(os.path.join(data_dir, "pathways_data.csv)"))

