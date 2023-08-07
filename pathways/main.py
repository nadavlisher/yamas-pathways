import os
from metaphlan import install_metaphlan, install_metaphlan_database, run_metaphlan
import multiprocessing
from findspecies import search_species
from extract_file_KEGG import download_file
from combine_fna import combine_fna_files
import rpy2.robjects as robjects
from reference_data2 import reference_data
from more_pathways import fetch_more_gene_pathway_data
from remove_zero import remove_zeros_from_count_matrix
import json

# Create a folder where you want all the data (not including database downloads) from this session to be stored
main_folder_name = "data_folder_2"
parent_directory = "/home/yisrael/Desktop"
main_folder_path = os.path.join(parent_directory, main_folder_name)
if not os.path.exists(main_folder_name):
    os.makedirs(main_folder_name)

# Create a folder where you want the databases (FNA/GTF) to be stored
# If you already have a database of FNA/GTF annotation files from a previous run,
# Then indicate its location and don't create a new file
reference_folder_name = "index_files"
parent_directory = "/home/yisrael/Desktop"
reference_folder_path = os.path.join(parent_directory, reference_folder_name)
if not os.path.exists(reference_folder_path):
    os.makedirs(reference_folder_path)

# Set a path to the dictionary of specie codes, the csv with species information, and the KEGG species list
dictionary_path = "/home/yisrael/Desktop/pathways_files/species_codes.json"
main_csv = "/home/yisrael/Desktop/pathways_files/species2.csv"
kegg_species = "/home/yisrael/Desktop/pathways_files/all_KEGG_species.txt"

# Indicate the current location of the input fastq files below (if they aren't paired, then remove the second)
fastq_file_1 = "/home/yisrael/Desktop/fastq_files/SRR2142192_1.fastq"
fastq_file_2 = "/home/yisrael/Desktop/fastq_files/SRR2142192_2.fastq"

# If you don't have metaphlan and its databases installed, then use the following functions
# Specify a location to hold the database regardless (or where they are already being held)
#metaphlan_database_folder = "/home/yisrael/metadata"
#install_metaphlan()
#install_metaphlan_database(metaphlan_database_folder)

# Run metaphlan using the specified folders
# If running the program from pycharm takes up too much memory, you can run it from the terminal with this command
# metaphlan </path/to/fastq> --bowtie2out metagenome.bowtie2.bz2 --nproc 2 --input_type fastq -o </path/to/main/folder/profiled_metagenome.txt> --bowtie2db </path/to/data_base>
# If you install from terminal, use the second function and specify the location of the output
# metaphlan_output = run_metaphlan(fastq_file_1, metaphlan_database_folder, main_folder_path)
metaphlan_output = "/home/yisrael/Desktop/data_folder_2/profiled_metagenome.txt"

# Create a list containing the names of the unique species identified
species_list = search_species(metaphlan_output, main_folder_path, kegg_species)
species_list.sort()

# Build a reference file for reading and alignment
# The FNA (fasta) and GTF file of each species will be downloaded from the NCBI database
# If they have been downloaded from a previous session, they will not be downloaded again
for species in species_list:
    species_directory = os.path.join(reference_folder_path, species)
    if not os.path.exists(species_directory):
        # Download the fna/gtf file only if the species folder doesn't exist
        download_file(species, "fna", reference_folder_path, main_csv)
        download_file(species, "gtf", reference_folder_path, main_csv)
        download_file(species, "gff", reference_folder_path, main_csv)
    else:
        print(f"Skipping {species}. Folder already exists in the output directory.")

# The FNA files of each species will be combined into one large annotation file for alignment
annotation_file = combine_fna_files(species_list, reference_folder_path, main_folder_path)

# The final steps of the processing are done in R
r_code = f'''
    library(Rsubread)
    library(clusterProfiler)

    # Set a folder path on the desktop
    folder_path <- "{main_folder_path}"

    # Set the paths to your input files
    buildindex(basename = "my_index", reference = "{annotation_file}")
    readfile1 <- "{fastq_file_1}"  # Path to the first FASTQ file
    readfile2 <- "{fastq_file_2}"  # Path to the second FASTQ file
    annotation_directory <- "{reference_folder_path}"  # Path to the annotation directory
    dictionary_path <- "{dictionary_path}"  # Path to the json dictionary file

    # Load species_list from a .txt file in the folder path
    species_list <- scan(file.path(folder_path, "species_list.txt"), what = "character")

    # Align reads (execute only once)
    align(index = "my_index", readfile1 = readfile1, readfile2 = readfile2, type = "dna", output_file = "aligned_reads.bam")

    # Loop through species_list
    for (species_name in species_list) {{
        # Create a subfolder for the current species
        species_folder <- file.path(folder_path, species_name)
        dir.create(species_folder, showWarnings = FALSE)
    
        # Build the gtf_file path for each species
        gtf_file_subdirectory <- file.path(annotation_directory, species_name)
        gtf_file <- file.path(gtf_file_subdirectory, paste0(species_name, ".gtf"))
  
        # Load the aligned reads and create gene id list
        countData <- featureCounts(files = "aligned_reads.bam", annot.ext = gtf_file, isGTFAnnotationFile = TRUE, isPairedEnd = TRUE)
        countMatrix <- countData$counts
        print(countMatrix)
  
        # Specify the file path and name for the output CSV file
        output_file <- file.path(species_folder, paste0(species_name, "_countMatrix.csv"))
  
        # Save the count matrix as a CSV file
        write.csv(countMatrix, file = output_file, row.names = TRUE)
  
        countMatrix2 <- read.csv(output_file)
        
        # Remove all entries with a value of 0
        filtered_matrix <- countMatrix2[countMatrix2[,2] != 0, ]
  
        # Convert count matrix to Entrez gene IDs
        genes <- filtered_matrix[, 1]
        genes
  
        # Perform pathway analysis
        # Read the json dictionary file to get the organism code
        json_data <- jsonlite::fromJSON(dictionary_path)
        organism <- json_data[[species_name]]
  
        pwRes <- enrichKEGG(gene = genes,
                organism = organism,
                keyType = "kegg",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)
  
        # Get the enriched pathways and associated genes
        pathways <- pwRes$ID
        genes <- pwRes$geneID
  
        output_file <- file.path(species_folder, paste0(species_name, "_pathway_data.txt"))
        cat("Enriched Pathways:\n", file = output_file)
        for (i in seq_along(pathways)) {{
            cat("\nPathway:", pathways[i], "\n", "Associated Genes:", genes[i], "\n", file = output_file, append = TRUE)
        }}
  
        # Save the enriched pathways to a file
        output_file_pwRes <- file.path(species_folder, paste0(species_name, "_pwRes_output.txt"))
        pwRes_text <- capture.output(print(pwRes))
        writeLines(pwRes_text, con = output_file_pwRes)
    }}
'''
robjects.r(r_code)

# The following section will use parallel processing to download the reference species data faster
def process_species_subset(species_list, annotation_folder_path, dictionary_path, start_idx, end_idx):
    subset_species_list = species_list[start_idx:end_idx]
    reference_data(subset_species_list, annotation_folder_path, dictionary_path)

if __name__ == "__main__":
    num_processes = 3    # Number of parallel processes you want to run

    # Divide the species list into equal parts for each process
    chunk_size = len(species_list) // num_processes
    processes = []

    for i in range(num_processes):
        start_idx = i * chunk_size
        end_idx = (i + 1) * chunk_size if i < num_processes - 1 else len(species_list)

        # Create a process for each subset of the species list
        p = multiprocessing.Process(target=process_species_subset,
                                    args=(species_list, reference_folder_path , main_csv, start_idx, end_idx))
        processes.append(p)
        p.start()

    # Wait for all processes to finish
    for p in processes:
        p.join()

    print("All processes have finished.")