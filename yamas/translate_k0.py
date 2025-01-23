import requests
from bs4 import BeautifulSoup
import math
import re
from .extract_file_KEGG import extract_csv
import pandas as pd
from ete3 import NCBITaxa
import os
from collections import defaultdict


def generate_pathway_csv_from_sample(sample_dir):
    """
    Extract pathways for each species in the given sample directory.
    Generates a CSV file with species as rows, pathways as columns, and binary presence/absence values.

    Parameters:
    - sample_dir (str): Path to the sample directory.

    Returns:
    - None: Writes a CSV file with species as rows, pathways as columns, and binary presence/absence values.
    """
    # Initialize a dictionary to hold species-pathway mapping
    pathway_dict = defaultdict(set)

    # Loop through each species directory in the sample directory
    for species_dir in os.listdir(sample_dir):
        species_path = os.path.join(sample_dir, species_dir)
        if not os.path.isdir(species_path):
            continue  # Skip if not a directory

        # Define the pathway data file for each species
        pathway_file = os.path.join(species_path, f"{species_dir}_pwRes_output.txt")
        species_name = species_dir  # the species name matches the directory name

        # Check if the pathway data file exists and extract pathways
        if os.path.isfile(pathway_file):
            # Skip empty files
            if os.stat(pathway_file).st_size == 0:
                print(f"Skipping empty file: {pathway_file}")
                continue

            # Read the file as a DataFrame
            try:
                df = pd.read_csv(pathway_file, sep="\t", header=0)
            except Exception as e:
                continue  # Skip this file and move to the next one

            # Extract the pathway name (before the " - ")
            for description in df["Description"]:
                pathway_name = description.split(" - ")[0]  # Get the pathway name without the species name
                pathway_dict[species_name].add(pathway_name)

        # Check if any pathways were found
    if not pathway_dict:
        print(f"No pathway data found in the sample directory: {sample_dir}")
        return

        # Consolidate all pathway names across all species
    all_pathways = sorted({pathway for pathways in pathway_dict.values() for pathway in pathways})

    # Create a DataFrame with species as rows and pathway names as columns
    data = []
    for species, pathways in pathway_dict.items():
        row = {pathway: 1 if pathway in pathways else 0 for pathway in all_pathways}
        row["Species"] = species
        data.append(row)

    # Convert to DataFrame and reorder columns
    df = pd.DataFrame(data)
    if "Species" not in df.columns:
        print("Error: 'Species' column is missing from the DataFrame. Check the data processing logic.")
        return

    df = df.set_index("Species").reindex(columns=all_pathways, fill_value=0)

    # Save the DataFrame to a CSV file
    output_file = os.path.join(sample_dir, f"pathways_data.csv")
    df.to_csv(output_file)



# Define a function to get full taxonomy hierarchy
def get_full_hierarchy(species_name):
    """
    This function uses ete3 in order to get the full hierarchy of a given species.
    :param species_name: species name
    :return: the full hierarchy of the given species
    """
    ncbi = NCBITaxa()
    try:
        # Convert the species name to a format suitable for searching
        formatted_name = species_name.replace("_", " ")

        # Look up taxid for the species
        taxid_dict = ncbi.get_name_translator([formatted_name])

        # Check if the taxid was found
        if formatted_name in taxid_dict:
            taxid = taxid_dict[formatted_name][0]

            # Get the lineage of the taxid
            lineage = ncbi.get_lineage(taxid)

            # Get the names and ranks for the lineage
            names = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)

            # Construct the taxonomy in the desired format
            hierarchy = []
            for rank in ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]:
                if rank in ranks.values():
                    rank_taxid = next(key for key, value in ranks.items() if value == rank)
                    rank_code = rank[0].upper()
                    hierarchy.append(f"{rank_code}__{names[rank_taxid]}")
            return ";".join(hierarchy)
        else:
            print(f"Taxonomy not found for species: {species_name}")
            return species_name
    except Exception as e:
        print(f"Error processing species '{species_name}': {e}")
        return "Error"


def extract_genes_and_k0(t_code, species_code, info_type):
    """
    Extracts gene information and KEGG Orthology (KO) numbers from the KEGG genome database.
    :param t_code: The taxonomic code used in the KEGG database for the organism of interest.
    :param species_code: The species-specific identifier code used in the gene IDs. Can be either 3 or 4 characters long.
    :param info_type: Specifies the type of information to return. Must be one of:
        - 'genes': Returns list of gene IDs
        - 'K0s': Returns list of KO numbers
        - 'dictionary': Returns a dictionary mapping genes to their KO numbers
    """
    code = f"{species_code}"
    code_length = len(code)

    # Create the base link using the species code
    genes_link = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:{t_code}"

    response = requests.get(genes_link)
    soup = BeautifulSoup(response.text, 'html.parser')
    hits_tag = soup.find(string=lambda text: 'Hits:' in text)
    num_hits = int(re.search(r'\d+', hits_tag).group())
    num_pages = math.ceil(num_hits / 1000)

    extracted_genes = []
    extracted_k0s = []

    # Extract the text from all pages and save it
    for page_number in range(1, num_pages + 1):
        url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+-p+{page_number}+gn:{t_code}"
        response = requests.get(url)
        response.raise_for_status()  # Check for any errors in the request

        # Assuming the website uses UTF-8 encoding for text
        text_content = response.text

        # Remove lines containing "no KO assigned"
        lines = [line for line in text_content.splitlines() if "<a href=\"/entry/" in line
                 and "no KO assigned" not in line]

        for line in lines:
            if code_length == 3:
                start_index_gene = line.find(f"{species_code}:") + 4
            else:
                start_index_gene = line.find(f"{species_code}:") + 5
            end_index_gene = line.find('"', start_index_gene)
            if 0 <= start_index_gene < end_index_gene:
                extracted_gene = line[start_index_gene:end_index_gene]
                extracted_gene = extracted_gene.strip()
                extracted_genes.append(extracted_gene)

            start_index_k0 = line.find("  K") + 2
            end_index_k0 = line.find(" ", start_index_k0)
            if 0 <= start_index_k0 < end_index_k0:
                extracted_k0 = line[start_index_k0:end_index_k0]
                extracted_k0 = extracted_k0.strip()
                extracted_k0s.append(extracted_k0)

    if info_type == "genes":
        return extracted_genes
    if info_type == "K0s":
        return extracted_k0s
    if info_type == "dictionary":
        gene_k0_dict = {}  # Initialize an empty dictionary

        if len(extracted_genes) == len(extracted_k0s):
            for i in range(len(extracted_genes)):
                gene = extracted_genes[i]
                k0 = extracted_k0s[i]
                gene_k0_dict[gene] = k0  # Map gene to K0 in the dictionary
            return gene_k0_dict
        else:
            print("Lists of genes and K0s have different lengths.")
            return None


def translate_genes_to_k0s(species_name, count_matrix_filename, csv_filename, verbose_print):
    """
    Translates gene names to K0 identifiers and aggregates counts based on K0 identifiers.

    Parameters:
        species_name (str): Name of the species to process.
        count_matrix_filename (str): File containing the gene count matrix.
        csv_filename (str): File containing the gene-K0 mapping information.

    Returns:
        list: A list of translated lines in the format "k0_name total_count".
    """
    t_code = extract_csv(csv_filename, species_name, "t_number")
    species_code = extract_csv(csv_filename, species_name, "code")

    gene_k0_dict = extract_genes_and_k0(t_code, species_code, "dictionary")

    if gene_k0_dict is None:
        print("Failed to create the gene-K0 dictionary.")
        return None

    translated_lines = []
    invalid_lines = []

    with open(count_matrix_filename, 'r') as matrix_file:
        gene_k0_counts = {}

        for line in matrix_file:
            # Split the line and check if it has exactly two parts
            parts = line.strip().split()
            if len(parts) != 2:
                # Log invalid lines for debugging
                invalid_lines.append(line.strip())
                continue

            gene_name, count = parts
            if gene_name in gene_k0_dict:
                k0_name = gene_k0_dict[gene_name]
                if k0_name in gene_k0_counts:
                    gene_k0_counts[k0_name] += int(count)
                else:
                    gene_k0_counts[k0_name] = int(count)
            else:
                verbose_print(f"Skipped gene: {gene_name} due to lack of translation to K0.")

        # Generate the translated lines
        for k0_name, total_count in gene_k0_counts.items():
            translated_line = f"{k0_name} {total_count}"
            translated_lines.append(translated_line)

    # Save the translated lines to a new file
    translated_filename = count_matrix_filename.replace(".txt", "_translated.txt")
    with open(translated_filename, 'w') as translated_file:
        for line in translated_lines:
            translated_file.write(line + '\n')

    print(f"Translated count matrix saved to: {translated_filename}")
    return translated_lines


def translate_k0_to_gene(input_csv, output_csv):
    """
    This function translates K0 names to gene names in the column headers of the CSV file
    generated by 'central_pipeline.py' (large_spreadsheet.csv).
    It takes the input CSV file location and the output CSV file location as parameters.
    """

    def translator(k0_name):
        base_url = "https://www.genome.jp/dbget-bin/www_bfind_sub"
        query_params = {
            'mode': 'bfind',
            'max_hit': '1000',
            'locale': 'en',
            'serv': 'kegg',
            'dbkey': 'orthology',
            'keywords': k0_name
        }
        # Step 1: Fetch the webpage content
        headers = {'User-Agent': 'Mozilla/5.0'}
        response = requests.get(base_url, params=query_params, headers=headers)

        if response.status_code == 200:
            # Step 2: Parse the HTML content
            soup = BeautifulSoup(response.text, 'html.parser')

            # Step 3: Find all <div> elements with style 'margin-left:2em'
            div_elements = soup.find_all('div', style="margin-left:2em")

            # Step 4: Extract the text from each matching <div> and return the first match
            div_texts = [div.get_text(strip=True) for div in div_elements]

            # Return the first matched text or the original K0 name if not found
            if div_texts:
                div_text = div_texts[0]  # Get the first match
                if ";" in div_text:
                    # Split by semicolon and return the part after it
                    return div_text.split(";", 1)[1].strip()
                else:
                    # If there's no semicolon, return the original text
                    return div_text
            else:
                # Return the K0 name if no matching div is found
                return k0_name
        else:
            return f"Failed to retrieve data for {k0_name}. Status code: {response.status_code}"

    df = pd.read_csv(input_csv)

    # Translate each column header (K0 name)
    translated_columns = []
    for col in df.columns:
        gene_info = translator(col)
        if gene_info:
            translated_columns.append(gene_info)
        else:
            translated_columns.append(col)  # If no translation, keep the original K0 name

    # Update the column headers with translated names
    df.columns = translated_columns

    # Save the DataFrame with the translated column headers and data to the new CSV
    df.to_csv(output_csv, index=False)
    print(f"Translated CSV saved to: {output_csv}")
