import os
import wget
import re
from bs4 import BeautifulSoup
import requests
import csv


# This function extracts specified information for a given species from the main csv file
def extract_csv(csv_file, species_name, info_type):
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['species_name'].lower() == species_name.lower():
                if info_type == 't_number':
                    return row['t_number']
                elif info_type == 'code':
                    return row['code']
                elif info_type == 'genome_link':
                    return row['genome_link']
                else:
                    return "Invalid info_type provided. Please choose from 't_number', 'species_kegg_code', or 'species_archive_link'."


# This program extracts FNA, GTF, and general table files from the NCBI database
# It extracts only the reference files used by KEGG and that contain the local KEGG gene IDs
def find_gtf_link(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        gtf_link = soup.find(href=re.compile(r'^GCA_.*?\.gtf\.gz$'))

        if gtf_link:
            return gtf_link['href']

    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")

    return None


def find_gff_link(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        gff_link = soup.find(href=re.compile(r'^GCA_.*?\.gff\.gz$'))

        if gff_link:
            return gff_link['href']

    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")

    return None


def find_table_link(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        table_link = soup.find(href=re.compile(r'^GCA_.*?_table\.txt\.gz$'))
        if table_link:
            return table_link['href']

    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")

    return None


def find_second_fna_link(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        fna_links = soup.find_all(href=re.compile(r'genomic\.fna\.gz$'))

        if len(fna_links) >= 2:
            return fna_links[1]['href']

    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")

    return None

# This function downloads the specified file type for a given species
def download_file(species_name, file_type, output_directory, main_csv):
    base_link = extract_csv(main_csv, species_name, 'genome_link')

    # Extract the file links and download the file
    if file_type == "fna":
        file_link = find_second_fna_link(base_link)
    elif file_type == "table":
        file_link = find_table_link(base_link)
    elif file_type == "gtf":
        file_link = find_gtf_link(base_link)
    elif file_type == "gff":
        file_link = find_gff_link(base_link)
    else:
        print("Invalid file type. Use 'fna' or 'gtf' or 'table'.")
        return None

    # Step 5: Create the species folder if it doesn't exist
    species_folder = os.path.join(output_directory, species_name)
    if not os.path.exists(species_folder):
        os.makedirs(species_folder)

    # Step 6: Download the file and extract it
    if file_type == "fna":
        output_file = f"{species_name}.fna.gz"  # Set the filename with .fna.gz suffix
    elif file_type == "table":
        output_file = f"{species_name}.txt.gz"  # Set the filename with .txt.gz suffix
    elif file_type == "gtf":
        output_file = f"{species_name}.gtf.gz"  # Set the filename with .gtf.gz suffix
    elif file_type == "gff":
        output_file = f"{species_name}.gff.gz"  # Set the filename with .gtf.gz suffix
    else:
        print("Invalid file type. Use 'fna' or 'gtf' or 'table.")
        return None

    file_path = os.path.join(species_folder, output_file)
    file_url = f"{base_link}/{file_link}"

    try:
        wget.download(file_url, out=file_path)
    except Exception as e:
        print(f"couldn't download fna\gtf. this is the error {e}")
        return

    # Use the full filename for the gunzip command
    gz_file_path = file_path

    # Change working directory to the species folder before extraction
    os.chdir(species_folder)

    # Extract the gz file using gunzip command
    os.system(f"gunzip {gz_file_path}")

    if file_type == "fna":
        print(f"{species_name} FNA file downloaded and extracted")
    elif file_type == "table":
        print(f"{species_name} table file downloaded and extracted")
    elif file_type == "gtf":
        print(f"{species_name} GTF file downloaded and extracted")
    elif file_type == "gff":
        print(f"{species_name} GFF file downloaded and extracted")

    if file_link:
        return f"{base_link}/{file_link}"