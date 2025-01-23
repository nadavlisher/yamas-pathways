import os
from bs4 import BeautifulSoup
import time
from .extract_file_KEGG import extract_csv
import csv
import random
import requests


def extract_paths_from_pages(t_number, info_type):
    """
    This function extracts a list of pathways for a given species using its t_number
    :param t_number: t number of the species
    :param info_type: what info you want to return. must be:
        - 'pathways': Returns list of pathway identifiers
        - 'functions': Returns list of pathway functional descriptions
        - 'name': Alternative for 'pathways' (maintained for backwards compatibility)
        - 'function': Alternative for 'functions' (maintained for backwards compatibility)
    :return:
    """
    url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+pathway+gn:{t_number}"

    try:
        response = requests.get(url)
        response.raise_for_status()  # Check for any errors in the request

        # Assuming the website uses UTF-8 encoding for text
        text_content = response.text

        # Remove lines containing "no KO assigned"
        lines = [line for line in text_content.splitlines() if "<a href=\"/pathway/" in line]

        extracted_paths = []
        extracted_functions = []

        for line in lines:
            start_index_path = line.find("y/") + 2
            end_index_path = line.find('"', start_index_path)
            if 0 <= start_index_path < end_index_path:
                extracted_path = line[start_index_path:end_index_path]
                extracted_path = extracted_path.strip()
                extracted_paths.append(extracted_path)

            start_index_function = line.find("            ") + 12
            end_index_function = line.find(" -", start_index_function)
            if 0 <= start_index_function < end_index_function:
                extracted_function = line[start_index_function:end_index_function]
                extracted_function = extracted_function.strip()
                extracted_functions.append(extracted_function)


        if not extracted_paths:
            print("No data found for the given t_number.")
            return

        if info_type == "pathways":
            return extracted_paths
        elif info_type == "functions":
            return extracted_functions
        else:
            print("Invalid info_type provided. Please use 'pathways' or 'functions'.")

    except requests.exceptions.RequestException as e:
        print("Error occurred while fetching the URL:", e)
        return None

    # except requests.exceptions.RequestException as e:
       # print("Error occurred while fetching the URL:", e)

    if info_type == "name":
        return extracted_paths

    if info_type == "function":
        return extracted_functions

# This program creates a csv containing all pathways within a species and the genes associated with them
def reference_data(species_list, annotation_folder_path, main_csv):
    for species in species_list:
        print(f"Starting: {species}")
        species_directory = os.path.join(annotation_folder_path, species)

        t_number = extract_csv(main_csv, species, "t_number")

        # Create a list of pathways and their functions
        pathway_name_list = extract_paths_from_pages(t_number, "pathways")
        pathway_function_list = extract_paths_from_pages(t_number, "functions")
        pathway_data_file = os.path.join(species_directory, f"{species}_pathway_data.csv")

        # Fetch the actual data from the KEGG database
        fetch_gene_pathway_data(species_directory, pathway_name_list, pathway_function_list, pathway_data_file)
        print(f"Finished: {species}")


# This function fetches the gene data for a given pathway
def fetch_gene_pathway_data(output_directory, path_list, function_list, output_file):
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Create a new CSV file to store the pathway data
    output_file_path = os.path.join(output_directory, output_file)

    existing_paths = set()

    if os.path.exists(output_file_path):
        # Read the existing CSV file and extract the pathways
        with open(output_file_path, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            try:
                next(csv_reader)  # Skip the header row
            except StopIteration:
                # CSV file is empty, and there's no header row
                pass
            else:
                existing_paths = set(row[0] for row in csv_reader)

    with open(output_file_path, 'a', newline='') as output_file:
        csv_writer = csv.writer(output_file, delimiter=',')

        # If the file is empty or there's no header row, write the header row
        if output_file.tell() == 0:
            csv_writer.writerow(["Pathway", "Functions", "Genes"])  # Updated header row

        error_links = []

        for path, function in zip(path_list, function_list):
            if path in existing_paths:
                continue

            # Downloading quickly from the KEGG server can sometimes result in a 403 error
            # So the following code is used to negate the issue and ensure no data is lost
            retries = 10
            for attempt in range(retries):
                url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:{path}"
                headers = {'User-Agent': 'Mozilla/5.0'}  # Use a popular browser's User-Agent
                response = requests.get(url, headers=headers)

                if response.status_code == 200:
                    soup = BeautifulSoup(response.content, 'html.parser')
                    path_info = soup.get_text()

                    separator = "--------------------------------------------------"
                    index = path_info.find(separator)
                    new_text = path_info[index + len(separator):]

                    lines = [line for line in new_text.splitlines()]

                    extracted_genes = []

                    for line in lines:
                        start_index_path = line.find(":") + 1
                        end_index_path = line.find('  ', start_index_path)
                        if 0 <= start_index_path < end_index_path:
                            extracted_gene = line[start_index_path:end_index_path]
                            extracted_gene = extracted_gene.strip()
                            extracted_genes.append(extracted_gene)

                    # Write pathway, functions, and genes to the CSV file
                    csv_writer.writerow([path, function, ", ".join(extracted_genes)])
                    print(f"Finished: {path}")
                    break
                elif response.status_code == 403:
                    if attempt < retries - 1:
                        # Retry with exponential backoff: wait longer after each attempt
                        print(f"retrying: {url}")
                        wait_time = 30 * attempt + random.random()  # Add random jitter
                        time.sleep(wait_time)
                        print(f"Wait time: {wait_time}")
                        continue
                    else:
                        print(f"Request for {path} failed after {retries} attempts.")
                        error_links.append(path)
                else:
                    error_links.append(path)

        if error_links:
            print("Some links encountered errors while fetching data:", error_links)

    return output_file_path
